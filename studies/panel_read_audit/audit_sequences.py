#!/usr/bin/env python3
"""Sequence audit of a frozen post-hoc cohort, using archived pairwise PAFs.

No real input is consumed unless --execute-frozen-protocol is supplied. PAF
distances reproduce upstream gt_dist v0.17.3's edit/alignment-size convention
with lexical representative-ID ties. The archive is not an exact edit-distance oracle.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

from audit_common import (Paf, alignment_segments, canonical, clean_marker,
                          donor, encode, isolated, junction_markers, opendata,
                          read_aliases, read_fasta, read_tsv, rolling_kmers, write_tsv)


def load_distance_index(path, groups):
    """Only small numeric distances in RAM; chosen CIGARs are read in a second pass."""
    distances, pair_lines, lengths = defaultdict(dict), {}, {}
    with opendata(path) as inp:
        for line_no, line in enumerate(inp, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            q, t = f[0], f[5]
            shift = int(f[9] == "+")
            size, matches = int(f[10+shift]), int(f[9+shift])
            if not 0 <= matches <= size or size == 0:
                raise ValueError("Invalid pairwise numeric distance")
            for h, n in ((q, int(f[1])), (t, int(f[6]))):
                for alias in groups[h]:
                    lengths.setdefault(alias, n)
            for qa in groups[q]:
                for ta in groups[t]:
                    distances[qa][ta] = distances[ta][qa] = (size - matches, size)
            pair_lines[tuple(sorted((q, t)))] = line_no
    # Preserve the upstream order: diagonal groups are added after PAF records.
    for hap, length in lengths.items():
        for h1 in groups.get(hap, [hap]):
            for h2 in groups.get(hap, [hap]):
                distances[h1][h2] = distances[h2][h1] = (0, length)
    return distances, pair_lines


def closest(distances, truth, sample, reverse):
    choices, ties = [], []
    for hap in truth:
        candidates = [(available, edit, size) for available, (edit, size) in distances[hap].items()
                      if donor(available) != sample]
        if not candidates:
            raise ValueError(f"No non-donor available haplotype for {hap}")
        from fractions import Fraction
        minimum = min(Fraction(edit, size) for _, edit, size in candidates)
        equal = [(available, edit, size) for available, edit, size in candidates
                 if Fraction(edit, size) == minimum]
        chosen = min(equal, key=lambda x: (reverse[x[0]], x[0]))
        choices.append(chosen)
        ties.append(len({reverse[x[0]] for x in equal}))
    return choices, ties


def marker_rows(event, truth, alt, other_truth, called_seqs, segments):
    rows = []
    for kind, seq, axis in (("truth", truth, "q"), ("alternative", alt, "a")):
        boundaries = [event[f"{axis}start"]]
        if event[f"{axis}end"] != boundaries[0]:
            boundaries.append(event[f"{axis}end"])
        for jix, junction in enumerate(boundaries, 1):
            for start, word in junction_markers(seq, junction):
                key = f'{event["event_id"]}|{kind}|j{jix}|{start}'
                rows.append(dict(marker_id=key, event_id=event["event_id"], case_id=event["case_id"],
                    sample=event["sample"], locus=event["locus"], truth_haplotype=event["truth_haplotype"],
                    kind=kind, junction=jix, sequence_start=start, junction_position=junction,
                    canonical_31mer=word, snp_clean=int(clean_marker(start, start+31, event, segments, axis)),
                    present_other_truth=int(word in other_truth or canonical_reverse(word) in other_truth),
                    present_called_haplotypes=sum(word in s or canonical_reverse(word) in s for s in called_seqs),
                    absent_both_truth=int(not any(word in s or canonical_reverse(word) in s for s in (truth, other_truth))),
                    within_locus_non_donor_present="", other_loci_present="", panel_private="", primary_eligible=""))
    for side, start in (("left", event["qstart"]-62), ("right", event["qend"]+31)):
        if start < 0 or start+31 > len(truth):
            event[f"flank_{side}_status"] = "out_of_bounds"
            continue
        seq = truth[start:start+31]
        if set(seq) - set("ACGT"):
            event[f"flank_{side}_status"] = "ambiguous_bases"
            continue
        word = canonical(seq)
        shared = word in alt or canonical_reverse(word) in alt
        event[f"flank_{side}_status"] = "available" if shared else "not_shared_with_alternative"
        if shared:
            rows.append(dict(marker_id=f'{event["event_id"]}|flank|{side}|{start}', event_id=event["event_id"],
                case_id=event["case_id"], sample=event["sample"], locus=event["locus"],
                truth_haplotype=event["truth_haplotype"], kind="flank", junction=side,
                sequence_start=start, junction_position="", canonical_31mer=word,
                snp_clean=int(clean_marker(start, start+31, event, segments)),
                present_other_truth=int(word in other_truth or canonical_reverse(word) in other_truth),
                present_called_haplotypes=sum(word in s or canonical_reverse(word) in s for s in called_seqs),
                absent_both_truth=0, within_locus_non_donor_present="", other_loci_present="", panel_private=0,
                primary_eligible=1))
    return rows


def canonical_reverse(word):
    from audit_common import revcomp
    return revcomp(word)


def run(args):
    cohort = read_tsv(args.cases) + read_tsv(args.controls)
    for row in cohort:
        row["case_id"] = row.get("case_id") or row["control_id"]
    if len({r["case_id"] for r in cohort}) != len(cohort):
        raise ValueError("Cohort cases/controls must be unique")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    panels, groups_by_locus = {}, {}
    choices_rows, resolution_rows, statuses, events, markers, gap_counts, cigar_records = [], [], [], [], [], [], []
    for locus in sorted({r["locus"] for r in cohort}):
        print(f"Sequence audit: {locus}", flush=True)
        directory = Path(args.database) / "loci" / locus
        try:
            seqs = read_fasta(directory / "haplotypes.fa.gz")
            groups, reverse = read_aliases(directory / "discarded_haplotypes.txt", seqs)
            distances, pair_lines = load_distance_index(directory / "haplotypes.paf.gz", groups)
            panels[locus], groups_by_locus[locus] = seqs, groups
        except Exception as exc:
            for row in cohort:
                if row["locus"] == locus:
                    statuses.append(dict(case_id=row["case_id"], sample=row["sample"], locus=locus, role=row["role"], status="panel_unresolved", reason=str(exc)))
            continue
        pending, required = [], set()
        for row in cohort:
            if row["locus"] != locus:
                continue
            status = dict(case_id=row["case_id"], sample=row["sample"], locus=locus, role=row["role"], status="pending", reason="")
            statuses.append(status)
            try:
                truth = row["expected_truth_haplotype_ids"].split(",")
                if len(truth) != 2 or any(h not in reverse or h not in distances for h in truth):
                    raise ValueError("Expected diploid truth aliases absent from panel distance index")
                for h in truth:
                    resolution_rows.append(dict(case_id=row["case_id"], sample=row["sample"], locus=locus,
                        truth_haplotype=h, representative=reverse[h], sequence_length=len(seqs[reverse[h]]),
                        sequence_sha256=hashlib.sha256(seqs[reverse[h]].encode()).hexdigest(),
                        sequence_aliases=",".join(groups[reverse[h]]),
                        non_donor_aliases=",".join(x for x in groups[reverse[h]] if donor(x) != row["sample"])))
                chosen, tie_counts = closest(distances, truth, row["sample"], reverse)
                total_edit, total_size = sum(x[1] for x in chosen), sum(x[2] for x in chosen)
                numeric_match = abs(total_edit/total_size - float(row["loo_avail_div"])) <= 1e-8
                comparison_matches = True
                for setting in ("full", "loo"):
                    predicted = row[f"{setting}_genotype"].split(",")
                    possibilities = []
                    for order in (predicted, predicted[::-1]):
                        d = [distances[h][g] for h, g in zip(truth, order)]
                        numerator, denominator = sum(x[0] for x in d), sum(x[1] for x in d)
                        possibilities.append((numerator / denominator, numerator, denominator))
                    calc_div, calc_edit, calc_size = min(possibilities)
                    compatible = abs(calc_div - float(row[f"{setting}_div"])) <= 1e-8
                    comparison_matches &= compatible
                    status.update({f"{setting}_comparison_compatible": int(compatible),
                                   f"{setting}_recomputed_edit": calc_edit,
                                   f"{setting}_recomputed_size": calc_size,
                                   f"{setting}_recomputed_div": f"{calc_div:.12g}"})
                status.update(recomputed_avail_edit=total_edit, recomputed_avail_size=total_size,
                              recomputed_avail_div=f"{total_edit/total_size:.12g}",
                              archived_avail_div=row["loo_avail_div"], numeric_compatibility=int(numeric_match),
                              numerator_denominator_exact_match=int(total_edit == int(row["loo_avail_edit"]) and total_size == int(row["loo_avail_size"])))
                for hap, (available, edit, size), tie_count in zip(truth, chosen, tie_counts):
                    rep, alt_rep = reverse[hap], reverse[available]
                    line_number = 0 if rep == alt_rep else pair_lines[tuple(sorted((rep, alt_rep)))]
                    choices_rows.append(dict(case_id=row["case_id"], sample=row["sample"], locus=locus,
                        truth_haplotype=hap, truth_representative=rep, alternative_haplotype=available,
                        alternative_representative=alt_rep, edit=edit, alignment_size=size,
                        divergence=f"{edit/size:.12g}", paf_line=line_number,
                        archived_numeric_compatibility=int(numeric_match), tied_representative_count=tie_count))
                    if line_number:
                        required.add(line_number)
                if not numeric_match or not comparison_matches:
                    status.update(status="archive_incompatible", reason="Closest or full/LOO comparisons do not reproduce archived divergence; ties_present="+str(any(n > 1 for n in tie_counts)))
                    continue
                pending.append((row, status, truth, chosen))
            except Exception as exc:
                status.update(status="truth_or_distance_unresolved", reason=str(exc))
        selected_pafs = {}
        with opendata(directory / "haplotypes.paf.gz") as inp:
            for n, line in enumerate(inp, 1):
                if n in required:
                    try:
                        selected_pafs[n] = Paf.parse(line, n)
                    except Exception as exc:
                        selected_pafs[n] = exc
        for row, status, truth_haps, chosen in pending:
            try:
                case_fully_aligned = True
                for hix, (hap, (available, _, _)) in enumerate(zip(truth_haps, chosen), 1):
                    rep, alt_rep = reverse[hap], reverse[available]
                    if rep == alt_rep:
                        truth, alt = seqs[rep], seqs[alt_rep]
                        qstart = astart = 0
                        cigar, line_number, strand = [(len(truth), "=")], 0, "+"
                    else:
                        line_number = pair_lines[tuple(sorted((rep, alt_rep)))]
                        paf = selected_pafs[line_number]
                        if isinstance(paf, Exception):
                            raise paf
                        truth, alt, qstart, astart, cigar = paf.orient(rep, seqs)
                        strand = paf.strand
                    segments = alignment_segments(cigar, qstart, astart)
                    qend = segments[-1]["qend"]
                    aend = segments[-1]["aend"]
                    fully_aligned = qstart == 0 and astart == 0 and qend == len(truth) and aend == len(alt)
                    case_fully_aligned &= fully_aligned
                    cigar_records.append(dict(case_id=row["case_id"], locus=locus, truth_haplotype=hap,
                        alternative_haplotype=available, paf_line=line_number, original_strand=strand,
                        truth_start=qstart, truth_end=qend, truth_length=len(truth), alternative_start=astart,
                        alternative_end=aend, alternative_length=len(alt), fully_aligned=int(fully_aligned),
                        normalized_cigar="".join(f"{n}{op}" for n, op in cigar)))
                    gaps = [s for s in segments if s["op"] in "ID"]
                    gap_counts.append(dict(case_id=row["case_id"], locus=locus, truth_haplotype=hap,
                        gaps_lt50=sum(s["length"] < 50 for s in gaps),
                        gaps_50_1000=sum(50 <= s["length"] <= 1000 for s in gaps),
                        gaps_gt1000=sum(s["length"] > 1000 for s in gaps),
                        mismatch_bases=sum(s["length"] for s in segments if s["op"] == "X"),
                        **{f"{name}_{label}_bases": sum(s["length"] for s in gaps if s["op"] == op and low <= s["length"] <= high)
                           for name, op in (("inserted", "I"), ("deleted", "D"))
                           for label, low, high in (("lt50", 1, 49), ("50_1000", 50, 1000), ("gt1000", 1001, float("inf")))},
                        fully_aligned=int(fully_aligned), unaligned_truth_bases=qstart+len(truth)-qend,
                        unaligned_alternative_bases=astart+len(alt)-aend))
                    other_truth = seqs[reverse[truth_haps[1 if hix == 1 else 0]]]
                    called_seqs = [seqs[reverse[h]] for h in row["loo_genotype"].split(",")]
                    for gap in gaps:
                        if not 50 <= gap["length"] <= 1000:
                            continue
                        event = dict(event_id=f'{row["case_id"]}|h{hix}|c{gap["cigar_index"]}',
                            case_id=row["case_id"], sample=row["sample"], locus=locus, role=row["role"],
                            truth_haplotype=hap, alternative_haplotype=available, paf_line=line_number,
                            **gap, isolated=int(isolated(gap, segments)), fully_aligned=int(fully_aligned),
                            local_mismatch_bases=sum(max(0, min(s["qend"], gap["qend"]+31)-max(s["qstart"], gap["qstart"]-31))
                                for s in segments if s["op"] == "X"))
                        event["primary_event_eligible"] = int(event["isolated"] and fully_aligned)
                        events.append(event)
                        markers.extend(marker_rows(event, truth, alt, other_truth, called_seqs, segments))
                status["status"] = "sequence_resolved" if case_fully_aligned else "partial_alignment_unresolved"
                if not case_fully_aligned:
                    status["reason"] = "At least one selected comparison does not align both complete sequences"
            except Exception as exc:
                status.update(status="alignment_unresolved", reason=str(exc))
                # Avoid presenting partially processed diploid case as eligible.
                for event in events:
                    if event["case_id"] == row["case_id"]:
                        event["primary_event_eligible"] = 0
        # Numeric index is large; release it before moving to next locus.
        del distances
    # Screen only candidate words rather than retaining whole-genome-sized kmer sets.
    wanted = {encode(r["canonical_31mer"]) for r in markers}
    presence = defaultdict(lambda: defaultdict(set))
    for locus, seqs in panels.items():
        print(f"Marker background screen: {locus}", flush=True)
        for rep, seq in seqs.items():
            found = {word for _, word in rolling_kmers(seq) if word in wanted}
            for word in found:
                presence[word][locus].add(rep)
    event_by_id = {r["event_id"]: r for r in events}
    background_complete = set(panels) == {r["locus"] for r in cohort}
    for marker in markers:
        loc, sample = marker["locus"], marker["sample"]
        hits = presence[encode(marker["canonical_31mer"])]
        non_donor = sorted({alias for rep in hits.get(loc, ()) for alias in groups_by_locus[loc][rep]
                            if donor(alias) != sample})
        other_loci = sorted(set(hits) - {loc})
        marker["within_locus_non_donor_present"] = ",".join(non_donor)
        marker["other_loci_present"] = ",".join(other_loci)
        marker["panel_private"] = int(not non_donor)
        event = event_by_id[marker["event_id"]]
        if marker["kind"] == "truth":
            marker["assay_eligible"] = int(not non_donor and not other_loci and marker["snp_clean"] and event["fully_aligned"] and background_complete)
            marker["primary_eligible"] = int(marker["assay_eligible"] and event["primary_event_eligible"])
            if not non_donor and marker["present_called_haplotypes"]:
                raise ValueError("Panel-private marker occurs in called LOO haplotype; alias/exclusion inconsistency")
        elif marker["kind"] == "alternative":
            marker["assay_eligible"] = int(marker["absent_both_truth"] and not other_loci and marker["snp_clean"] and event["fully_aligned"] and background_complete)
            marker["primary_eligible"] = int(marker["assay_eligible"] and event["primary_event_eligible"])
        else:
            marker["assay_eligible"] = 1
        marker["background_complete"] = int(background_complete)
    for filename, rows in (("case_status.tsv", statuses), ("resolved_haplotypes.tsv", resolution_rows),
                           ("pairwise_choices.tsv", choices_rows), ("events.tsv", events),
                           ("markers.tsv.gz", markers), ("gap_counts.tsv", gap_counts),
                           ("normalized_alignments.tsv.gz", cigar_records)):
        write_tsv(out / filename, rows)
    summary = dict(cases_attempted=len(cohort), sequence_resolved=sum(s["status"] == "sequence_resolved" for s in statuses),
                   statuses={k: sum(s["status"] == k for s in statuses) for k in sorted({s["status"] for s in statuses})},
                   events_50_1000=len(events), primary_events=sum(e["primary_event_eligible"] for e in events),
                   raw_marker_rows=len(markers), primary_truth_markers=sum(m["kind"] == "truth" and m["primary_eligible"] for m in markers),
                   background_loci=sorted(panels), background_scope="Selected 15 locus panels; not whole genome",
                   coordinate_convention="0-based; truth-forward; alternative oriented to truth", argv=sys.argv)
    (out / "sequence_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--database", required=True)
    p.add_argument("--cases", required=True)
    p.add_argument("--controls", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--execute-frozen-protocol", action="store_true")
    args = p.parse_args()
    if not args.execute_frozen_protocol:
        p.error("Freeze the protocol and source checksums before supplying --execute-frozen-protocol")
    run(args)


if __name__ == "__main__":
    main()
