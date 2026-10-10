#!/usr/bin/env python3
"""Count exact Q20 junction words in primary, nonduplicate regional BAM reads.

This is an exploratory positive-evidence audit. Regional read slices cannot
establish absence of an allele, and fragment names/signatures are not UMI counts.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import hashlib
from pathlib import Path
import sys

from audit_common import encode, read_tsv, rolling_kmers, write_tsv

EXCLUDE_FLAGS = 0x100 | 0x800 | 0x200 | 0x400


def fragment_signature(read):
    """Canonical mapped read/mate starts; both records of a pair yield one key."""
    positions = set()
    if not read.is_unmapped and read.reference_id >= 0 and read.reference_start >= 0:
        positions.add((read.reference_name, read.reference_start, int(read.is_reverse)))
    mate_mapped = (read.is_paired and not read.mate_is_unmapped and
                   read.next_reference_id >= 0 and read.next_reference_start >= 0)
    if mate_mapped:
        positions.add((read.next_reference_name, read.next_reference_start, int(read.mate_is_reverse)))
    complete = int(not read.is_unmapped and mate_mapped)
    return json.dumps(sorted(positions), separators=(",", ":")), complete


def support_summary(witnesses):
    fragments = {w["qname"] for w in witnesses}
    # An empty signature supplies no independent start evidence.
    starts = {w["fragment_start_signature"] for w in witnesses if w["fragment_start_signature"] != "[]"}
    high = [w for w in witnesses if int(w["mapq"]) >= 20]
    high_fragments = {w["qname"] for w in high}
    high_starts = {w["fragment_start_signature"] for w in high if w["fragment_start_signature"] != "[]"}
    return dict(distinct_fragments=len(fragments), distinct_fragment_start_signatures=len(starts),
                fragments_mapq20=len(high_fragments), fragment_start_signatures_mapq20=len(high_starts),
                exploratory_supported=int(len(fragments) >= 3 and len(starts) >= 2),
                exploratory_supported_mapq20=int(len(high_fragments) >= 3 and len(high_starts) >= 2),
                fragments_with_complete_start_signature=len({w["qname"] for w in witnesses if w["signature_complete"]}))


def run(args):
    import pysam
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    markers, events, cases = read_tsv(args.markers), read_tsv(args.events), read_tsv(args.case_status)
    source = json.loads(Path(args.read_manifest).read_text())
    records = source["records"] if isinstance(source, dict) else source
    sources = {(r["sample"], r["locus"]): r for r in records}
    by_case = defaultdict(list)
    for marker in markers:
        by_case[marker["case_id"]].append(marker)
    all_witnesses, marker_support, junction_support, event_support, exposures = [], [], [], [], []
    events_by_case = defaultdict(list)
    for event in events:
        events_by_case[event["case_id"]].append(event)
    for case in cases:
        key, case_id = (case["sample"], case["locus"]), case["case_id"]
        print(f"Read scan: {case_id}", flush=True)
        source_rec = sources.get(key, {})
        bam_path = source_rec.get("bam_path")
        exposure = dict(case_id=case_id, sample=key[0], locus=key[1], role=case["role"],
                        sequence_status=case["status"], read_status="pending", reason="",
                        bam_sha256=source_rec.get("bam_sha256", ""),
                        bam_path=bam_path or "", reads_seen=0, eligible_reads=0,
                        excluded_flag_reads=0, missing_quality_reads=0, reads_with_wanted_q20_word=0,
                        eligible_fragments=0, marker_rows=len(by_case[case_id]))
        exposures.append(exposure)
        wanted = defaultdict(list)
        for marker in by_case[case_id]:
            if int(marker["assay_eligible"]):
                wanted[encode(marker["canonical_31mer"])].append(marker)
        per_marker, seen_witness, eligible_names = defaultdict(list), set(), set()
        try:
            if (source_rec.get("status") != "verified" or not source_rec.get("acquisition_verified")
                    or source_rec.get("counts", {}).get("all_records", 0) <= 0):
                raise ValueError("Read acquisition is not independently verified with positive record count")
            if not bam_path or not Path(bam_path).is_file():
                raise FileNotFoundError(source_rec.get("error") or "No completed regional BAM slice")
            digest = hashlib.sha256()
            with open(bam_path, "rb") as raw:
                for block in iter(lambda: raw.read(4*1024*1024), b""):
                    digest.update(block)
            if not source_rec.get("bam_sha256") or digest.hexdigest() != source_rec["bam_sha256"]:
                raise ValueError("BAM SHA256 absent from or inconsistent with source manifest")
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                for read in bam.fetch(until_eof=True):
                    exposure["reads_seen"] += 1
                    if read.flag & EXCLUDE_FLAGS:
                        exposure["excluded_flag_reads"] += 1
                        continue
                    if read.query_sequence is None or read.query_qualities is None:
                        exposure["missing_quality_reads"] += 1
                        continue
                    exposure["eligible_reads"] += 1
                    eligible_names.add(read.query_name)
                    signature, complete = fragment_signature(read)
                    read_words = {}
                    for start, word in rolling_kmers(read.query_sequence.upper(), read.query_qualities):
                        if word in wanted:
                            read_words.setdefault(word, start)
                    exposure["reads_with_wanted_q20_word"] += bool(read_words)
                    for word, start in read_words.items():
                        for marker in wanted[word]:
                            witness_key = (marker["marker_id"], read.query_name, read.flag,
                                           read.reference_id, read.reference_start)
                            if witness_key in seen_witness:
                                continue
                            seen_witness.add(witness_key)
                            witness = dict(marker_id=marker["marker_id"], event_id=marker["event_id"],
                                case_id=case_id, kind=marker["kind"], junction=marker["junction"],
                                primary_eligible=int(marker["primary_eligible"]), qname=read.query_name,
                                flag=read.flag, mapq=read.mapping_quality,
                                reference=read.reference_name or "*", reference_start=read.reference_start,
                                is_reverse=int(read.is_reverse), fragment_start_signature=signature,
                                signature_complete=complete, word_read_start=start,
                                minimum_word_quality=min(read.query_qualities[start:start+31]))
                            per_marker[marker["marker_id"]].append(witness)
                            all_witnesses.append(witness)
            exposure["eligible_fragments"] = len(eligible_names)
            if exposure["reads_seen"] != source_rec["counts"]["all_records"]:
                raise ValueError("Read scan count differs from verified acquisition manifest")
            exposure["read_status"] = "scanned"
        except Exception as exc:
            exposure.update(read_status="unresolved", reason=f"{type(exc).__name__}: {exc}")
        grouped, group_markers = defaultdict(list), defaultdict(list)
        for marker in by_case[case_id]:
            witnesses = per_marker[marker["marker_id"]]
            marker_support.append(dict(marker_id=marker["marker_id"], event_id=marker["event_id"],
                case_id=case_id, kind=marker["kind"], junction=marker["junction"],
                primary_eligible=marker["primary_eligible"], assay_eligible=marker["assay_eligible"],
                read_status=exposure["read_status"] if int(marker["assay_eligible"]) else "marker_ineligible_not_assayed",
                **support_summary(witnesses)))
            group = (marker["event_id"], marker["kind"], marker["junction"])
            group_markers[group].append(marker)
            if int(marker["assay_eligible"]):
                grouped[group].extend(witnesses)
        group_metrics = {}
        for group, group_rows in group_markers.items():
            metrics = support_summary(grouped[group])
            group_metrics[group] = metrics
            junction_support.append(dict(event_id=group[0], case_id=case_id, kind=group[1], junction=group[2],
                raw_marker_count=len({m["canonical_31mer"] for m in group_rows}),
                snp_clean_marker_count=len({m["canonical_31mer"] for m in group_rows if int(m["snp_clean"])}),
                eligible_marker_count=len({m["canonical_31mer"] for m in group_rows if int(m["assay_eligible"])}),
                read_status=exposure["read_status"], **metrics))
        for event in events_by_case[case_id]:
            ntruth = 2 if event["op"] == "I" else 1
            nalternative = 1 if event["op"] == "I" else 2
            row = dict(event_id=event["event_id"], case_id=case_id, sample=key[0], locus=key[1],
                       role=case["role"], op=event["op"], length=event["length"],
                       isolated=event["isolated"], primary_event_eligible=event["primary_event_eligible"],
                       read_status=exposure["read_status"],
                       flank_left_status=event.get("flank_left_status", "unresolved"),
                       flank_right_status=event.get("flank_right_status", "unresolved"))
            for kind, number in (("truth", ntruth), ("alternative", nalternative)):
                available, supported, high = 0, 0, 0
                for j in range(1, number+1):
                    group = (event["event_id"], kind, str(j))
                    markable = any(int(m["assay_eligible"]) for m in group_markers.get(group, []))
                    available += markable
                    metrics = group_metrics.get(group, {})
                    supported += metrics.get("exploratory_supported", 0)
                    high += metrics.get("exploratory_supported_mapq20", 0)
                row.update({f"{kind}_required_junctions": number,
                            f"{kind}_markable_junctions": available,
                            f"{kind}_supported_junctions": supported,
                            f"{kind}_full_support": int(supported == number and exposure["read_status"] == "scanned"),
                            f"{kind}_partial_support": int(0 < supported < number and exposure["read_status"] == "scanned"),
                            f"{kind}_full_support_mapq20": int(high == number and exposure["read_status"] == "scanned")})
                row[f"{kind}_assay_status"] = (
                    "read_unresolved" if exposure["read_status"] != "scanned" else
                    "complete_event_unassayable" if available < number else
                    "both_or_all_junction_words_supported" if supported == number else
                    "partial_junction_word_support" if supported else "no_support_flag")
            for side in ("left", "right"):
                group = (event["event_id"], "flank", side)
                row[f"{side}_flank_available"] = int(bool(group_markers.get(group)))
                row[f"{side}_flank_fragments"] = group_metrics.get(group, {}).get("distinct_fragments", 0)
            event_support.append(row)
    for filename, rows, fields in (
        ("read_witnesses.tsv.gz", all_witnesses, None),
        ("marker_support.tsv.gz", marker_support, None),
        ("junction_support.tsv", junction_support, None),
        ("event_support.tsv", event_support, None),
        ("case_exposure.tsv", exposures, None),
    ):
        write_tsv(out / filename, rows, fields)
    primary = [e for e in event_support if int(e["primary_event_eligible"])]
    summary = dict(cases_attempted=len(cases), cases_read_scanned=sum(x["read_status"] == "scanned" for x in exposures),
        cases_read_unresolved=sum(x["read_status"] != "scanned" for x in exposures),
        primary_events=len(primary), primary_events_truth_supported=sum(x["truth_full_support"] for x in primary),
        primary_events_truth_supported_mapq20=sum(x["truth_full_support_mapq20"] for x in primary),
        primary_case_ids_truth_supported=sorted({x["case_id"] for x in primary if x["truth_full_support"]}),
        primary_case_ids_truth_supported_mapq20=sorted({x["case_id"] for x in primary if x["truth_full_support_mapq20"]}),
        read_records_scanned=sum(x["reads_seen"] for x in exposures), witness_rows=len(all_witnesses),
        pysam_version=pysam.__version__, htslib_version=pysam.__samtools_version__,
        limits=["Truth-selected case-control audit", "Regional slices cannot support negative allele claims",
                "Exact words demonstrate sequence evidence, not phase, dosage, structural-class absence, or genotyping accuracy"],
        argv=sys.argv)
    (out / "read_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("markers", "events", "case-status", "read-manifest", "output"):
        p.add_argument("--"+name, required=True)
    p.add_argument("--execute-frozen-protocol", action="store_true")
    args = p.parse_args()
    if not args.execute_frozen_protocol:
        p.error("Freeze the protocol and source checksums before supplying --execute-frozen-protocol")
    run(args)


if __name__ == "__main__":
    main()
