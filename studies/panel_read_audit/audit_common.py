"""Auditable sequence primitives; no caller-confidence or phenotype inference."""
from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path

K = 31
BASES = "ACGT"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
ENC = {ord(c): i for i, c in enumerate(BASES)}
MASK = (1 << (2 * K)) - 1


def opendata(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def canonical(seq):
    return min(seq, revcomp(seq))


def encode(seq):
    x = 0
    for c in seq.encode():
        x = (x << 2) | ENC[c]
    return x


def rolling_kmers(seq, qualities=None, minimum_quality=20):
    """Yield (zero-based read start, canonical 2-bit 31mer), reset on N/low Q."""
    fwd = rev = valid = 0
    for i, c in enumerate(seq.encode()):
        b = ENC.get(c)
        if b is None or (qualities is not None and qualities[i] < minimum_quality):
            fwd = rev = valid = 0
            continue
        fwd = ((fwd << 2) | b) & MASK
        rev = (rev >> 2) | ((3 - b) << (2 * (K - 1)))
        valid += 1
        if valid >= K:
            yield i - K + 1, min(fwd, rev)


def read_fasta(path):
    records, name, seq = {}, None, []
    with opendata(path) as inp:
        for line in inp:
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise ValueError(f"Duplicate FASTA identifier {name}")
                    records[name] = "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        if name in records:
            raise ValueError(f"Duplicate FASTA identifier {name}")
        records[name] = "".join(seq).upper()
    return records


def read_aliases(path, seqs):
    groups = {hap: [hap] for hap in seqs}
    if Path(path).exists():
        with opendata(path) as inp:
            for line in inp:
                if not line.strip() or line.startswith("#"):
                    continue
                rep, aliases = line.strip().split("=")
                rep = rep.strip()
                if rep not in seqs:
                    raise ValueError(f"Alias representative absent from FASTA: {rep}")
                groups[rep].extend(x.strip() for x in aliases.split(","))
    reverse = {}
    for rep, aliases in groups.items():
        for alias in aliases:
            if alias in reverse and reverse[alias] != rep:
                raise ValueError(f"Conflicting alias {alias}")
            reverse[alias] = rep
    return groups, reverse


def donor(hap):
    return re.sub(r"[._][1-9]$", "", hap)


def parse_cigar(cigar):
    parts = [(int(n), op) for n, op in re.findall(r"(\d+)([=XMID])", cigar)]
    if not parts or "".join(f"{n}{op}" for n, op in parts) != cigar or any(n < 1 for n, _ in parts):
        raise ValueError(f"Unsupported/malformed CIGAR: {cigar[:100]}")
    return parts


@dataclass
class Paf:
    query: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    target: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    size: int
    cigar: list
    line_number: int = 0
    legacy_shift: bool = False

    @classmethod
    def parse(cls, line, line_number=0):
        f = line.rstrip().split("\t")
        shift = int(f[9] == "+")
        tags = {x[:2]: x[5:] for x in f[12 + shift:] if len(x) > 5 and x[2] == ":"}
        obj = cls(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5],
                  int(f[6]), int(f[7]), int(f[8]), int(f[9 + shift]),
                  int(f[10 + shift]), parse_cigar(tags["cg"]), line_number, bool(shift))
        if obj.strand not in ("+", "-"):
            raise ValueError(f"Invalid strand {obj.strand}")
        if not (0 <= obj.qstart <= obj.qend <= obj.qlen and 0 <= obj.tstart <= obj.tend <= obj.tlen):
            raise ValueError("PAF coordinates out of bounds")
        qcons = sum(n for n, op in obj.cigar if op in "=XMI")
        tcons = sum(n for n, op in obj.cigar if op in "=XMD")
        if qcons != obj.qend - obj.qstart or tcons != obj.tend - obj.tstart:
            raise ValueError("PAF CIGAR/query/target consumption mismatch")
        if sum(n for n, _ in obj.cigar) != obj.size:
            raise ValueError("PAF alignment size differs from CIGAR columns")
        if not 0 <= obj.matches <= obj.size:
            raise ValueError("PAF matches invalid")
        return obj

    @property
    def edit(self):
        return self.size - self.matches

    def orient(self, truth_rep, seqs):
        """Truth always forward; alternative may reverse-complement. Coordinates are
        zero-based in these oriented sequences. Minus-strand PAF CIGAR proceeds
        target-forward, as defined by PAF; validate actual bases below.
        """
        if len(seqs[self.query]) != self.qlen or len(seqs[self.target]) != self.tlen:
            raise ValueError("FASTA/PAF length mismatch")
        if truth_rep == self.query:
            if self.strand == "+":
                cigar, qstart, astart = self.cigar, self.qstart, self.tstart
                alt = seqs[self.target]
            else:
                cigar, qstart, astart = self.cigar[::-1], self.qstart, self.tlen - self.tend
                alt = revcomp(seqs[self.target])
        elif truth_rep == self.target:
            cigar = [(n, {"I": "D", "D": "I"}.get(op, op)) for n, op in self.cigar]
            qstart, astart = self.tstart, self.qstart if self.strand == "+" else self.qlen - self.qend
            alt = seqs[self.query] if self.strand == "+" else revcomp(seqs[self.query])
        else:
            raise ValueError("Truth is neither PAF endpoint")
        truth = seqs[truth_rep]
        validate_bases(truth, alt, qstart, astart, cigar, self.matches)
        # Convert M to verified =/X runs for the SNP-clean marker guard.
        expanded, q, a = [], qstart, astart
        for n, op in cigar:
            if op == "M":
                for x, y in zip(truth[q:q+n], alt[a:a+n]):
                    curr = "=" if x == y else "X"
                    if expanded and expanded[-1][1] == curr:
                        expanded[-1] = (expanded[-1][0] + 1, curr)
                    else:
                        expanded.append((1, curr))
            else:
                expanded.append((n, op))
            q += n if op in "=XMI" else 0
            a += n if op in "=XMD" else 0
        cigar = expanded
        return truth, alt, qstart, astart, cigar


def validate_bases(truth, alt, qstart, astart, cigar, expected_matches=None):
    q, a, matches = qstart, astart, 0
    for n, op in cigar:
        if op in "=XM":
            left, right = truth[q:q+n], alt[a:a+n]
            if len(left) != n or len(right) != n:
                raise ValueError("Alignment exceeds sequence")
            equal = sum(x == y for x, y in zip(left, right))
            if (op == "=" and equal != n) or (op == "X" and equal):
                raise ValueError(f"CIGAR {op} disagrees with actual sequence at {q}/{a}")
            matches += equal
        if op in "=XMI":
            q += n
        if op in "=XMD":
            a += n
    if q > len(truth) or a > len(alt):
        raise ValueError("Alignment consumes beyond sequence")
    if expected_matches is not None and matches != expected_matches:
        raise ValueError("PAF match count disagrees with actual sequence")


def alignment_segments(cigar, qstart=0, astart=0):
    q, a = qstart, astart
    out = []
    for ix, (n, op) in enumerate(cigar):
        qend, aend = q + (n if op in "=XMI" else 0), a + (n if op in "=XMD" else 0)
        out.append(dict(cigar_index=ix, length=n, op=op, qstart=q, qend=qend, astart=a, aend=aend))
        q, a = qend, aend
    return out


def interval_distance(x, left, right):
    return max(left - x, x - right, 0)


def isolated(event, segments):
    return not any(s["op"] in "ID" and s["cigar_index"] != event["cigar_index"]
                   and min(interval_distance(event[b], s["qstart"], s["qend"]) for b in ("qstart", "qend")) <= K
                   for s in segments)


def clean_marker(start, end, event, segments, axis="q"):
    """Every outside-event base must come from exact-match (=) columns."""
    covered = 0
    for s in segments:
        if (s["op"] in "ID" and s["cigar_index"] != event["cigar_index"]
                and s[f"{axis}start"] == s[f"{axis}end"]
                and start < s[f"{axis}start"] < end):
            return False
        a, b = max(start, s[f"{axis}start"]), min(end, s[f"{axis}end"])
        if a < b:
            if s["op"] == "=" or s["cigar_index"] == event["cigar_index"]:
                covered += b - a
            else:
                return False
    return covered == end - start


def junction_markers(sequence, junction):
    # Both sides contribute at least ten bases: starts j-21 through j-10.
    for start in range(junction - 21, junction - 9):
        if start < 0 or start + K > len(sequence):
            continue
        word = sequence[start:start+K]
        if set(word) <= set(BASES):
            yield start, canonical(word)


def read_tsv(path):
    with opendata(path) as inp:
        return list(csv.DictReader(inp, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with opendata(path, "wt") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
