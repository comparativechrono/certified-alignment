"""Toy checks cover orientation, marker specificity, aliases and fragment counting."""
import gzip
from pathlib import Path
import random
import tempfile
import types
import unittest
import contextlib
import io
import json
import hashlib

from audit_common import (Paf, alignment_segments, canonical, clean_marker, encode,
                          isolated, junction_markers, read_aliases, revcomp, rolling_kmers,
                          read_tsv, write_tsv)
from audit_sequences import closest, load_distance_index
from scan_reads import fragment_signature, support_summary


def dna(n, seed=1):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


class AuditTests(unittest.TestCase):
    def test_kmer_strand_quality_and_ambiguity(self):
        seq = dna(70)
        expected = [(i, encode(canonical(seq[i:i+31]))) for i in range(40)]
        self.assertEqual(list(rolling_kmers(seq)), expected)
        self.assertEqual({x[1] for x in rolling_kmers(seq)}, {x[1] for x in rolling_kmers(revcomp(seq))})
        q = [30]*70
        q[35] = 19
        self.assertEqual([i for i, _ in rolling_kmers(seq, q)], list(range(5))+list(range(36, 40)))
        self.assertEqual([i for i, _ in rolling_kmers(seq[:35]+"N"+seq[36:])], list(range(5))+list(range(36, 40)))

    def test_paf_orientation_swap_and_minus(self):
        left, insertion, right = dna(30, 1), dna(50, 2), dna(50, 3)
        truth, alt = left+insertion+right, left+right
        cases = [
            (f"Q\t130\t0\t130\t+\tA\t80\t0\t80\t80\t130\t60\tcg:Z:30=50I50=", {"Q":truth, "A":alt}),
            (f"A\t80\t0\t80\t+\tQ\t130\t0\t130\t80\t130\t60\tcg:Z:30=50D50=", {"Q":truth, "A":alt}),
            (f"Q\t130\t0\t130\t-\tA\t80\t0\t80\t80\t130\t60\tcg:Z:50=50I30=", {"Q":truth, "A":revcomp(alt)}),
            (f"A\t80\t0\t80\t-\tQ\t130\t0\t130\t80\t130\t60\tcg:Z:30=50D50=", {"Q":truth, "A":revcomp(alt)}),
        ]
        for line, seqs in cases:
            parsed = Paf.parse(line)
            q, a, qs, ast, cigar = parsed.orient("Q", seqs)
            self.assertEqual((q, a, qs, ast, cigar), (truth, alt, 0, 0, [(30,"="),(50,"I"),(50,"=")]))

    def test_paf_legacy_and_invalid_consumption(self):
        p = Paf.parse("Q\t40\t0\t40\t+\tA\t40\t0\t40\t+\t40\t40\t60\tcg:Z:40=")
        self.assertTrue(p.legacy_shift)
        with self.assertRaisesRegex(ValueError, "consumption"):
            Paf.parse("Q\t40\t0\t40\t+\tA\t40\t0\t40\t40\t40\t60\tcg:Z:39=")
        with self.assertRaisesRegex(ValueError, "actual sequence"):
            p.orient("Q", {"Q":"A"*40, "A":"C"*40})

    def test_verified_m(self):
        seq = dna(40)
        p = Paf.parse("Q\t40\t0\t40\t+\tA\t40\t0\t40\t40\t40\t60\tcg:Z:40M")
        self.assertEqual(p.orient("Q", {"Q":seq,"A":seq})[-1], [(40,"=")])

    def test_isolation_includes_small_gaps_and_boundary(self):
        seg = alignment_segments([(50,"="),(1,"D"),(31,"="),(50,"I"),(100,"=")])
        self.assertFalse(isolated(seg[3], seg))
        seg = alignment_segments([(50,"="),(1,"D"),(32,"="),(50,"I"),(100,"=")])
        self.assertTrue(isolated(seg[3], seg))

    def test_snp_clean_guard(self):
        seg = alignment_segments([(95,"="),(1,"X"),(4,"="),(50,"I"),(100,"=")])
        self.assertFalse(clean_marker(90,121,seg[3],seg))
        self.assertTrue(clean_marker(96,127,seg[3],seg))
        # Deletion junction has no query bases in the focal gap.
        seg = alignment_segments([(100,"="),(50,"D"),(100,"=")])
        self.assertTrue(clean_marker(80,111,seg[1],seg))
        seg = alignment_segments([(90,"="),(1,"D"),(10,"="),(50,"I"),(100,"=")])
        self.assertFalse(clean_marker(85,116,seg[3],seg))
        self.assertEqual(len(list(junction_markers(dna(200),100))),12)

    def test_alias_non_donor_and_lexical_ties(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"alias.txt"
            p.write_text("D.1 = Z.1\n")
            groups, reverse = read_aliases(p,{"D.1":"A"*100,"B.1":"A"*100,"D.2":"A"*100})
            self.assertEqual(reverse["Z.1"],"D.1")
            d={"D.1":{"D.1":(0,100),"D.2":(1,100),"Z.1":(0,100),"B.1":(0,100)},
               "D.2":{"D.2":(0,100),"Z.1":(2,100),"B.1":(1,50)}}
            choices,ties=closest(d,["D.1","D.2"],"D",reverse)
            self.assertEqual([x[0] for x in choices],["B.1","B.1"])
            self.assertEqual(ties,[2,2])

    def test_fragment_mates_collapse_starts(self):
        r = types.SimpleNamespace(is_unmapped=False,reference_id=0,reference_start=100,
            reference_name="chr1",is_reverse=False,is_paired=True,mate_is_unmapped=False,
            next_reference_id=0,next_reference_start=300,next_reference_name="chr1",mate_is_reverse=True)
        mate = types.SimpleNamespace(**vars(r))
        mate.reference_start,mate.next_reference_start=300,100
        mate.is_reverse,mate.mate_is_reverse=True,False
        self.assertEqual(fragment_signature(r),fragment_signature(mate))
        sig,complete=fragment_signature(r)
        witnesses=[dict(qname=f"r{i}",fragment_start_signature=sig,signature_complete=complete,mapq=60) for i in range(3)]
        self.assertEqual(support_summary(witnesses)["exploratory_supported"],0)
        r.reference_start=101
        witnesses[-1]["fragment_start_signature"]=fragment_signature(r)[0]
        self.assertEqual(support_summary(witnesses)["exploratory_supported"],1)
        self.assertEqual(support_summary(witnesses*2)["distinct_fragments"],3)

    def test_toy_sequence_pipeline_and_read_support(self):
        from audit_sequences import run as sequence_run
        from scan_reads import run as read_run
        try:
            import pysam
        except ImportError:
            self.skipTest("pysam not installed; rerun with the pinned pysam PYTHONPATH")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            locus = root/"db"/"loci"/"L"
            locus.mkdir(parents=True)
            left, ins, right = dna(80,10), dna(50,11), dna(80,12)
            truth, alt = left+ins+right, left+right
            with gzip.open(locus/"haplotypes.fa.gz","wt") as out:
                out.write(f">D.1\n{truth}\n>B.1\n{alt}\n")
            (locus/"discarded_haplotypes.txt").write_text("B.1 = D.2\n")
            with gzip.open(locus/"haplotypes.paf.gz","wt") as out:
                out.write("D.1\t210\t0\t210\t+\tB.1\t160\t0\t160\t160\t210\t60\tcg:Z:80=50I80=\n")
            row=dict(case_id="D__L",sample="D",locus="L",role="toy",expected_truth_haplotype_ids="D.1,D.2",
                     loo_avail_edit=50,loo_avail_size=370,loo_avail_div=f"{50/370:.9f}",
                     loo_genotype="B.1,B.1",loo_div=f"{50/370:.9f}",full_genotype="D.1,D.2",full_div="0")
            write_tsv(root/"cases.tsv",[row])
            write_tsv(root/"controls.tsv",[],list(row))
            output=root/"results"
            with contextlib.redirect_stdout(io.StringIO()):
                sequence_run(types.SimpleNamespace(cases=root/"cases.tsv",controls=root/"controls.tsv",
                              database=root/"db",output=output))
            self.assertEqual(read_tsv(output/"case_status.tsv")[0]["status"],"sequence_resolved")
            event=read_tsv(output/"events.tsv")[0]
            self.assertEqual(event["isolated"],"1")
            self.assertEqual(event["flank_left_status"],"available")
            bam=root/"toy.bam"
            header={"HD":{"VN":"1.6"},"SQ":[{"SN":"chr1","LN":10000}]}
            with pysam.AlignmentFile(bam,"wb",header=header) as out:
                for i in range(3):
                    read=pysam.AlignedSegment(out.header)
                    read.query_name=f"fragment{i}"
                    read.query_sequence=truth
                    read.query_qualities=pysam.qualitystring_to_array("I"*len(truth))
                    read.flag=65
                    read.reference_id=0
                    read.reference_start=100+i
                    read.mapping_quality=60
                    read.cigarstring=f"{len(truth)}M"
                    read.next_reference_id=0
                    read.next_reference_start=400+i
                    out.write(read)
            (root/"reads.json").write_text(json.dumps({"records":[{"sample":"D","locus":"L","bam_path":str(bam),
                "bam_sha256":hashlib.sha256(bam.read_bytes()).hexdigest(), "status":"verified",
                "acquisition_verified":True,"counts":{"all_records":3}}]}))
            with contextlib.redirect_stdout(io.StringIO()):
                read_run(types.SimpleNamespace(markers=output/"markers.tsv.gz",events=output/"events.tsv",
                           case_status=output/"case_status.tsv",read_manifest=root/"reads.json",output=output))
            supported=read_tsv(output/"event_support.tsv")[0]
            self.assertEqual(supported["truth_full_support"],"1")
            self.assertEqual(supported["alternative_full_support"],"0")
            self.assertEqual(supported["truth_supported_junctions"],"2")


if __name__ == "__main__":
    unittest.main()
