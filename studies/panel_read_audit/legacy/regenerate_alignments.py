#!/usr/bin/env python3
"""Regenerate all pairwise PAFs with the frozen release-matched runtime.

Run only after the protocol amendment is frozen. All FASTA inputs must match
the panel source manifest before any alignment starts. No outcome-dependent
parameters or pair selection occur. Existing PAFs are never overwritten.
"""

import argparse
import datetime
import gzip
import json
import pathlib
import resource
import shutil
import subprocess
import time

from prepare_legacy_aligner import ALIGN_FLAGS, BINARY_SHA256, sha256


def fasta_headers(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        names = [line[1:].split()[0] for line in handle if line.startswith('>')]
    assert len(names) == len(set(names)), f'Duplicate FASTA names: {path}'
    assert len(names) > 1, f'Insufficient sequences: {path}'
    return names


def verify_pairs(path, names):
    """Every unordered pair must have one complete PAF record and CIGAR."""
    expected = {frozenset((a, b)) for i, a in enumerate(names) for b in names[:i]}
    observed = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            assert len(fields) >= 12, 'Short PAF record'
            pair = frozenset((fields[0], fields[5]))
            assert pair in expected and pair not in observed, 'Unknown or duplicate PAF pair'
            assert any(tag.startswith('cg:Z:') for tag in fields[12:]), 'Missing CIGAR'
            observed.add(pair)
    assert observed == expected, 'Incomplete all-pairs PAF'
    return len(observed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=pathlib.Path)
    parser.add_argument('--database', required=True, type=pathlib.Path,
                        help='Extracted db directory containing loci/')
    parser.add_argument('--panel-manifest', required=True, type=pathlib.Path)
    parser.add_argument('--output-database', required=True, type=pathlib.Path,
                        help='PAFs are written under this directory/loci/LOCUS/')
    parser.add_argument('--manifest-output', required=True, type=pathlib.Path)
    parser.add_argument('--timeout-seconds', type=int, default=1800)
    parser.add_argument('--max-address-space-gb', type=int, default=16,
                        help='Per subprocess virtual-memory ceiling (GiB)')
    args = parser.parse_args()
    assert args.timeout_seconds > 0 and args.max_address_space_gb > 0
    binary = args.binary.resolve()
    if sha256(binary) != BINARY_SHA256:
        raise ValueError('Executable checksum mismatch')
    source = json.loads(args.panel_manifest.read_text())
    loci = sorted(source['selected_loci'])
    entries = {row['archive_member']: row for row in source['files']}
    planned = []
    for locus in loci:
        assert pathlib.Path(locus).name == locus
        source_entry = entries[f'db/loci/{locus}/haplotypes.fa.gz']
        fasta = args.database.resolve() / 'loci' / locus / 'haplotypes.fa.gz'
        assert fasta.stat().st_size == source_entry['bytes'], f'FASTA size mismatch: {locus}'
        assert sha256(fasta) == source_entry['sha256'], f'FASTA checksum mismatch: {locus}'
        names = fasta_headers(fasta)
        target = args.output_database.resolve() / 'loci' / locus / 'haplotypes.paf.gz'
        assert not target.exists(), f'Output already exists; refusing overwrite: {target}'
        planned.append((locus, fasta, target, names, source_entry))
    max_bytes = args.max_address_space_gb * 1024**3

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))

    report = {
        'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'binary_sha256': BINARY_SHA256,
        'panel_manifest_sha256': sha256(args.panel_manifest),
        'alignment_flags': ALIGN_FLAGS,
        'timeout_seconds_per_locus': args.timeout_seconds,
        'max_address_space_bytes_per_locus': max_bytes,
        'interpretation': 'Regenerated release-matched alignments; archived outcome compatibility remains to be checked.',
        'loci': [],
    }
    args.manifest_output.parent.mkdir(parents=True, exist_ok=True)

    def save_report():
        args.manifest_output.write_text(json.dumps(report, indent=2) + '\n')

    save_report()
    for locus, fasta, target, names, source_entry in planned:
        target.parent.mkdir(parents=True, exist_ok=True)
        raw = target.with_name('regenerating.paf')
        log_path = target.with_name('regeneration.log')
        command = [str(binary), 'align', '-i', str(fasta), '-o', str(raw), '-A'] + ALIGN_FLAGS
        row = {'locus': locus, 'input_sha256': source_entry['sha256'],
               'sequence_count': len(names), 'expected_pairs': len(names)*(len(names)-1)//2,
               'command': command, 'status': 'running'}
        report['loci'].append(row)
        save_report()
        started = time.monotonic()
        try:
            with log_path.open('w') as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                        timeout=args.timeout_seconds, preexec_fn=limits)
            row['returncode'] = result.returncode
            if result.returncode:
                row['status'] = 'runtime_failed'
            else:
                row['observed_pairs'] = verify_pairs(raw, names)
                partial = target.with_name('regenerating.paf.gz')
                with raw.open('rb') as src, partial.open('wb') as out:
                    with gzip.GzipFile(filename='', mode='wb', fileobj=out, mtime=0) as compressed:
                        shutil.copyfileobj(src, compressed)
                row['uncompressed_sha256'] = sha256(raw)
                row['output_sha256'] = sha256(partial)
                row['output_bytes'] = partial.stat().st_size
                partial.replace(target)
                raw.unlink()
                row['status'] = 'complete_pending_outcome_compatibility'
        except subprocess.TimeoutExpired:
            row['status'] = 'timeout'
        except (ValueError, AssertionError, OSError) as error:
            row['status'] = 'validation_failed'
            row['error'] = str(error)
        row['elapsed_seconds'] = round(time.monotonic() - started, 3)
        row['log_sha256'] = sha256(log_path)
        save_report()
        print(json.dumps({'locus': locus, 'status': row['status'],
                          'elapsed_seconds': row['elapsed_seconds']}), flush=True)
    report['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    report['all_complete'] = all(row['status'] == 'complete_pending_outcome_compatibility'
                                 for row in report['loci'])
    save_report()
    if not report['all_complete']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
