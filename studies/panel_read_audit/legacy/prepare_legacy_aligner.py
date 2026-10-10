#!/usr/bin/env python3
"""Verify a release-matched Locityper 0.17.3 runtime; align synthetic DNA only.

This downloads a 2.90 MB Bioconda package and extracts it without installing
Conda or changing the system. Requires Python >=3.12 and Linux x86_64.
It does not reconstruct the authors' original benchmark build: the packaged
recipe did not pin the WFA2 dependency's commit.
"""

import argparse
import hashlib
import json
import pathlib
import random
import re
import subprocess
import tarfile
import urllib.request

PACKAGE_URL = (
    'https://api.anaconda.org/download/bioconda/locityper/0.17.3/'
    'linux-64/locityper-0.17.3-h919a2d8_0.tar.bz2'
)
PACKAGE_SIZE = 2904704
PACKAGE_SHA256 = '8a9d28e7ba8a9605fbab6a18d370dc41948822365eadbd64413a7aadfff3db0a'
BINARY_SHA256 = 'cc73d9e285584e1969f68637ce95f789e77c985fbd76c9c818cce632fc7b371d'
SOURCE_COMMIT = '146cfd42e9179bafd63d1f41db7fc449940989de'
ALIGN_FLAGS = ['-m', '15', '15', '-D', '1', '-k', '25,51,101',
               '-g', '500', '-a', '9', '-M', '4', '-O', '6', '-E', '1', '-@', '8']


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify_paf(path, sequences):
    """Validate every synthetic alignment against its actual sequence bases."""
    records = []
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        tags = dict(item.split(':', 2)[::2] for item in fields[12:])
        query, reference = sequences[fields[0]], sequences[fields[5]]
        assert fields[4] == '+'
        assert list(map(int, fields[1:4])) == [len(query), 0, len(query)]
        assert list(map(int, fields[6:9])) == [len(reference), 0, len(reference)]
        assert 'cg' in tags, 'Synthetic pair was unexpectedly skipped'
        ops = re.findall(r'(\d+)([=XID])', tags['cg'])
        assert ''.join(n + op for n, op in ops) == tags['cg']
        qpos = rpos = matches = errors = score = 0
        for ntext, op in ops:
            n = int(ntext)
            if op == '=':
                assert query[qpos:qpos+n] == reference[rpos:rpos+n]
                qpos += n
                rpos += n
                matches += n
            elif op == 'X':
                assert all(a != b for a, b in zip(query[qpos:qpos+n], reference[rpos:rpos+n]))
                qpos += n
                rpos += n
                errors += n
                score -= 4 * n
            elif op == 'I':
                qpos += n
                errors += n
                score -= 6 + n
            elif op == 'D':
                rpos += n
                errors += n
                score -= 6 + n
        assert (qpos, rpos) == (len(query), len(reference))
        assert [int(fields[9]), int(fields[10])] == [matches, matches + errors]
        assert int(tags['NM']) == errors
        assert int(tags['AS']) == score
        assert abs(float(tags['dv']) - errors / (matches + errors)) < 1e-9
        records.append({'query': fields[0], 'reference': fields[5],
                        'cigar': tags['cg'], 'edit': errors, 'score': score})
    assert len(records) == 6
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', required=True, type=pathlib.Path)
    parser.add_argument('--manifest-output', type=pathlib.Path)
    args = parser.parse_args()
    dest = args.destination.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    archive = dest / PACKAGE_URL.rsplit('/', 1)[1]
    if not archive.exists():
        partial = archive.with_suffix(archive.suffix + '.partial')
        with urllib.request.urlopen(PACKAGE_URL, timeout=120) as response, partial.open('wb') as out:
            for block in iter(lambda: response.read(1024 * 1024), b''):
                out.write(block)
        if partial.stat().st_size != PACKAGE_SIZE or sha256(partial) != PACKAGE_SHA256:
            raise ValueError('Downloaded package checksum/size mismatch')
        partial.replace(archive)
    if archive.stat().st_size != PACKAGE_SIZE or sha256(archive) != PACKAGE_SHA256:
        raise ValueError('Existing package checksum/size mismatch')
    with tarfile.open(archive) as handle:
        handle.extractall(dest / 'prefix', filter='data')
    binary = dest / 'prefix/bin/locityper'
    if sha256(binary) != BINARY_SHA256:
        raise ValueError('Executable checksum mismatch')
    version = subprocess.run([str(binary), 'version'], capture_output=True, text=True, check=True)
    assert version.stdout.splitlines()[0] == 'locityper v0.17.3'
    help_result = subprocess.run([str(binary), 'align', '-h'], capture_output=True, text=True, check=True)
    (dest / 'align_help.txt').write_text(help_result.stdout)

    # A deterministic random template, with three controlled edits. These are
    # synthetic sequences and contain no actual study candidate haplotypes.
    rng = random.Random(173)
    reference = ''.join(rng.choice('ACGT') for _ in range(1200))
    insertion = ''.join(rng.choice('ACGT') for _ in range(60))
    substitution = next(base for base in 'ACGT' if base != reference[500])
    sequences = {
        'synthetic_reference': reference,
        'synthetic_insertion60': reference[:600] + insertion + reference[600:],
        'synthetic_deletion80': reference[:400] + reference[480:],
        'synthetic_substitution1': reference[:500] + substitution + reference[501:],
    }
    fasta, paf = dest / 'synthetic.fa', dest / 'synthetic.paf'
    fasta.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in sequences.items()))
    command = [str(binary), 'align', '-i', str(fasta), '-o', str(paf), '-A'] + ALIGN_FLAGS
    run = subprocess.run(command, capture_output=True, text=True, check=True)
    (dest / 'synthetic.log').write_text(run.stdout + run.stderr)
    records = verify_paf(paf, sequences)
    manifest = {
        'package_url': PACKAGE_URL,
        'package_sha256': PACKAGE_SHA256,
        'package_size_bytes': PACKAGE_SIZE,
        'binary_sha256': BINARY_SHA256,
        'source_tag': 'v0.17.3',
        'source_commit': SOURCE_COMMIT,
        'published_alignment_command': 'locityper align -i loci/$locus/haplotypes.fa.gz -o loci/$locus/haplotypes.paf.gz -A -D 1',
        'published_command_source': {
            'repository': 'https://github.com/tprodanov/locityper-bench',
            'commit': 'ed5c5b4c6f7495287a9fff6b62eb7666409ef2e2',
            'file': 'README.md',
            'file_sha256': '363f616250840f94fe062c7f01b57763c3d20be7869b536616aea540992b1f02',
        },
        'source_archive_sha256_from_recipe': 'd9ecc464a4394a0523afd6ba4c5f7bf80fe44ca15926cca1075ac70cd14b9255',
        'version_output': version.stdout.strip(),
        'alignment_flags': ALIGN_FLAGS,
        'synthetic_input_sha256': sha256(fasta),
        'synthetic_paf_sha256': sha256(paf),
        'synthetic_records': records,
        'limitations': [
            'Release-matched Bioconda build; not the authors original benchmark executable.',
            'The packaged recipe cloned WFA2 without pinning its source commit.',
            'Synthetic smoke validation does not establish identical outputs on study haplotypes.',
            'Default backbone/max-gap alignment is heuristic; accuracy 9 does not certify global optimality.',
        ],
    }
    output = args.manifest_output or dest / 'runtime_verification.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'runtime': str(binary), 'synthetic_alignments_validated': len(records),
                      'manifest': str(output)}, indent=2))


if __name__ == '__main__':
    main()
