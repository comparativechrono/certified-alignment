#!/usr/bin/env python3
"""Download/verify the pinned official runtime and demo; never run genotyping.

Requires Python >=3.12 (tar extraction filtering). The default downloads about
145 MB compressed. Files remain outside the repository at an explicit destination.
The native runtime requires Linux x86_64. No Docker daemon is used.
"""

import argparse
import hashlib
import json
import pathlib
import shlex
import subprocess
import tarfile
import urllib.request

MANIFEST = 'sha256:a9575354804c47be19b4cb5da74f551d3f018ebaaef86da501808ae20d4abbae'
REGISTRY = 'https://registry-1.docker.io/v2/tprodanov/locityper/'
DATA_URL = 'https://zenodo.org/records/21903460/files/data.tar.xz?download=1'
DATA_SHA = 'f5c199e31723b361a5200e409a35b65e0df98bc7cc3e0fb84322ed48f72fe495'
DATA_SIZE = 107320572


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(url, path, size, digest, headers=None):
    if path.exists():
        if path.stat().st_size == size and sha256(path) == digest:
            return
        raise ValueError(f'Existing file has wrong checksum: {path}')
    request = urllib.request.Request(url, headers=headers or {})
    temp = path.with_suffix(path.suffix + '.partial')
    with urllib.request.urlopen(request, timeout=120) as response, temp.open('wb') as out:
        for block in iter(lambda: response.read(1024 * 1024), b''):
            out.write(block)
    if temp.stat().st_size != size or sha256(temp) != digest:
        raise ValueError(f'Download has wrong size or checksum: {url}')
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=pathlib.Path, required=True)
    parser.add_argument('--skip-data', action='store_true')
    args = parser.parse_args()
    dest = args.destination.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    image = dest / 'locityper_image'
    image.mkdir(exist_ok=True)
    token_url = ('https://auth.docker.io/token?service=registry.docker.io'
                 '&scope=repository:tprodanov/locityper:pull')
    with urllib.request.urlopen(token_url, timeout=60) as response:
        token = json.load(response)['token']
    headers = {'Authorization': 'Bearer ' + token,
               'Accept': 'application/vnd.oci.image.manifest.v1+json'}
    request = urllib.request.Request(REGISTRY + 'manifests/' + MANIFEST, headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        raw = response.read()
    if hashlib.sha256(raw).hexdigest() != MANIFEST.split(':', 1)[1]:
        raise ValueError('Manifest checksum mismatch')
    (image / 'manifest_pinned.json').write_bytes(raw)
    manifest = json.loads(raw)
    rootfs = image / 'rootfs'
    for layer in manifest['layers']:
        digest = layer['digest'].split(':', 1)[1]
        archive = image / (digest + '.tar.gz')
        download(REGISTRY + 'blobs/' + layer['digest'], archive,
                 layer['size'], digest, headers)
        with tarfile.open(archive, 'r:gz') as handle:
            members = [member for member in handle.getmembers()
                       if member.name.startswith(('lib/', 'usr/lib/', 'usr/local/'))
                       and (member.isfile() or member.issym())]
            handle.extractall(rootfs, members=members, filter='data')
        print('Verified image layer', digest, flush=True)

    bindir = image / 'bin'
    bindir.mkdir(exist_ok=True)
    loader = rootfs / 'lib/ld-musl-x86_64.so.1'
    libpath = ':'.join(str(rootfs / part) for part in ('lib', 'usr/lib', 'usr/local/lib'))
    versions = []
    for name in ('locityper', 'samtools', 'strobealign', 'jellyfish', 'minimap2'):
        binary = rootfs / 'usr/local/bin' / name
        command = [str(binary)] if name == 'strobealign' else [
            str(loader), '--library-path', libpath, str(binary)]
        wrapper = bindir / name
        wrapper.write_text('#!/bin/sh\nexec ' + shlex.join(command) + ' "$@"\n')
        wrapper.chmod(0o755)
        argv = [str(wrapper), 'version' if name == 'locityper' else '--version']
        result = subprocess.run(argv, text=True, capture_output=True, check=True)
        versions.append({'executable': name, 'sha256': sha256(binary),
                         'returncode': result.returncode,
                         'output': (result.stdout + result.stderr).strip()})
        print(versions[-1]['output'].splitlines()[0], flush=True)
    (image / 'runtime_verification.json').write_text(json.dumps(versions, indent=2) + '\n')

    if not args.skip_data:
        data = dest / 'locityper_test_v2'
        data.mkdir(exist_ok=True)
        archive = data / 'data.tar.xz'
        download(DATA_URL, archive, DATA_SIZE, DATA_SHA)
        with tarfile.open(archive, 'r:xz') as handle:
            handle.extractall(data, filter='data')
        print('Verified official HG00733 demo archive', flush=True)
    print('Preparation complete. No genotypes were called.', flush=True)


if __name__ == '__main__':
    main()
