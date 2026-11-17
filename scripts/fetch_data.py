#!/usr/bin/env python3
"""Download the frozen benchmark data, verify it and extract regular files."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

REPO = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--archive', type=Path, help='Use an already downloaded release archive')
parser.add_argument('--data', type=Path, default=REPO/'data')
args = parser.parse_args()
asset = json.loads((REPO/'metadata/data_asset.json').read_text())
archive = args.archive or REPO/'.cache'/asset['asset']
if not archive.exists():
    archive.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--output', str(archive), asset['url']], check=True)
assert hashlib.sha256(archive.read_bytes()).hexdigest() == asset['sha256'], 'Archive SHA256 mismatch'
args.data.mkdir(parents=True, exist_ok=True)
with tarfile.open(archive, 'r:gz') as handle:
    for member in handle.getmembers():
        assert member.isfile() and not Path(member.name).is_absolute() and '..' not in Path(member.name).parts
    handle.extractall(args.data, filter='data')
for line in (args.data/'DATA_SHA256SUMS').read_text().splitlines():
    digest, name = line.split('  ', 1)
    assert hashlib.sha256((args.data/name).read_bytes()).hexdigest() == digest, name
print(f"Verified {asset['files']} source files in {args.data}")
