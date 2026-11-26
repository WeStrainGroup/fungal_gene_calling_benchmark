#!/usr/bin/env python3
"""CTC: Combined Tiberius and Carbon-A, routed by FASTA record length."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shlex
import subprocess
import time


def fasta(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        header, parts = None, []
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header, parts = line[1:].strip(), []
                if not header:
                    raise ValueError('Empty FASTA identifier')
            elif line.strip():
                if header is None:
                    raise ValueError('Sequence before first FASTA header')
                parts.append(''.join(line.split()))
        if header is not None:
            yield header, ''.join(parts)


def split_fasta(source, output, minimum, cutoff):
    counts = {t:dict(records=0, bp=0) for t in ['tiberius', 'carbon', 'skipped']}
    seen = set()
    with (output/'tiberius/input.fna').open('w') as short, \
         (output/'carbon/input.fna').open('w') as long, \
         (output/'routing.tsv').open('w') as report:
        writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['sequence', 'length_bp', 'tool'])
        for header, seq in fasta(source):
            name = header.split()[0]
            if name in seen:
                raise ValueError(f'Duplicate FASTA identifier: {name}')
            seen.add(name)
            tool = 'skipped' if len(seq) < minimum else ('tiberius' if len(seq) <= cutoff else 'carbon')
            counts[tool]['records'] += 1
            counts[tool]['bp'] += len(seq)
            writer.writerow([name, len(seq), tool])
            if tool != 'skipped':
                (short if tool == 'tiberius' else long).write(f'>{header}\n{seq}\n')
    if not seen:
        raise ValueError('Empty input FASTA')
    return counts


def merge_gffs(sources, output):
    """Namespace model IDs; leave coordinates, phase and other attributes intact."""
    cds = 0
    temporary = output.with_suffix('.tmp')
    with temporary.open('w') as target:
        target.write('##gff-version 3\n')
        for tool, path in sources:
            with path.open() as source:
                for line in source:
                    if line.startswith('##FASTA'):
                        break
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) != 9:
                        raise ValueError(f'Malformed GFF row in {path}')
                    attributes = []
                    for attribute in fields[8].split(';'):
                        key, sep, value = attribute.partition('=')
                        if key in {'ID', 'Parent', 'Derives_from'}:
                            value = ','.join(tool+'_'+x for x in value.split(','))
                        attributes.append(key+sep+value)
                    fields[8] = ';'.join(attributes)
                    target.write('\t'.join(fields)+'\n')
                    cds += fields[2] == 'CDS'
    temporary.replace(output)
    return cds


def extra_args(text, reserved):
    args = shlex.split(text)
    if any(x.split('=')[0] in reserved for x in args):
        raise ValueError('Input/output/model options are managed by CTC; use its named options')
    return args


def run(a):
    started = time.monotonic()
    output = a.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for tool in ['tiberius', 'carbon']:
        (output/tool).mkdir()
    report = dict(strategy='CTC', input=str(a.input.resolve()), min_length=a.min_length,
                  cutoff=a.cutoff, status='running', components=[])
    try:
        split_start = time.monotonic()
        report['routing'] = split_fasta(a.input, output, a.min_length, a.cutoff)
        report['split_seconds'] = time.monotonic()-split_start
        sources = []
        for tool in ['tiberius', 'carbon']:
            if not report['routing'][tool]['records']:
                continue
            dest = output/tool
            env = os.environ.copy()
            python = (a.tiberius_python if tool == 'tiberius' else a.carbon_python).absolute()
            env['PATH'] = str(python.parent)+os.pathsep+env.get('PATH', '')
            env['LD_LIBRARY_PATH'] = str(python.parent.parent/'lib')+os.pathsep+env.get('LD_LIBRARY_PATH', '')
            if tool == 'tiberius':
                command = [str(python), '-m', 'tiberius.main', '--genome', str(dest/'input.fna'),
                    '--model', str(a.tiberius_model.resolve()), '--out', str(dest/'prediction.gff3'),
                    '--codingseq', str(dest/'native.cds.fa'), '--protseq', str(dest/'native.proteins.fa'),
                    '--min_genome_seqlen', str(a.min_length)] + extra_args(a.tiberius_args,
                    {'--genome','--model','--out','--codingseq','--protseq','--min_genome_seqlen'})
                gff = dest/'prediction.gff3'
            else:
                command = ['bash', str(a.carbon_dir.resolve()/'annotation_pipeline/run_annotation.sh'),
                    '--input', str(dest/'input.fna'), '--model', str((a.carbon_model or a.carbon_dir).resolve()),
                    '--output-dir', str(dest/'work')] + extra_args(a.carbon_args,
                    {'--input','--model','--output-dir','--codon-table'})
                if not any(x.split('=')[0] == '--min-record-length' for x in command):
                    command += ['--min-record-length', '0']
                if not any(x.split('=')[0] == '--gpu-count' for x in command):
                    command += ['--gpu-count', '1']
                if a.codon_table:
                    command += ['--codon-table', str(a.codon_table.resolve())]
                gff = dest/'work/prediction.gff'
            part = dict(tool=tool, command=command)
            report['components'].append(part)
            part_start = time.monotonic()
            with (dest/'run.log').open('w') as log:
                completed = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            part.update(elapsed_seconds=time.monotonic()-part_start, exit_code=completed.returncode)
            if completed.returncode:
                raise RuntimeError(f'{tool} failed; see {dest}/run.log')
            if tool == 'carbon' and not gff.is_file():
                gff = dest/'work_WARNING/prediction.gff'
            if not gff.is_file() or not gff.stat().st_size:
                raise RuntimeError(f'Missing native GFF: {gff}')
            sources.append((tool, gff))
        merge_start = time.monotonic()
        report['cds_features'] = merge_gffs(sources, output/'prediction.gff3')
        report['merge_seconds'] = time.monotonic()-merge_start
        report['status'] = 'complete'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['elapsed_seconds'] = time.monotonic()-started
        (output/'run.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version='CTC 0.1.0')
    p.add_argument('--input', required=True, type=Path, help='FASTA or FASTA.gz')
    p.add_argument('--output', required=True, type=Path, help='New output directory')
    p.add_argument('--min-length', type=int, default=1000)
    p.add_argument('--cutoff', type=int, default=10000, help='Inclusive upper length for Tiberius; longer records go to Carbon-A')
    p.add_argument('--tiberius-python', required=True, type=Path)
    p.add_argument('--tiberius-model', required=True, type=Path)
    p.add_argument('--tiberius-args', default='', help='Quoted additional native options; use --tiberius-args="..."')
    p.add_argument('--carbon-python', required=True, type=Path, help='Python in the Carbon-A environment; its bin directory is placed first in PATH')
    p.add_argument('--carbon-dir', required=True, type=Path)
    p.add_argument('--carbon-model', type=Path, help='Defaults to --carbon-dir')
    p.add_argument('--carbon-args', default='', help='Quoted additional native options; use --carbon-args="..."')
    p.add_argument('--codon-table', type=Path, help='Native Carbon-A codon-table JSON; not passed to Tiberius')
    a = p.parse_args()
    if not 1 <= a.min_length <= a.cutoff:
        p.error('Require 1 <= min-length <= cutoff')
    try:
        result = run(a)
    except (ValueError, OSError, RuntimeError) as error:
        p.exit(1, f'CTC: {error}\n')
    print(f"CTC complete: {a.output}/prediction.gff3 ({result['elapsed_seconds']:.1f} s)")


if __name__ == '__main__':
    main()
