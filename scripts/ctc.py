#!/usr/bin/env python3
"""CTC: Combined Tiberius and Carbon-A, routed by FASTA record length."""
import argparse
import csv
import gzip
import hashlib
from itertools import product
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


def codon_table(path=None):
    codons = [''.join(c) for c in product('TCAG', repeat=3)]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    table = dict(zip(codons, amino_acids)) if path is None else json.loads(path.read_text())
    if (not isinstance(table, dict) or set(table) != set(codons) or
        any(not isinstance(aa, str) or len(aa) != 1 or aa not in 'ACDEFGHIKLMNPQRSTVWYOU*'
            for aa in table.values())):
        raise ValueError('Codon table must map all 64 uppercase DNA codons to one amino acid or *')
    return table


def translate_cds(source, output, diagnostics, table, transcript_ids=None):
    """Translate frame 1, retain defects, and remove only one terminal stop."""
    iupac = dict(A='A', C='C', G='G', T='T', R='AG', Y='CT', S='GC', W='AT',
                 K='GT', M='AC', B='CGT', D='AGT', H='ACT', V='ACG', N='ACGT')
    seen = set()
    summary = dict(records=0, internal_stops=0, unknown_residues=0, partial_codons=0)
    with output.open('w') as proteins, diagnostics.open('w') as report:
        writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['sequence', 'native_header', 'cds_bp', 'protein_aa', 'ambiguous_codons',
                         'trailing_bases', 'internal_stops', 'unknown_residues'])
        for header, seq in fasta(source):
            name = header.split()[0]
            protein_header = header
            if transcript_ids is not None:
                parts = name.split('|')
                transcript = name if name in transcript_ids else (parts[1] if len(parts) > 1 else name)
                if transcript not in transcript_ids or transcript in seen:
                    raise ValueError(f'Unmatched or duplicate Tiberius CDS identifier: {name}')
                seen.add(transcript)
                name = 'tiberius_'+transcript
                protein_header = name+' '+header
            seq = seq.upper().replace('U', 'T')
            residues, ambiguous = [], 0
            for i in range(0, len(seq)-2, 3):
                codon = seq[i:i+3]
                if codon in table:
                    residues.append(table[codon])
                else:
                    ambiguous += 1
                    possibilities = {table[''.join(c)] for c in product(*(iupac.get(b, '') for b in codon))}
                    residues.append(next(iter(possibilities)) if len(possibilities) == 1 else 'X')
            trailing = len(seq) % 3
            if trailing:
                residues.append('X')
            protein = ''.join(residues)
            if protein.endswith('*'):
                protein = protein[:-1]
            stops, unknown = protein.count('*'), protein.count('X')
            proteins.write(f'>{protein_header}\n{protein}\n')
            writer.writerow([name, header, len(seq), len(protein), ambiguous, trailing, stops, unknown])
            summary['records'] += 1
            summary['internal_stops'] += stops
            summary['unknown_residues'] += unknown
            summary['partial_codons'] += bool(trailing)
    if transcript_ids is not None and seen != transcript_ids:
        raise ValueError('Some native Tiberius CDS models have no matching FASTA record')
    return summary


def cds_parents(gff):
    ids = set()
    with gff.open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if len(fields) == 9 and fields[2] == 'CDS':
                attrs = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
                ids.update(attrs['Parent'].split(','))
    return ids


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
    report = dict(strategy='CTC', version='0.2.0', input=str(a.input.resolve()), min_length=a.min_length,
                  cutoff_length=a.cutoff_length, status='running', components=[])
    try:
        table = codon_table(a.codon_table)
        effective_table = output/'codon_table.json'
        effective_table.write_text(json.dumps(table, indent=2, sort_keys=True)+'\n')
        report['codon_table'] = dict(source=str(a.codon_table.resolve()) if a.codon_table else 'NCBI 1',
            sha256=hashlib.sha256(effective_table.read_bytes()).hexdigest())
        split_start = time.monotonic()
        report['routing'] = split_fasta(a.input, output, a.min_length, a.cutoff_length)
        report['split_seconds'] = time.monotonic()-split_start
        sources = []
        for tool in ['tiberius', 'carbon']:
            if not report['routing'][tool]['records']:
                continue
            dest = output/tool
            env = os.environ.copy()
            python = (a.tiberius_path if tool == 'tiberius' else a.carbon_a_path).absolute()
            env['PATH'] = str(python.parent)+os.pathsep+env.get('PATH', '')
            env['LD_LIBRARY_PATH'] = str(python.parent.parent/'lib')+os.pathsep+env.get('LD_LIBRARY_PATH', '')
            if tool == 'tiberius':
                command = [str(python), '-m', 'tiberius.main', '--genome', str(dest/'input.fna'),
                    '--model', str(a.tiberius_model.resolve()), '--out', str(dest/'prediction.gff3'),
                    '--codingseq', str(dest/'native.cds.fa'),
                    '--min_genome_seqlen', str(a.min_length)] + extra_args(a.tiberius_args,
                    {'--genome','--model','--out','--codingseq','--protseq','--min_genome_seqlen'})
                gff = dest/'prediction.gff3'
            else:
                command = ['bash', str(a.carbon_a_model.resolve()/'annotation_pipeline/run_annotation.sh'),
                    '--input', str(dest/'input.fna'), '--model', str(a.carbon_a_model.resolve()),
                    '--output-dir', str(dest/'work')] + extra_args(a.carbon_a_args,
                    {'--input','--model','--output-dir','--codon-table'})
                if not any(x.split('=')[0] == '--min-record-length' for x in command):
                    command += ['--min-record-length', '0']
                if not any(x.split('=')[0] == '--gpu-count' for x in command):
                    command += ['--gpu-count', '1']
                command += ['--codon-table', str(effective_table)]
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
            if tool == 'tiberius':
                translation_start = time.monotonic()
                part['translation'] = translate_cds(dest/'native.cds.fa', dest/'proteins.fa',
                    dest/'translation.tsv', table, cds_parents(gff))
                part['translation_seconds'] = time.monotonic()-translation_start
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
    p.add_argument('--version', action='version', version='CTC 0.2.0')
    p.add_argument('--input', required=True, type=Path, help='FASTA or FASTA.gz')
    p.add_argument('--output', required=True, type=Path, help='New output directory')
    p.add_argument('--min-length', type=int, default=1000)
    p.add_argument('--cutoff-length', type=int, default=10000, help='Inclusive upper length for Tiberius; longer records go to Carbon-A')
    p.add_argument('--tiberius-path', required=True, type=Path, help='Python executable in the Tiberius environment')
    p.add_argument('--tiberius-model', required=True, type=Path)
    p.add_argument('--tiberius-args', default='', help='Quoted additional native options; use --tiberius-args="..."')
    p.add_argument('--carbon-a-path', required=True, type=Path, help='Python in the Carbon-A environment; its bin directory is placed first in PATH')
    p.add_argument('--carbon-a-model', required=True, type=Path, help='Carbon-A directory containing model weights and annotation_pipeline/')
    p.add_argument('--carbon-a-args', default='', help='Quoted additional native options; use --carbon-a-args="..."')
    p.add_argument('--codon-table', type=Path, help='64-codon JSON shared by Carbon-A and Tiberius CDS translation; default: standard code (NCBI 1)')
    a = p.parse_args()
    if not 1 <= a.min_length <= a.cutoff_length:
        p.error('Require 1 <= min-length <= cutoff-length')
    try:
        result = run(a)
    except (ValueError, OSError, RuntimeError) as error:
        p.exit(1, f'CTC: {error}\n')
    print(f"CTC complete: {a.output}/prediction.gff3 ({result['elapsed_seconds']:.1f} s)")


if __name__ == '__main__':
    main()
