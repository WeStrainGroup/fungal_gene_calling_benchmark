#!/usr/bin/env python3
"""CTC: Combined Tiberius and Carbon-A, routed by FASTA record length."""
import argparse
import csv
import gzip
import hashlib
from itertools import product
import json
import math
import os
from pathlib import Path
import shlex
import shutil
from urllib.parse import unquote
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
    table = dict(zip(codons, amino_acids))
    if path is not None:
        if str(path) in {'1','12','26'}:
            table['CTG'] = {'1':'L','12':'S','26':'A'}[str(path)]
        else:
            table = json.loads(Path(path).read_text())
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
    counts = {t:dict(records=0, bp=0) for t in ['tiberius', 'carbon_a', 'skipped']}
    seen = set()
    with (output/'tiberius/input.fna').open('w') as short, \
         (output/'carbon_a/input.fna').open('w') as long, \
         (output/'combined/routing.tsv').open('w') as report:
        writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['sequence', 'length_bp', 'tool'])
        for header, seq in fasta(source):
            name = header.split()[0]
            if name in seen:
                raise ValueError(f'Duplicate FASTA identifier: {name}')
            seen.add(name)
            tool = 'skipped' if len(seq) < minimum else ('tiberius' if len(seq) <= cutoff else 'carbon_a')
            counts[tool]['records'] += 1
            counts[tool]['bp'] += len(seq)
            writer.writerow([name, len(seq), tool])
            if tool != 'skipped':
                (short if tool == 'tiberius' else long).write(f'>{header}\n{seq}\n')
    if not seen:
        raise ValueError('Empty input FASTA')
    return counts


def read_gff(path):
    rows, features = [], {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Malformed GFF3 row in {path}')
            attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
            row = dict(fields=fields, attrs=attrs, index=len(rows))
            rows.append(row)
            if 'ID' in attrs:
                features.setdefault(unquote(attrs['ID']), []).append(row)
    for key, parts in features.items():
        signatures = {(r['fields'][0], r['fields'][2], r['fields'][6], r['attrs'].get('Parent')) for r in parts}
        coordinates = {(r['fields'][3], r['fields'][4]) for r in parts}
        if len(signatures) != 1 or len(coordinates) != len(parts):
            raise ValueError(f'Duplicate independent feature ID: {key}')
        if len(parts) > 1 and parts[0]['fields'][2] in {'gene', 'mRNA', 'transcript'}:
            raise ValueError(f'Duplicate gene/transcript ID: {key}')
    for row in rows:
        for attribute in ['Parent', 'Derives_from']:
            for ref in references(row, attribute):
                if ref not in features:
                    raise ValueError(f'Dangling {attribute}: {ref}')
    return rows, features


def references(row, attribute='Parent'):
    return [unquote(x) for x in row['attrs'].get(attribute, '').split(',') if x]


def fasta_index(path):
    records = {}
    for header, seq in fasta(path):
        key = header.split()[0]
        if key in records:
            raise ValueError(f'Duplicate FASTA ID: {key}')
        records[key] = (header, seq)
    return records


def native_models(tool, directory, include_incomplete=False, confidence_threshold=0.1):
    rows, features = read_gff(directory/'prediction.gff3')
    genes = {key:parts[0] for key,parts in features.items() if parts[0]['fields'][2] == 'gene'}
    transcripts = {key:parts[0] for key,parts in features.items() if parts[0]['fields'][2] in {'mRNA','transcript'}}
    cds_ids = {p for row in rows if row['fields'][2] == 'CDS' for p in references(row)}
    if cds_ids != set(transcripts):
        raise ValueError(f'{tool}: CDS and transcript model sets differ')
    owner = {}
    def owners(key, visiting=()):
        if key in visiting:
            raise ValueError('Cyclic GFF3 Parent references')
        if key not in owner:
            owner[key] = {key} if key in genes else set().union(
                *(owners(p, visiting+(key,)) for p in references(features[key][0])))
        return owner[key]
    for key in features:
        owners(key)
    gene_transcripts = {}
    for transcript in transcripts:
        for gene in owners(transcript):
            gene_transcripts.setdefault(gene, []).append(transcript)
    cds, proteins = fasta_index(directory/'cds.fa'), fasta_index(directory/'proteins.fa')
    if set(cds) != set(proteins):
        raise ValueError(f'{tool}: CDS/protein FASTA IDs differ')
    model_records = {}
    for original, record in cds.items():
        if tool == 'tiberius':
            parts = original.split('|')
            candidates = [original] if original in transcripts else parts[1:2]
        else:
            candidates = gene_transcripts.get(original, [])
            if original in transcripts:
                candidates = [original]
        if len(candidates) != 1 or candidates[0] not in transcripts or candidates[0] in model_records:
            raise ValueError(f'{tool}: ambiguous or unmatched native FASTA ID: {original}')
        model_records[candidates[0]] = (original, record, proteins[original])
    if set(model_records) != cds_ids:
        raise ValueError(f'{tool}: missing FASTA records for native CDS models')
    allowed = set(genes)
    incomplete_excluded = confidence_excluded = 0
    if tool == 'carbon_a':
        for key, gene in genes.items():
            status = gene['attrs'].get('status')
            if status not in {'complete', 'incomplete'}:
                raise ValueError('Carbon-A gene lacks an explicit complete/incomplete status')
            raw = gene['attrs'].get('confidence')
            if raw is None:
                raise ValueError(f'Carbon-A gene {key} lacks confidence (use explicit NA for unscored models)')
            score = None if raw == 'NA' else float(raw)
            if score is not None and (not math.isfinite(score) or not 0 <= score <= 1):
                raise ValueError(f'Carbon-A gene {key} has invalid confidence: {raw}')
            if status == 'incomplete' and not include_incomplete:
                allowed.remove(key)
                incomplete_excluded += 1
            elif score is not None and score < confidence_threshold:
                allowed.remove(key)
                confidence_excluded += 1
    selected = []
    for row in rows:
        roots = owners(unquote(row['attrs']['ID'])) if 'ID' in row['attrs'] else set().union(
            *(owners(p) for p in references(row)))
        if not roots:
            raise ValueError(f'{tool}: feature has no gene ancestor')
        if roots & allowed and not roots <= allowed:
            raise ValueError('A shared feature spans retained and excluded genes')
        if roots <= allowed:
            selected.append(row)
    return dict(tool=tool, rows=selected, genes={g:genes[g] for g in allowed},
                transcripts={t:r for t,r in transcripts.items() if owners(t) <= allowed},
                records=model_records, incomplete_excluded=incomplete_excluded,
                confidence_excluded=confidence_excluded)


def combine_models(sources, output, include_incomplete=False, confidence_threshold=0.1):
    """Merge a single selected model set, preserving native IDs, headers and sequences."""
    if not math.isfinite(confidence_threshold) or not 0 <= confidence_threshold <= 1:
        raise ValueError('confidence-threshold must be between 0 and 1')
    datasets = [native_models(tool, path, include_incomplete, confidence_threshold) for tool,path in sources]
    feature_ids, sequence_ids, merged_rows, records = set(), set(), [], []
    for data in datasets:
        current = {unquote(r['attrs']['ID']) for r in data['rows'] if 'ID' in r['attrs']}
        if current & feature_ids:
            raise ValueError('Native GFF IDs collide between tools; CTC does not rename them')
        feature_ids.update(current)
        for row in data['rows']:
            for attribute in ['Parent','Derives_from']:
                if any(ref not in current for ref in references(row,attribute)):
                    raise ValueError(f'{attribute} points outside the selected model set')
        merged_rows.extend(data['rows'])
        for transcript in data['transcripts']:
            original, cds, protein = data['records'][transcript]
            if original in sequence_ids:
                raise ValueError('Native FASTA IDs collide between tools; CTC does not rename them')
            sequence_ids.add(original)
            records.append((cds, protein))
    with (output/'prediction.gff3').open('w') as handle:
        handle.write('##gff-version 3\n')
        for row in merged_rows:
            handle.write('\t'.join(row['fields'])+'\n')
    read_gff(output/'prediction.gff3')
    with (output/'cds.fa').open('w') as cds_file, (output/'proteins.fa').open('w') as proteins_file:
        for (cds_header, dna), (protein_header, aa) in records:
            cds_file.write(f'>{cds_header}\n{dna}\n')
            proteins_file.write(f'>{protein_header}\n{aa}\n')
    return dict(genes=sum(len(d['genes']) for d in datasets),models=len(records),
                cds_features=sum(r['fields'][2]=='CDS' for r in merged_rows),
                carbon_a_incomplete_excluded=sum(d['incomplete_excluded'] for d in datasets),
                carbon_a_confidence_excluded=sum(d['confidence_excluded'] for d in datasets))


def extra_args(text, reserved):
    args = shlex.split(text)
    if any(x.split('=')[0] in reserved for x in args):
        raise ValueError('Input/output/model options are managed by CTC; use its named options')
    return args


def run(a):
    started = time.monotonic()
    output = a.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for tool in ['tiberius', 'carbon_a']:
        (output/tool).mkdir()
    (output/'combined').mkdir()
    report = dict(strategy='CTC', version='0.4.0', input=str(a.input.resolve()), min_length=a.min_length,
                  cutoff_length=a.cutoff_length, include_incomplete=a.include_incomplete,
                  confidence_threshold=a.confidence_threshold, status='running', components=[])
    try:
        table = codon_table(a.codon_table)
        effective_table = output/'combined/codon_table.json'
        effective_table.write_text(json.dumps(table, indent=2, sort_keys=True)+'\n')
        report['codon_table'] = dict(source=str(a.codon_table) if a.codon_table else 'NCBI 1',
            sha256=hashlib.sha256(effective_table.read_bytes()).hexdigest())
        split_start = time.monotonic()
        report['routing'] = split_fasta(a.input, output, a.min_length, a.cutoff_length)
        report['split_seconds'] = time.monotonic()-split_start
        sources = []
        for tool in ['tiberius', 'carbon_a']:
            if not report['routing'][tool]['records']:
                (output/tool/'prediction.gff3').write_text('##gff-version 3\n')
                (output/tool/'cds.fa').touch()
                (output/tool/'proteins.fa').touch()
                continue
            dest = output/tool
            env = os.environ.copy()
            python = (a.tiberius_path if tool == 'tiberius' else a.carbon_a_path).absolute()
            env['PATH'] = str(python.parent)+os.pathsep+env.get('PATH', '')
            env['LD_LIBRARY_PATH'] = str(python.parent.parent/'lib')+os.pathsep+env.get('LD_LIBRARY_PATH', '')
            if tool == 'tiberius':
                command = [str(python), '-m', 'tiberius.main', '--genome', str(dest/'input.fna'),
                    '--model', str(a.tiberius_model.resolve()), '--out', str(dest/'prediction.gff3'),
                    '--codingseq', str(dest/'cds.fa'),
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
            if tool == 'carbon_a' and not gff.is_file():
                gff = dest/'work_WARNING/prediction.gff'
            if not gff.is_file() or not gff.stat().st_size:
                raise RuntimeError(f'Missing native GFF: {gff}')
            if tool == 'tiberius':
                translation_start = time.monotonic()
                part['translation'] = translate_cds(dest/'cds.fa', dest/'proteins.fa',
                    dest/'translation.tsv', table, cds_parents(gff))
                part['translation_seconds'] = time.monotonic()-translation_start
            else:
                for source_name, target in [('prediction.gff','prediction.gff3'),('prediction.fna','cds.fa'),('prediction.faa','proteins.fa')]:
                    shutil.move(str(gff.parent/source_name), str(dest/target))
            sources.append((tool, dest))
        merge_start = time.monotonic()
        report['combined'] = combine_models(sources, output/'combined', a.include_incomplete, a.confidence_threshold)
        report['cds_features'] = report['combined']['cds_features']
        report['merge_seconds'] = time.monotonic()-merge_start
        report['status'] = 'complete'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['elapsed_seconds'] = time.monotonic()-started
        (output/'combined/run.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--include-incomplete', action='store_true', help='Also merge Carbon-A status=incomplete models; native tool outputs always retain them')
    p.add_argument('--confidence-threshold', type=float, default=0.1,
                   help='Minimum Carbon-A confidence for combined outputs (inclusive; default: 0.1); NA is retained, 0 disables score filtering')
    p.add_argument('--version', action='version', version='CTC 0.4.0')
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
    p.add_argument('--codon-table', help='NCBI code 1/12/26 or 64-codon JSON shared by Carbon-A and Tiberius CDS translation; default: standard code (NCBI 1)')
    a = p.parse_args()
    if not 1 <= a.min_length <= a.cutoff_length:
        p.error('Require 1 <= min-length <= cutoff-length')
    if not math.isfinite(a.confidence_threshold) or not 0 <= a.confidence_threshold <= 1:
        p.error('confidence-threshold must be between 0 and 1')
    try:
        result = run(a)
    except (ValueError, OSError, RuntimeError) as error:
        p.exit(1, f'CTC: {error}\n')
    print(f"CTC complete: {a.output}/combined/prediction.gff3 ({result['elapsed_seconds']:.1f} s)")


if __name__ == '__main__':
    main()
