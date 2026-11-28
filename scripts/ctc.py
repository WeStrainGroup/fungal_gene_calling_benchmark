#!/usr/bin/env python3
"""CTC: Combined Tiberius and Carbon-A, routed by FASTA record length."""
import argparse
import csv
import gzip
import hashlib
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


def codon_table(code=1):
    """Return an NCBI code that the native Carbon-A pipeline can represent."""
    from Bio.Data import CodonTable
    code = int(code)
    if code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f'Unknown or discontinued NCBI genetic code: {code}')
    native = CodonTable.unambiguous_dna_by_id[code]
    if set(native.stop_codons) & set(native.forward_table):
        raise ValueError(f'NCBI code {code} has context-dependent stop codons; Carbon-A requires a fixed codon map')
    unsupported = sorted(set(native.stop_codons) - {'TAA','TAG','TGA'})
    if unsupported:
        raise ValueError(f'NCBI code {code} uses stop codons unsupported by Carbon-A: {", ".join(unsupported)}')
    return {**native.forward_table, **dict.fromkeys(native.stop_codons, '*')}


def translate_cds(source, output, diagnostics, code=1, phases=None):
    """Translate joined CDS once, applying only its initial GFF phase."""
    from Bio.Seq import Seq
    seen = set()
    summary = dict(records=0, internal_stops=0, unknown_residues=0, partial_codons=0)
    with output.open('w') as proteins, diagnostics.open('w') as report:
        writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['sequence', 'native_header', 'cds_bp', 'initial_phase', 'protein_aa',
                         'ambiguous_codons', 'trailing_bases', 'internal_stops', 'unknown_residues'])
        for header, seq in fasta(source):
            name = header.split()[0]
            if name in seen or (phases is not None and name not in phases):
                raise ValueError(f'Duplicate or unmatched CDS identifier: {name}')
            seen.add(name)
            phase = 0 if phases is None else phases[name]
            if phase not in {0, 1, 2} or phase > len(seq):
                raise ValueError(f'Invalid initial CDS phase for {name}: {phase}')
            coding = seq.upper().replace('U', 'T')[phase:]
            trailing = len(coding) % 3
            coding = coding[:len(coding)-trailing]
            protein = str(Seq(coding).translate(table=int(code), cds=False, to_stop=False))
            if protein.endswith('*'):
                protein = protein[:-1]
            ambiguous = sum(any(b not in 'ACGT' for b in coding[i:i+3]) for i in range(0,len(coding),3))
            stops, unknown = protein.count('*'), protein.count('X')
            proteins.write(f'>{header}\n{protein}\n')
            writer.writerow([name, header, len(seq), phase, len(protein), ambiguous, trailing, stops, unknown])
            summary['records'] += 1
            summary['internal_stops'] += stops
            summary['unknown_residues'] += unknown
            summary['partial_codons'] += bool(trailing)
    if phases is not None and seen != set(phases):
        raise ValueError('Some selected models have no matching CDS FASTA record')
    return summary


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


def native_models(tool, directory, confidence_threshold=0.1):
    rows, features = read_gff(directory/'prediction.gff3')
    genes = {key:parts[0] for key,parts in features.items() if parts[0]['fields'][2] == 'gene'}
    transcripts = {key:parts[0] for key,parts in features.items() if parts[0]['fields'][2] in {'mRNA','transcript'}}
    cds_rows = {}
    for row in rows:
        if row['fields'][2] == 'CDS':
            for parent in references(row):
                cds_rows.setdefault(parent, []).append(row)
    cds_ids = set(cds_rows)
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
        transcript = candidates[0]
        strand = transcripts[transcript]['fields'][6]
        if strand not in {'+', '-'}:
            raise ValueError(f'{tool}: invalid CDS strand for {transcript}')
        first = min(cds_rows[transcript], key=lambda r:int(r['fields'][3])) if strand == '+' else max(
            cds_rows[transcript], key=lambda r:int(r['fields'][4]))
        if first['fields'][7] not in {'0','1','2'}:
            raise ValueError(f'{tool}: missing or invalid initial phase for {transcript}')
        model_records[transcript] = (original, record, int(first['fields'][7]))
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
            if status == 'complete' and score is None:
                raise ValueError(f'Carbon-A complete gene {key} lacks numeric confidence')
            if status == 'incomplete':
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
    core = [r for r in selected if r['fields'][2] in {'gene','mRNA','transcript','CDS'}]
    return dict(tool=tool, rows=core, auxiliary_rows_omitted=len(selected)-len(core), genes={g:genes[g] for g in allowed},
                transcripts={t:r for t,r in transcripts.items() if owners(t) <= allowed},
                records=model_records, owners=owner, incomplete_excluded=incomplete_excluded,
                confidence_excluded=confidence_excluded)


def rename_models(datasets, record_order):
    """Assign gene/transcript IDs, then remap every GFF3 reference consistently."""
    genes = [(i,key,row) for i,d in enumerate(datasets) for key,row in d['genes'].items()]
    def location(row):
        f = row['fields']
        if record_order is not None:
            if f[0] not in record_order:
                raise ValueError(f'Predicted sequence absent from input FASTA: {f[0]}')
            contig = record_order[f[0]]
        else:
            contig = f[0]
        return contig, int(f[3]), int(f[4]), f[6]
    genes.sort(key=lambda item:(location(item[2]),datasets[item[0]]['tool'],item[1]))
    gene_rank = {(i,key):rank for rank,(i,key,_) in enumerate(genes,1)}
    merged_rows, records, mapping, phases = [], [], [], {}
    shared_counts = {}
    for i,data in enumerate(datasets):
        ids = {key:f'g{gene_rank[i,key]}' for key in data['genes']}
        def roots(row):
            if 'ID' in row['attrs']:
                return data['owners'][unquote(row['attrs']['ID'])]
            return set().union(*(data['owners'][p] for p in references(row)))
        def rank(row):
            return min(gene_rank[i,g] for g in roots(row))
        transcripts = sorted(data['transcripts'], key=lambda t:(rank(data['transcripts'][t]),
            location(data['transcripts'][t]),t))
        counters = {}
        transcript_rank = {}
        for t in transcripts:
            gene = rank(data['transcripts'][t])
            counters[gene] = counters.get(gene,0)+1
            ids[t] = f'g{gene}.t{counters[gene]}'
            transcript_rank[t] = (gene,counters[gene])
            original, (_,dna), phase = data['records'][t]
            records.append((transcript_rank[t],ids[t],dna))
            mapping.append((transcript_rank[t],ids[t],data['tool']+':'+original))
            phases[ids[t]] = phase

        # CDS segments with the same Parent set form one discontinuous feature.
        groups = {}
        row_ids = {}
        for row in data['rows']:
            old = unquote(row['attrs']['ID']) if 'ID' in row['attrs'] else None
            if old in ids:
                row_ids[row['index']] = ids[old]
                continue
            key = ('CDS',tuple(sorted(references(row))))
            groups.setdefault(key,[]).append(row)
        for parents,parts in sorted((key[1],parts) for key,parts in groups.items()):
            if len(parents)>1:
                shared_counts['cds'] = shared_counts.get('cds',0)+1
                new = f"shared_cds{shared_counts['cds']}"
            else:
                new = ids[parents[0]]+'.cds'
            for row in parts:
                row_ids[row['index']] = new
        for row in data['rows']:
            attrs = {'ID':row_ids[row['index']]}
            parents = references(row)
            if parents:
                if any(parent not in ids for parent in parents):
                    raise ValueError('Parent points outside the selected model set')
                attrs['Parent'] = ','.join(sorted({ids[parent] for parent in parents}))
            fields = row['fields'][:]
            if fields[2] == 'transcript':
                fields[2] = 'mRNA'
            fields[8] = ';'.join(f'{key}={value}' for key,value in attrs.items())
            priority = 0 if fields[2]=='gene' else 1 if fields[2] in {'mRNA','transcript'} else 2
            merged_rows.append(((rank(row),priority,location(row),fields[2],attrs['ID']),fields))
    return [fields for _,fields in sorted(merged_rows)], sorted(records), sorted(mapping), phases


def combine_models(sources, output, confidence_threshold=0.1, code=1, record_order=None):
    """Select complete models, unify IDs and translate their native CDS."""
    if not math.isfinite(confidence_threshold) or not 0 <= confidence_threshold <= 1:
        raise ValueError('confidence-threshold must be between 0 and 1')
    datasets = [native_models(tool,path,confidence_threshold) for tool,path in sorted(sources)]
    rows,records,mapping,phases = rename_models(datasets,record_order)
    with (output/'prediction.gff3').open('w') as handle:
        handle.write('##gff-version 3\n')
        for fields in rows:
            handle.write('\t'.join(fields)+'\n')
    read_gff(output/'prediction.gff3')
    with (output/'cds.fa').open('w') as handle:
        for _,name,dna in records:
            handle.write(f'>{name}\n{dna}\n')
    with (output/'id_map.tsv').open('w') as handle:
        writer = csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['new_id','original_id'])
        writer.writerows((name,original) for _,name,original in mapping)
    translation = translate_cds(output/'cds.fa',output/'proteins.fa',output/'translation.tsv',code,phases)
    return dict(genes=sum(len(d['genes']) for d in datasets),models=len(records),translation=translation,
                cds_features=sum(r[2]=='CDS' for r in rows),
                auxiliary_rows_omitted=sum(d['auxiliary_rows_omitted'] for d in datasets),
                carbon_a_incomplete_excluded=sum(d['incomplete_excluded'] for d in datasets),
                carbon_a_confidence_excluded=sum(d['confidence_excluded'] for d in datasets))


def extra_args(text, reserved):
    args = shlex.split(text)
    if any(x.split('=')[0] in reserved for x in args):
        raise ValueError('Input/output/model options are managed by CTC; use its named options')
    return args


def run(a):
    started = time.monotonic()
    table = codon_table(a.codon_table)
    from Bio import __version__ as biopython_version
    output = a.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for tool in ['tiberius', 'carbon_a']:
        (output/tool).mkdir()
    (output/'combined').mkdir()
    report = dict(strategy='CTC', version='0.8.0', input=str(a.input.resolve()), min_length=a.min_length,
                  cutoff_length=a.cutoff_length,
                  confidence_threshold=a.confidence_threshold, translation_engine=f'Biopython {biopython_version}',
                  status='running', components=[])
    try:
        effective_table = output/'combined/codon_table.json'
        effective_table.write_text(json.dumps(table, indent=2, sort_keys=True)+'\n')
        report['codon_table'] = dict(source=f'NCBI {a.codon_table}',
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
                    '--codingseq', str(dest/'cds.fa'), '--protseq', str(dest/'proteins.fa'),
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
            if tool == 'carbon_a':
                for source_name, target in [('prediction.gff','prediction.gff3'),('prediction.fna','cds.fa'),('prediction.faa','proteins.fa')]:
                    shutil.move(str(gff.parent/source_name), str(dest/target))
            sources.append((tool, dest))
        merge_start = time.monotonic()
        with (output/'combined/routing.tsv').open() as handle:
            record_order = {r['sequence']:i for i,r in enumerate(csv.DictReader(handle,delimiter='\t'))}
        report['combined'] = combine_models(sources, output/'combined', a.confidence_threshold, a.codon_table, record_order)
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
    p.add_argument('--confidence-threshold', type=float, default=0.1,
                   help='Minimum Carbon-A complete-model confidence for combined outputs (inclusive; default: 0.1); 0 disables score filtering, incomplete models are always excluded')
    p.add_argument('--version', action='version', version='CTC 0.8.0')
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
    p.add_argument('--codon-table', type=int, default=1, help='NCBI code for combined proteins and native Carbon-A (default: 1); Tiberius native proteins use its own default; Carbon-A cannot represent codes 2, 22, 23, 27, 28 or 31')
    a = p.parse_args()
    if not 1 <= a.min_length <= a.cutoff_length:
        p.error('Require 1 <= min-length <= cutoff-length')
    if not math.isfinite(a.confidence_threshold) or not 0 <= a.confidence_threshold <= 1:
        p.error('confidence-threshold must be between 0 and 1')
    try:
        result = run(a)
    except (ValueError, OSError, RuntimeError, ImportError) as error:
        p.exit(1, f'CTC: {error}\n')
    print(f"CTC complete: {a.output}/combined/prediction.gff3 ({result['elapsed_seconds']:.1f} s)")


if __name__ == '__main__':
    main()
