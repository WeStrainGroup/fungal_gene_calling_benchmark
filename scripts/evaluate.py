#!/usr/bin/env python3
"""Reevaluate existing annotations with the frozen official Carbon-A evaluator."""
import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / 'data'
sys.path.insert(0, str(REPO / 'vendor/carbon_a'))
import evaluation as official

SPECIES = ('Saccharomyces_cerevisiae', 'Candida_albicans')
TOOLS = ('carbon', 'annevo', 'tiberius')
GROUPS = ('genome', '50-100kb', '5-10kb', '1-2kb')
LEVELS = ('nucleotide', 'exon', 'gene')


def write_table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def fasta(path):
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, 'fasta')}


def read_models(path, include_stops=False):
    transcripts, blocks, stops = {}, collections.defaultdict(list), collections.defaultdict(list)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip().split('\t')
            a = dict(x.split('=', 1) for x in c[8].split(';') if '=' in x)
            if c[2] in ('mRNA', 'transcript'):
                transcripts[a['ID']] = a.get('Parent', a['ID']).split(',')[0]
            if c[2] in ('CDS', 'stop_codon'):
                for parent in a.get('Parent', '').split(','):
                    assert parent
                    block = (c[0], int(c[3]), int(c[4]), c[6], c[7], a.get('partial') == 'true')
                    (blocks if c[2] == 'CDS' else stops)[parent].append(block)
    result = []
    for tx, parts in blocks.items():
        parts.sort(key=lambda b: b[1])
        assert len({(b[0], b[3]) for b in parts}) == 1
        first = parts[0] if parts[0][3] == '+' else parts[-1]
        if include_stops: assert first[4] in ('0', '1', '2')
        intervals = [(b[1], b[2]) for b in parts]
        if include_stops:
            assert all((b[0], b[3]) == (first[0], first[3]) for b in stops[tx])
            intervals = official.merge_intervals(intervals + [(b[1], b[2]) for b in stops[tx]])
        result.append(dict(id=tx, gene=transcripts.get(tx, tx), record=first[0], strand=first[3],
                           blocks=intervals, phase=int(first[4]) if first[4] != '.' else 0,
                           partial=any(b[5] for b in parts)))
    if include_stops:
        assert len(result) == len(transcripts), (path, len(result), len(transcripts))
    return result


def sequence(model, genome):
    assert all(1 <= a <= b <= len(genome[model['record']]) for a, b in model['blocks'])
    seq = ''.join(genome[model['record']][a-1:b] for a, b in model['blocks'])
    return str(Seq(seq).reverse_complement()) if model['strand'] == '-' else seq


def qualified_reference(species, genome):
    source = next((ROOT / 'provenance' / species).glob('*_genomic.gff.gz'))
    table = json.loads((ROOT / 'inputs' / species / 'codon_table.json').read_text())
    kept, audit = [], []
    for m in read_models(source):
        reasons = []
        if m['record'] not in genome:
            reasons.append('excluded_record')
            seq = ''
        else:
            seq = sequence(m, genome)
            if m['phase'] != 0: reasons.append('nonzero_initial_phase')
            if m['partial']: reasons.append('partial_annotation')
            b = m['blocks']
            if any(b[i+1][0] - b[i][1] - 1 < 2 for i in range(len(b)-1)):
                reasons.append('invalid_intron')
            if not seq.startswith('ATG'): reasons.append('missing_start')
            if len(seq) % 3: reasons.append('length_not_multiple_of_3')
            if table.get(seq[-3:]) != '*': reasons.append('missing_stop')
            aa = ''.join(table.get(seq[i:i+3], 'X') for i in range(0, len(seq)-2, 3))
            if '*' in aa[:-1]: reasons.append('internal_stop')
        audit.append(dict(transcript=m['id'], gene=m['gene'], record=m['record'],
                          cds_bp=len(seq), status=','.join(reasons) or 'PASS'))
        if not reasons: kept.append(m)
    dest = ROOT / 'reference' / species
    dest.mkdir(parents=True, exist_ok=True)
    write_table(dest / 'qualification.tsv', audit)
    return kept, collections.Counter(r['status'] for r in audit)


def clipped_reference(models, species, group, original, genome):
    fragments = list(csv.DictReader((ROOT / 'inputs' / species / f'{group}.fragments.tsv').open(), delimiter='\t'))
    by_chrom = collections.defaultdict(list)
    for r in fragments:
        r['start'], r['end'] = int(r['start']), int(r['end'])
        assert genome[r['record']] == original[r['chrom']][r['start']-1:r['end']]
        assert hashlib.sha256(genome[r['record']].encode()).hexdigest() == r['sequence_sha256']
        by_chrom[r['chrom']].append(r)
    assert len(fragments) == len(genome)
    for chrom, rows in by_chrom.items():
        assert rows[0]['start'] == 1 and rows[-1]['end'] == len(original[chrom])
        assert all(a['end'] + 1 == b['start'] for a, b in zip(rows, rows[1:]))
    assert set(by_chrom) == set(original)
    ends = {chrom: [r['end'] for r in rows] for chrom, rows in by_chrom.items()}
    clipped = []
    for m in models:
        rows = by_chrom[m['record']]
        start_index = bisect.bisect_left(ends[m['record']], m['blocks'][0][0])
        reconstructed = []
        for r in rows[start_index:]:
            start, end = r['start'], r['end']
            if start > m['blocks'][-1][1]: break
            blocks = [(max(a, start)-start+1, min(b, end)-start+1)
                      for a, b in m['blocks'] if a <= end and b >= start]
            if not blocks: continue
            removed = sum(max(0, min(b, start-1)-a+1) if m['strand'] == '+'
                          else max(0, b-max(a, end+1)+1) for a, b in m['blocks'])
            cm = dict(m, id=r['record']+'__'+m['id'], record=r['record'], blocks=blocks,
                      phase=(m['phase']-removed) % 3,
                      partial=not (m['blocks'][0][0] >= start and m['blocks'][-1][1] <= end))
            clipped.append(cm)
            reconstructed.append(sequence(cm, genome))
        if m['strand'] == '-': reconstructed.reverse()
        assert ''.join(reconstructed) == sequence(m, original), m['id']
    return clipped


def export_models(models, genome, prefix, reference=False):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    gff_path = prefix.parent / 'genomic.gff' if reference else prefix.with_suffix('.gff')
    fa_path = prefix.parent / 'cds_from_genomic.fna' if reference else prefix.with_suffix('.fna')
    with gff_path.open('w') as gf, fa_path.open('w') as ff:
        gf.write('##gff-version 3\n')
        for m in models:
            record, tx, gene, strand = (m[k] for k in ('record', 'id', 'gene', 'strand'))
            gf.write(f'{record}\tbenchmark\tmRNA\t{m["blocks"][0][0]}\t{m["blocks"][-1][1]}\t.\t{strand}\t.\tID={tx};Parent={gene}\n')
            ordered = m['blocks'] if strand == '+' else list(reversed(m['blocks']))
            phase = m['phase']
            for a, b in ordered:
                gf.write(f'{record}\tbenchmark\tCDS\t{a}\t{b}\t.\t{strand}\t{phase}\tParent={tx}\n')
                phase = (phase-(b-a+1)) % 3
            if reference:
                header = f'lcl|{record}_cds_{tx} [locus_tag={gene}] [frame={m["phase"]+1}]'
            else:
                header = f'{tx} record_id={record} phase={m["phase"]}'
            ff.write(f'>{header}\n{sequence(m, genome)}\n')


def independent_matching(predictions, refs):
    keys = list(refs)
    candidates = collections.defaultdict(list)
    for j, key in enumerate(keys):
        for seq in set(refs[key]): candidates[(key[0], seq)].append(j)
    rr, cc = [], []
    for i, pred in enumerate(predictions):
        for j in candidates[pred]: rr.append(i); cc.append(j)
    matrix = csr_matrix(([1]*len(rr), (rr, cc)), shape=(len(predictions), len(keys)))
    return int((maximum_bipartite_matching(matrix, perm_type='column') >= 0).sum())


def validate_evaluator():
    refs = {('r', 'a'): ['AAA', 'CCC'], ('r', 'b'): ['AAA']}
    pred = [('r', 'AAA'), ('r', 'CCC')]
    assert len(official.gene_exact_match_assignments(pred, refs)) == 2
    assert independent_matching(pred, refs) == 2
    assert len(official.gene_exact_match_assignments([('other', 'AAA')], refs)) == 0
    assert len(official.gene_exact_match_assignments([('r', 'CCC')]*2, refs)) == 1
    assert official.merge_intervals([(0, 2), (3, 5), (8, 9), (8, 10)]) == [(0, 5), (8, 10)]
    counts = official.interval_metric_counts({('r', '+'): [(0, 2)]}, {('r', '-'): [(0, 2)]})
    assert counts['nucleotide_tp'] == 0 and counts['nucleotide_fp'] == 3


def main():
    validate_evaluator()
    version = json.loads((REPO / 'metadata/carbon_evaluator.json').read_text())
    for name, digest in version['sha256'].items():
        assert hashlib.sha256((REPO / 'vendor/carbon_a' / name).read_bytes()).hexdigest() == digest
    results, summaries, audits = [], [], []
    for species in SPECIES:
        original = fasta(ROOT / 'inputs' / species / 'genome.fna')
        source_dir = ROOT / 'provenance' / species
        checks = {Path(line.split(None, 1)[1]).name: line.split()[0] for line in (source_dir / 'md5checksums.txt').read_text().splitlines()}
        for path in source_dir.glob('*_genomic.*.gz'):
            assert hashlib.md5(path.read_bytes()).hexdigest() == checks[path.name]
        with gzip.open(next(source_dir.glob('*_genomic.fna.gz')), 'rt') as handle:
            source_genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, 'fasta')}
        assert all(seq == source_genome[record] for record, seq in original.items())
        refs, reasons = qualified_reference(species, original)
        print(species, 'qualified models', len(refs), dict(reasons), flush=True)
        for group in GROUPS:
            genome = fasta(ROOT / 'inputs' / species / f'{group}.fna')
            current = refs if group == 'genome' else clipped_reference(refs, species, group, original, genome)
            ref_dir = ROOT / 'reference' / species / group
            export_models(current, genome, ref_dir / 'reference', reference=True)
            ref_tracks = official.parse_reference_cds_tracks(ref_dir / 'genomic.gff')
            ref_genes = official.parse_reference_cds_fasta(ref_dir / 'cds_from_genomic.fna')
            assert len(ref_genes) == len({(m['record'], m['gene']) for m in current})
            assert sum(map(len, ref_genes.values())) == len(current)
            summaries.append(dict(species=species, group=group, records=len(genome), bp=sum(map(len, genome.values())),
                                  reference_genes=len(ref_genes), reference_isoforms=len(current),
                                  partial_isoforms=sum(m['partial'] for m in current)))
            for tool in TOOLS:
                out = ROOT / 'out' / tool / species / group
                source = out / ('prediction.gff' if tool == 'carbon' else 'prediction.gff3')
                models = read_models(source, include_stops=True)
                assert len({(m['record'], m['gene']) for m in models}) == len(models)
                prefix = out / 'evaluation/prediction'
                export_models(models, genome, prefix)
                preds = official.parse_predicted_cds_fasta(prefix.with_suffix('.fna'))
                tracks = official.parse_predicted_cds_tracks(prefix.with_suffix('.gff'))
                counts = official.interval_metric_counts(tracks, ref_tracks)
                assignments = official.gene_exact_match_assignments(preds, ref_genes)
                counts.update(official.gene_metric_counts(preds, ref_genes, assignments))
                assert len(assignments) == independent_matching(preds, ref_genes)
                full = official.evaluate_prediction(prefix, ref_dir, confidence_threshold=0,
                         codon_table_path=ROOT / 'inputs' / species / 'codon_table.json')
                for level in LEVELS:
                    tp, fp, fn = (counts[f'{level}_{suffix}'] for suffix in ('tp', 'fp', 'fn'))
                    assert all(full[f'{level}_{s}'] == counts[f'{level}_{s}'] for s in ('tp', 'fp', 'fn'))
                    p, r, f = official.precision_recall_f1(tp, fp, fn)
                    results.append(dict(species=species, group=group, tool=tool, level=level,
                                        tp=tp, fp=fp, fn=fn, predicted=tp+fp, reference=tp+fn,
                                        precision=100*p, recall=100*r, f1=100*f))
                audits.append(dict(species=species, group=group, tool=tool, raw_models=len(models),
                                   gene_matches=len(assignments), scipy_matching='PASS', official_entrypoint='PASS',
                                   annotation_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
                print(species, group, tool, len(models), len(assignments), flush=True)
    assert len(results) == len(SPECIES)*36 and len(audits) == len(SPECIES)*12
    write_table(ROOT / 'benchmark.tsv', results)
    write_table(ROOT / 'reference_summary.tsv', summaries)
    write_table(ROOT / 'evaluation_validation.tsv', audits)
    (ROOT / 'evaluation_validation.json').write_text(json.dumps(dict(status='PASS', cases=len(audits),
        evaluator_fixtures='PASS', fragment_sequence_reconstruction='PASS',
        official_entrypoint_agreement=f'{len(audits)}/{len(audits)}', independent_maximum_matching=f'{len(audits)}/{len(audits)}',
        external_min_cds_bp=0, confidence_threshold=0), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT, help='Benchmark data directory')
    parser.add_argument('--samples', type=Path, help='Optional TSV with a species column for a new set')
    parser.add_argument('--expected', type=Path, help='Optional published table to check')
    args = parser.parse_args()
    ROOT = args.data.resolve()
    if args.samples:
        SPECIES = tuple(r['species'] for r in csv.DictReader(args.samples.open(), delimiter='\t'))
    main()
    expected_path = args.expected or (REPO / 'results/benchmark.tsv' if not args.samples else None)
    if expected_path:
        expected = list(csv.DictReader(expected_path.open(), delimiter='\t'))
        actual = list(csv.DictReader((ROOT / 'benchmark.tsv').open(), delimiter='\t'))
        assert expected == actual, 'Recomputed metrics differ from the published table'
        print(f'All {len(actual)} published metric rows reproduced exactly.')
