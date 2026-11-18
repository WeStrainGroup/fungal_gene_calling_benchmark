#!/usr/bin/env python3
"""Independently reconcile v0.2 gene counts and all TP coordinates from raw GFF."""
import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote


def table(path):
    return list(csv.DictReader(path.open(), delimiter='\t'))


def union(intervals):
    result = []
    for a, b in sorted(intervals):
        if result and a <= result[-1][1]+1:
            result[-1] = (result[-1][0], max(b, result[-1][1]))
        else:
            result.append((a, b))
    return result


def raw_gff(path):
    genes, features, stops = {}, collections.defaultdict(list), collections.defaultdict(list)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip().split('\t')
            a = {k: unquote(v) for k, v in (x.split('=', 1) for x in c[8].split(';') if '=' in x)}
            if c[2] in ('mRNA', 'transcript'):
                genes[a['ID']] = a.get('Parent', a['ID']).split(',')[0]
            if c[2] in ('CDS', 'stop_codon'):
                for parent in a['Parent'].split(','):
                    (features if c[2] == 'CDS' else stops)[parent].append((c[0], c[6], int(c[3]), int(c[4]), c[7]))
    result = {}
    for tx, rows in features.items():
        ordered = sorted(rows, key=lambda r: r[2], reverse=rows[0][1] == '-')
        phase = 0 if ordered[0][4] == '.' else int(ordered[0][4])
        expected, consistent = phase, True
        for row in ordered:
            consistent &= row[4] != '.' and int(row[4]) == expected
            expected = (expected - (row[3]-row[2]+1)) % 3
        intervals = union([(r[2], r[3]) for r in rows+stops[tx]])
        result[tx] = dict(record=rows[0][0], strand=rows[0][1], blocks=intervals,
                          phase=phase, phase_consistent=consistent, gene=genes.get(tx, tx))
    return result


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def run(data, output):
    metrics = table(output/'metrics.tsv')
    case_rows = table(output/'validation.tsv')
    region_rows = table(output/'reference_regions.tsv')
    regions = collections.defaultdict(list)
    for r in region_rows:
        regions[(r['species'], r['group'])].append(r)
    reference_cache, examples, checked_tp = {}, [], 0
    for case in case_rows:
        species, group, tool = (case[k] for k in ('species', 'group', 'tool'))
        if species not in reference_cache:
            source = next((data/'provenance'/species).glob('*_genomic.gff.gz'))
            raw = raw_gff(source)
            passed = {r['transcript'] for r in table(output/'reference'/species/'qualification.tsv') if r['status']=='PASS'}
            by_gene = collections.defaultdict(list)
            for tx, m in raw.items():
                if tx in passed:
                    by_gene[m['gene']].append(m)
            reference_cache[species] = by_gene
        mapping = {} if group == 'genome' else {r['record']: r for r in table(data/'inputs'/species/f'{group}.fragments.tsv')}
        native = data/'out'/tool/species/group/('prediction.gff' if tool=='carbon' else 'prediction.gff3')
        predicted = raw_gff(native)
        audit = table(output/'cases'/species/group/tool/'prediction_audit.tsv')
        assert len(predicted) == len(audit) == int(case['predictions'])
        assert set(predicted) == {r['transcript'] for r in audit}
        matched_keys = []
        for row in audit:
            if not row['matched_reference']:
                continue
            key = tuple(row['matched_reference'].split('|'))
            matched_keys.append(key)
            assert row['phase_consistent'] == 'True'
            pred = predicted[row['transcript']]
            assert pred['phase_consistent'], 'Inconsistent native exon phase'
            candidates = []
            for ref in reference_cache[species][key[2]]:
                if group == 'genome':
                    if ref['record'] == key[0]:
                        candidates.append((ref['strand'], ref['blocks'], ref['phase']))
                    continue
                fragment = mapping[key[0]]
                if ref['record'] != fragment['chrom']:
                    continue
                start, end = int(fragment['start']), int(fragment['end'])
                blocks = [(max(a, start)-start+1, min(b, end)-start+1)
                          for a, b in ref['blocks'] if a <= end and b >= start]
                if not blocks:
                    continue
                removed = 0
                for a, b in ref['blocks']:
                    if ref['strand'] == '+':
                        removed += len(range(a, min(b+1, start)))
                    else:
                        removed += len(range(max(a, end+1), b+1))
                candidates.append((ref['strand'], blocks, (ref['phase']-removed)%3))
            assert pred['record'] == key[0] and pred['strand'] == key[1]
            assert (pred['strand'], pred['blocks'], pred['phase']) in candidates, (species, group, tool, row)
            checked_tp += 1
        assert len(matched_keys) == len(set(matched_keys)), 'Reference counted twice'
        for stratum in ('complete', 'partial'):
            if stratum == 'partial' and tool != 'carbon':
                continue
            expected = next(r for r in metrics if (r['species'],r['group'],r['tool'],r['stratum'],r['level']) ==
                            (species,group,tool,stratum,'gene'))
            refs = {(r['record'],r['strand'],r['gene']) for r in regions[(species,group)] if r['stratum']==stratum}
            pred_rows = [r for r in audit if r['bucket']==stratum]
            tp = sum(bool(r['matched_reference']) for r in pred_rows)
            assert all(tuple(r['matched_reference'].split('|')) in refs for r in pred_rows if r['matched_reference'])
            counts = (tp, len(pred_rows)-tp, len(refs)-tp)
            assert counts == tuple(int(expected[k]) for k in ('tp','fp','fn'))
            if refs:
                precision = tp/len(pred_rows) if pred_rows else 0
                recall = tp/len(refs)
                f1 = 2*precision*recall/(precision+recall) if precision+recall else 0
                for k,v in [('precision',precision),('recall',recall),('f1',f1)]:
                    assert abs(float(expected[k])-100*v)<1e-9
            else:
                assert all(expected[k]=='NA' for k in ('precision','recall','f1'))
        for category in ('complete','partial','cross_scope','unscored_reference'):
            sample = next((r for r in audit if r['bucket']==category and
                           (bool(r['matched_reference']) if category in ('complete','partial') else True)), None)
            if sample:
                p = predicted[sample['transcript']]
                examples.append(dict(species=species,group=group,tool=tool,category=category,
                    transcript=sample['transcript'],reference=sample['matched_reference'],record=p['record'],strand=p['strand'],
                    coordinates=','.join(f'{a}-{b}' for a,b in p['blocks']),phase=p['phase'],confidence=sample['confidence']))
    sources = json.loads((output/'source_sha256.json').read_text())
    for name, value in sources.items():
        assert sha256(data/name)==value, name
    with (output/'validation_examples.tsv').open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(examples[0]) if examples else [],delimiter='\t')
        w.writeheader();w.writerows(examples)
    receipt=dict(status='PASS',cases=len(case_rows),all_gene_true_positives_checked=checked_tp,
                 raw_gff_coordinate_and_phase_agreement='PASS',independent_gene_accounting='PASS',
                 source_files_rehashed=len(sources),source_files_unchanged=True,
                 examples=len(examples),note='Read native GFF and original reference directly; no evaluator imports.')
    (output/'independent_validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.data.resolve(),a.output.resolve())
