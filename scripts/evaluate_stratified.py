#!/usr/bin/env python3
"""v0.2: truth-defined CDS strata and exact coordinate/frame matching."""
import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import unquote

import common
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

VERSION = '0.2.0'
GROUPS = common.GROUPS
TOOLS = common.TOOLS
SUBTYPES = ('missing_5prime', 'missing_3prime', 'missing_both')
LEVELS = ('nucleotide', 'exon', 'gene')
MIN_RECORD = 1000
MIN_CDS = 200
CONFIDENCE = 0.1
AUDIT_FIELDS = ('transcript', 'gene', 'record', 'strand', 'cds_bp', 'phase',
                'phase_consistent', 'confidence', 'native_status', 'incomplete_reason',
                'bucket', 'background', 'matched_reference', 'subtype_unresolved')


def write_table(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def length(model):
    return sum(b - a + 1 for a, b in model['blocks'])


def span(model):
    return model['blocks'][0][0], model['blocks'][-1][1]


def gene_key(model):
    return model['record'], model['strand'], model['gene']


def structural_key(model):
    if not model.get('phase_consistent', True):
        return None
    return model['record'], model['strand'], tuple(model['blocks']), model['phase']


def read_models(path):
    """Preserve gene confidence/status and transcript/CDS phase metadata."""
    attrs, parents = {}, {}
    blocks, stops = collections.defaultdict(list), collections.defaultdict(list)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            c = line.rstrip().split('\t')
            a = dict(x.split('=', 1) for x in c[8].split(';') if '=' in x)
            a = {k: unquote(v) for k, v in a.items()}
            if c[2] in ('gene', 'mRNA', 'transcript'):
                attrs[a['ID']] = a
                if c[2] != 'gene':
                    parents[a['ID']] = a.get('Parent', a['ID']).split(',')[0]
            if c[2] in ('CDS', 'stop_codon'):
                for parent in a.get('Parent', '').split(','):
                    if not parent:
                        raise ValueError(f'CDS without Parent in {path}')
                    row = (c[0], int(c[3]), int(c[4]), c[6], c[7], a.get('partial') == 'true')
                    (blocks if c[2] == 'CDS' else stops)[parent].append(row)
    result = []
    for tx, parts in blocks.items():
        parts.sort(key=lambda b: b[1])
        assert len({(b[0], b[3]) for b in parts + stops[tx]}) == 1
        ordered = parts if parts[0][3] == '+' else parts[::-1]
        phase = int(ordered[0][4]) if ordered[0][4] != '.' else 0
        expected, consistent = phase, True
        for row in ordered:
            consistent &= row[4] != '.' and int(row[4]) == expected
            expected = (expected - (row[2] - row[1] + 1)) % 3
        gene = parents.get(tx, tx)
        attributes = dict(attrs.get(gene, {}), **attrs.get(tx, {}))
        value = attributes.get('confidence')
        confidence = None if value in (None, 'NA') else float(value)
        if confidence is not None and (not math.isfinite(confidence) or not 0 <= confidence <= 1):
            raise ValueError(f'Invalid confidence: {path}: {tx}: {value}')
        intervals = common.official.merge_intervals([(r[1], r[2]) for r in parts + stops[tx]])
        result.append(dict(id=tx, gene=gene, record=parts[0][0], strand=parts[0][3],
                           blocks=intervals, raw_blocks=[(r[1], r[2]) for r in parts], phase=phase, phase_consistent=consistent,
                           partial=any(r[5] for r in parts), confidence=confidence,
                           confidence_present=value is not None, status=attributes.get('status', ''),
                           incomplete_reason=attributes.get('incomplete_reason', '')))
    return result


def reference_reason(model, genome, code):
    """Retain the original benchmark's sequence qualification rules."""
    if model['record'] not in genome:
        return 'excluded_record'
    seq = common.sequence(model, genome)
    reasons = []
    if model['phase'] != 0:
        reasons.append('nonzero_initial_phase')
    if model['partial']:
        reasons.append('partial_annotation')
    raw_blocks = model.get('raw_blocks', model['blocks'])
    if any(c - b - 1 < 2 for (_, b), (c, _) in zip(raw_blocks, raw_blocks[1:])):
        reasons.append('invalid_intron')
    if not seq.startswith('ATG'):
        reasons.append('missing_start')
    if len(seq) % 3:
        reasons.append('length_not_multiple_of_3')
    if code.get(seq[-3:]) != '*':
        reasons.append('missing_stop')
    aa = ''.join(code.get(seq[i:i+3], 'X') for i in range(0, len(seq)-2, 3))
    if '*' in aa[:-1]:
        reasons.append('internal_stop')
    return ','.join(reasons) or 'PASS'


def clip_model(model, record, start, end):
    blocks = [(max(a, start)-start+1, min(b, end)-start+1)
              for a, b in model['blocks'] if a <= end and b >= start]
    if not blocks:
        return None
    left = sum(max(0, min(b, start-1)-a+1) for a, b in model['blocks'])
    right = sum(max(0, b-max(a, end+1)+1) for a, b in model['blocks'])
    five, three = (left, right) if model['strand'] == '+' else (right, left)
    kind = ('missing_both' if five and three else 'missing_5prime' if five
            else 'missing_3prime' if three else 'complete')
    return dict(model, record=record, blocks=blocks, phase=(model['phase']-five) % 3,
                phase_consistent=True, kind=kind, removed_5prime=five, removed_3prime=three)


def project_reference(models, original, genome, fragments):
    """Validate entire partitions and reconstruct every original CDS."""
    by_chrom = collections.defaultdict(list)
    for row in fragments:
        start, end = int(row['start']), int(row['end'])
        seq = genome[row['record']]
        assert seq == original[row['chrom']][start-1:end]
        assert hashlib.sha256(seq.encode()).hexdigest() == row['sequence_sha256']
        by_chrom[row['chrom']].append((start, end, row['record']))
    assert set(by_chrom) == set(original)
    assert {r[2] for rows in by_chrom.values() for r in rows} == set(genome)
    ends = {}
    for chrom, rows in by_chrom.items():
        assert rows[0][0] == 1 and rows[-1][1] == len(original[chrom])
        assert all(a[1]+1 == b[0] for a, b in zip(rows, rows[1:]))
        ends[chrom] = [r[1] for r in rows]
    projected = []
    for model in models:
        if model['record'] not in original:
            continue
        rows = by_chrom[model['record']]
        first = bisect.bisect_left(ends[model['record']], span(model)[0])
        visible = []
        for start, end, record in rows[first:]:
            if start > span(model)[1]:
                break
            clipped = clip_model(model, record, start, end)
            if clipped:
                projected.append(clipped)
                visible.append(common.sequence(clipped, genome))
        if model['strand'] == '-':
            visible.reverse()
        assert ''.join(visible) == common.sequence(model, original), model['id']
    return projected


class Intervals:
    """Sorted disjoint intervals for logarithmic overlap/containment queries."""
    def __init__(self, intervals=()):
        self.rows = common.official.merge_intervals(list(intervals))
        self.starts = [a for a, _ in self.rows]
        self.ends = [b for _, b in self.rows]

    def overlaps(self, interval):
        a, b = interval
        i = bisect.bisect_left(self.ends, a)
        return i < len(self.rows) and self.starts[i] <= b

    def contains(self, interval):
        a, b = interval
        i = bisect.bisect_right(self.starts, a)-1
        return i >= 0 and self.ends[i] >= b


def build_domains(models, genome):
    """Resolve same-strand overlap components before assigning scoring domains."""
    by_gene = collections.defaultdict(list)
    for model in models:
        if len(genome[model['record']]) >= MIN_RECORD:
            by_gene[gene_key(model)].append(model)
    units = []
    for key, items in sorted(by_gene.items()):
        qualified = [m for m in items if m['qualification'] == 'PASS']
        classes = {'complete' if m['kind'] == 'complete' else 'partial' for m in qualified}
        eligible = [m for m in qualified if length(m) >= MIN_CDS]
        reason = ('unqualified_reference' if not qualified else 'mixed_isoform_class' if len(classes) > 1
                  else 'reference_below_200bp' if not eligible else '')
        chosen = eligible if not reason else items
        subtypes = {m['kind'] for m in chosen}
        unit = dict(key=key, models=chosen, kind='ignore' if reason else next(iter(classes)),
                    subtype=next(iter(subtypes)) if len(subtypes) == 1 else 'mixed_ends',
                    interval=(min(span(m)[0] for m in chosen), max(span(m)[1] for m in chosen)),
                    reason=reason)
        units.append(unit)
        # Failed isoforms remain unknown, including where they overlap an eligible isoform.
        if not reason:
            for model in items:
                if model['qualification'] != 'PASS':
                    units.append(dict(key=key, models=[], kind='ignore', subtype='unknown',
                                      interval=span(model), reason='unqualified_isoform'))
                elif length(model) < MIN_CDS:
                    # Eligible isoforms already account for their shared gene region.
                    # Mask any short-isoform extension rather than calling it background.
                    a, b = span(model)
                    left, right = unit['interval']
                    for start, end in ((a, min(b, left-1)), (max(a, right+1), b)):
                        if start <= end:
                            units.append(dict(key=key, models=[], kind='ignore', subtype='unknown',
                                              interval=(start, end), reason='short_isoform_extension'))
    tracks = collections.defaultdict(list)
    for i, unit in enumerate(units):
        tracks[unit['key'][:2]].append(i)
    ignored = collections.defaultdict(list)
    for track, indices in tracks.items():
        indices.sort(key=lambda i: units[i]['interval'])
        components = []
        for i in indices:
            a, b = units[i]['interval']
            if not components or a > components[-1][1]:
                components.append([a, b, [i]])
            else:
                components[-1][1] = max(b, components[-1][1])
                components[-1][2].append(i)
        for a, b, members in components:
            types = {units[i]['kind'] for i in members}
            if 'ignore' in types or len(types) > 1:
                ignored[track].append((a, b))
                for i in members:
                    if units[i]['kind'] != 'ignore':
                        units[i]['reason'] = 'overlapping_unscored_reference' if 'ignore' in types else 'complete_partial_overlap'
                        units[i]['kind'] = 'ignore'
    domains = {}
    for track, indices in tracks.items():
        local_units = [units[i] for i in indices]
        domains[track] = {'ignore': Intervals(ignored[track])}
        for kind in ('complete', 'partial'):
            domains[track][kind] = Intervals(u['interval'] for u in local_units if u['kind'] == kind)
        for subtype in SUBTYPES + ('mixed_ends',):
            domains[track][subtype] = Intervals(u['interval'] for u in local_units if u['kind'] == 'partial'
                                                and u['subtype'] == subtype)
    return units, domains


def prediction_bucket(model, tool, genome, domains):
    if model['record'] not in genome:
        raise ValueError(f'Prediction on unknown input record: {model["record"]}')
    common.sequence(model, genome)  # Validate every coordinate, even for filtered models.
    if tool == 'carbon' and not model['confidence_present']:
        raise ValueError(f'Carbon-A model lacks an explicit confidence attribute: {model["id"]}')
    if len(genome[model['record']]) < MIN_RECORD:
        return 'input_below_1kb', False
    if length(model) < MIN_CDS:
        return 'cds_below_200bp', False
    if tool == 'carbon' and model['confidence'] is not None and model['confidence'] < CONFIDENCE:
        return 'low_confidence', False
    d = domains.get((model['record'], model['strand']))
    if d is None:
        return 'complete', True
    interval = span(model)
    if d['ignore'].overlaps(interval):
        return 'unscored_reference', False
    if d['partial'].overlaps(interval):
        return ('partial', False) if d['partial'].contains(interval) else ('cross_scope', False)
    return 'complete', not d['complete'].overlaps(interval)


def matching(predictions, references):
    """Maximum one-to-one matching to gene IDs, allowing qualified isoforms."""
    candidates = collections.defaultdict(set)
    for key, models in references.items():
        for model in models:
            value = structural_key(model)
            if value is not None:
                candidates[value].add(key)
    edges = [sorted(candidates.get(structural_key(p), ())) for p in predictions]
    owners = {}
    def augment(i, seen):
        for key in edges[i]:
            if key in seen:
                continue
            seen.add(key)
            if key not in owners or augment(owners[key], seen):
                owners[key] = i
                return True
        return False
    for i in range(len(predictions)):
        augment(i, set())
    result = {i: key for key, i in owners.items()}
    # Independent sparse maximum matching checks the TP count.
    keys = sorted(references)
    index = {key: j for j, key in enumerate(keys)}
    rr, cc = [], []
    for i, choices in enumerate(edges):
        for key in choices:
            rr.append(i); cc.append(index[key])
    if predictions and keys:
        matrix = csr_matrix(([1]*len(rr), (rr, cc)), shape=(len(predictions), len(keys)))
        independent = int((maximum_bipartite_matching(matrix, perm_type='column') >= 0).sum())
    else:
        independent = 0
    assert len(result) == independent
    assert len(set(result.values())) == len(result)
    return result


def tracks(models):
    result = collections.defaultdict(list)
    for model in models:
        result[(model['record'], model['strand'])].extend(model['blocks'])
    return result


def metric_rows(predictions, references):
    matched = matching(predictions, references)
    longest = [max(models, key=lambda m: (length(m), m['id'])) for models in references.values()]
    counts = common.official.interval_metric_counts(tracks(predictions), tracks(longest))
    counts.update(gene_tp=len(matched), gene_fp=len(predictions)-len(matched), gene_fn=len(references)-len(matched))
    rows = []
    for level in LEVELS:
        tp, fp, fn = [counts[f'{level}_{s}'] for s in ('tp', 'fp', 'fn')]
        if not references:
            p = r = f = 'NA'
        else:
            p, r, f = (100*x for x in common.official.precision_recall_f1(tp, fp, fn))
        rows.append(dict(level=level, tp=tp, fp=fp, fn=fn, predicted=tp+fp, reference=tp+fn,
                         precision=p, recall=r, f1=f, assessment='PASS' if references else 'NO_REFERENCE'))
    return rows, matched


def score_case(predictions, tool, genome, units, domains):
    buckets = collections.defaultdict(list)
    audit = []
    for model in predictions:
        bucket, background = prediction_bucket(model, tool, genome, domains)
        buckets[bucket].append(model)
        audit.append(dict(transcript=model['id'], gene=model['gene'], record=model['record'],
                          strand=model['strand'], cds_bp=length(model), phase=model['phase'],
                          phase_consistent=model['phase_consistent'], confidence='NA' if model['confidence'] is None else model['confidence'],
                          native_status=model['status'], incomplete_reason=model['incomplete_reason'],
                          bucket=bucket, background=background, matched_reference=''))
    references = {kind: {u['key']: u['models'] for u in units if u['kind'] == kind}
                  for kind in ('complete', 'partial')}
    results, matches = [], {}
    for kind in ('complete', 'partial'):
        rows, assignments = metric_rows(buckets[kind], references[kind])
        matches[kind] = assignments
        if kind == 'complete' or tool == 'carbon':
            results.extend(dict(stratum=kind, **r) for r in rows)
        for i, key in assignments.items():
            matches[(kind, buckets[kind][i]['id'])] = key
    # Subtype TP follows the globally matched reference; unmatched FP uses location.
    subtype_by_gene = {u['key']: u['subtype'] for u in units if u['kind'] == 'partial'}
    subtype_predictions = collections.defaultdict(list)
    unresolved_subtypes = set()
    if tool == 'carbon':
        for i, model in enumerate(buckets['partial']):
            if i in matches['partial']:
                subtype = subtype_by_gene[matches['partial'][i]]
            else:
                d = domains[(model['record'], model['strand'])]
                hit = [s for s in SUBTYPES + ('mixed_ends',) if d[s].overlaps(span(model))]
                subtype = hit[0] if len(hit) == 1 and d[hit[0]].contains(span(model)) else 'mixed_ends'
            subtype_predictions[subtype].append(model)
            if subtype == 'mixed_ends':
                unresolved_subtypes.add(model['id'])
        for subtype in SUBTYPES:
            refs = {key: models for key, models in references['partial'].items() if subtype_by_gene[key] == subtype}
            rows, _ = metric_rows(subtype_predictions[subtype], refs)
            results.extend(dict(stratum=subtype, **r) for r in rows)
    for row in audit:
        key = matches.get((row['bucket'], row['transcript']))
        if key:
            row['matched_reference'] = '|'.join(key)
        row['subtype_unresolved'] = row['transcript'] in unresolved_subtypes
    stats = []
    for bucket in ('input_below_1kb', 'cds_below_200bp', 'low_confidence', 'unscored_reference', 'cross_scope', 'complete', 'partial'):
        models = buckets[bucket]
        stats.append(dict(category=bucket, models=len(models), cds_bp=sum(map(length, models)),
                          percent_of_raw=100*len(models)/len(predictions) if predictions else 0))
    for name, rows in [('background_fp', [r for r in audit if r['background']]),
                       ('subtype_unresolved', [r for r in audit if r['subtype_unresolved']])]:
        assert name != 'background_fp' or all(not r['matched_reference'] for r in rows)
        stats.append(dict(category=name, models=len(rows), cds_bp=sum(r['cds_bp'] for r in rows),
                          percent_of_raw=100*len(rows)/len(predictions) if predictions else 0))
    assert sum(len(buckets[k]) for k in buckets) == len(predictions)
    assert sum(bool(r['matched_reference']) for r in audit) == len(matches['complete']) + len(matches['partial'])
    return results, stats, audit


def validate_metric_grid(rows, cases):
    expected = {(c['species'], c['group'], c['tool'], scope, level)
                for c in cases
                for scope in (('complete', 'partial') + SUBTYPES if c['tool'] == 'carbon' else ('complete',))
                for level in LEVELS}
    actual = [(r['species'], r['group'], r['tool'], r['stratum'], r['level']) for r in rows]
    assert len(actual) == len(set(actual)), 'Duplicate metric row'
    assert set(actual) == expected, 'Incomplete metric grid'


def run(args):
    data, output = args.data.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.samples:
        species = [r['species'] for r in csv.DictReader(args.samples.open(), delimiter='\t')]
    else:
        species = list(common.SPECIES)
    source_hashes, metrics, diagnostics, input_rows, ref_rows, pending, validations = {}, [], [], [], [], [], []
    for name in species:
        folder = data/'inputs'/name
        original = common.fasta(folder/'genome.fna')
        code = json.loads((folder/'codon_table.json').read_text())
        source = next((data/'provenance'/name).glob('*_genomic.gff.gz'))
        source_hashes[str(source.relative_to(data))] = digest(source)
        fasta_source = next((data/'provenance'/name).glob('*_genomic.fna.gz'))
        checks = {Path(line.split(None, 1)[1]).name: line.split()[0]
                  for line in (source.parent/'md5checksums.txt').read_text().splitlines()}
        for reference_file in (source, fasta_source):
            assert hashlib.md5(reference_file.read_bytes()).hexdigest() == checks[reference_file.name]
            source_hashes[str(reference_file.relative_to(data))] = digest(reference_file)
        with gzip.open(fasta_source, 'rt') as handle:
            source_genome = {r.id: str(r.seq).upper() for r in common.SeqIO.parse(handle, 'fasta')}
        assert all(seq == source_genome[record] for record, seq in original.items())
        source_hashes[str((folder/'codon_table.json').relative_to(data))] = digest(folder/'codon_table.json')
        models = read_models(source)
        qualification = []
        for m in models:
            m['qualification'] = reference_reason(m, original, code)
            qualification.append(dict(transcript=m['id'], gene=m['gene'], record=m['record'], status=m['qualification']))
        write_table(output/'reference'/name/'qualification.tsv', qualification)
        for group in args.groups:
            path = folder/f'{group}.fna'
            genome = common.fasta(path)
            source_hashes[str(path.relative_to(data))] = digest(path)
            if group == 'genome':
                current = [dict(m, kind='complete', phase_consistent=True, removed_5prime=0, removed_3prime=0)
                           for m in models if m['record'] in genome]
            else:
                mapping = folder/f'{group}.fragments.tsv'
                source_hashes[str(mapping.relative_to(data))] = digest(mapping)
                current = project_reference(models, original, genome, list(csv.DictReader(mapping.open(), delimiter='\t')))
            eligible = {k: v for k, v in genome.items() if len(v) >= MIN_RECORD}
            input_rows.append(dict(species=name, group=group, total_records=len(genome), total_bp=sum(map(len, genome.values())),
                                   scored_records=len(eligible), scored_bp=sum(map(len, eligible.values())),
                                   excluded_records=len(genome)-len(eligible), excluded_bp=sum(len(v) for v in genome.values() if len(v)<MIN_RECORD)))
            units, domains = build_domains(current, genome)
            group_ref_rows = []
            for unit in units:
                row = dict(species=name, group=group, record=unit['key'][0], strand=unit['key'][1], gene=unit['key'][2],
                           stratum=unit['kind'], subtype=unit['subtype'], start=unit['interval'][0], end=unit['interval'][1],
                           isoforms=len(unit['models']), reason=unit['reason'])
                group_ref_rows.append(row)
            ref_rows.extend(group_ref_rows)
            write_table(output/'reference'/name/group/'regions.tsv', group_ref_rows)
            for tool in TOOLS:
                native = data/'out'/tool/name/group/('prediction.gff' if tool == 'carbon' else 'prediction.gff3')
                if not native.is_file() or (args.completed_only and not (native.parent/'COMPLETE').exists()):
                    pending.append(dict(species=name, group=group, tool=tool))
                    continue
                source_hashes[str(native.relative_to(data))] = digest(native)
                predictions = read_models(native)
                rows, stats, audit = score_case(predictions, tool, genome, units, domains)
                prefix = dict(cohort=args.cohort, species=name, group=group, tool=tool)
                metrics.extend(dict(**prefix, **r) for r in rows)
                diagnostics.extend(dict(**prefix, **r) for r in stats)
                write_table(output/'cases'/name/group/tool/'prediction_audit.tsv', audit, AUDIT_FIELDS)
                assert digest(native) == source_hashes[str(native.relative_to(data))], 'Prediction changed while being read'
                validations.append(dict(**prefix, predictions=len(predictions), matching='SCIPY_AGREES',
                                        coordinate_accounting='PASS', source_sha256=source_hashes[str(native.relative_to(data))]))
                print(name, group, tool, 'PASS', len(predictions), flush=True)
    expected_hashes_path = data/'input_sha256.json'
    if expected_hashes_path.exists():
        expected_hashes = json.loads(expected_hashes_path.read_text())
        for path, value in source_hashes.items():
            if path in expected_hashes:
                assert value == expected_hashes[path], path
    assert len(validations) + len(pending) == len(species)*len(args.groups)*len(TOOLS)
    validate_metric_grid(metrics, validations)
    write_table(output/'metrics.tsv', metrics)
    write_table(output/'diagnostics.tsv', diagnostics)
    write_table(output/'inputs.tsv', input_rows)
    write_table(output/'reference_regions.tsv', ref_rows)
    write_table(output/'validation.tsv', validations)
    (output/'source_sha256.json').write_text(json.dumps(source_hashes, indent=2)+'\n')
    (output/'STATUS.json').write_text(json.dumps(dict(version=VERSION, cohort=args.cohort,
        status='PARTIAL' if pending else 'COMPLETE', expected_cases=len(species)*len(args.groups)*len(TOOLS),
        completed_cases=len(validations), pending=pending, min_record_bp=MIN_RECORD,
        min_cds_bp=MIN_CDS, confidence_threshold=CONFIDENCE,
        source_files_unchanged=True, fragment_reconstruction='PASS', independent_matching='PASS'), indent=2)+'\n')
    print(f'{len(validations)}/{len(species)*len(args.groups)*len(TOOLS)} cases evaluated.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=common.REPO/'data')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cohort', required=True)
    parser.add_argument('--groups', nargs='+', default=GROUPS, choices=(*GROUPS, '2-5kb'))
    parser.add_argument('--samples', type=Path)
    parser.add_argument('--completed-only', action='store_true', help='Require native COMPLETE markers for live inference directories')
    run(parser.parse_args())
