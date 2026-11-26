#!/usr/bin/env python3
"""Take within-input replicate medians, then equal-species means."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev
from math import isfinite

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--metrics', nargs='+', required=True, type=Path, help='One combined TSV per round')
p.add_argument('--output', required=True, type=Path)
p.add_argument('--expected-rounds', type=int, choices=[1, 3], default=3)
p.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1]/'metadata/efficiency_inputs.tsv')
a = p.parse_args()
assert len(a.metrics) == a.expected_rounds and len({x.resolve() for x in a.metrics}) == len(a.metrics)
groups = ['genome', '10-100kb', '5-10kb', '2-5kb', '1-2kb']
tools = ['annevo', 'tiberius', 'helixer', 'carbon']
with a.manifest.open() as handle:
    manifest = {(r['species'], r['group']): (r['records'], r['bp'])
                for r in csv.DictReader(handle, delimiter='\t')}
replicates = defaultdict(list)
inputs = defaultdict(set)
job_ids = set()
expected_species = None
for path in a.metrics:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    species = {r['species'] for r in rows}
    expected_species = species if expected_species is None else expected_species
    assert species == expected_species and len(species) == 17
    keys = {(r['species'], r['group'], r['tool']) for r in rows}
    assert len(rows) == len(keys) == 340
    assert keys == {(s, g, t) for s in species for g in groups for t in tools}
    assert all(r['passed'] == 'True' and r['exit_code'] == '0' for r in rows)
    current_jobs = {r['job_id'] for r in rows}
    assert len(current_jobs) == 17 and not job_ids.intersection(current_jobs)
    job_ids.update(current_jobs)
    for r in rows:
        assert isfinite(float(r['elapsed_seconds'])) and float(r['elapsed_seconds']) > 0
        assert isfinite(float(r['max_process_rss_kib']))
        assert (r['records'], r['bp']) == manifest[r['species'], r['group']]
        inputs[r['species'], r['group']].add((r['records'], r['bp']))
        replicates[r['species'], r['group'], r['tool']].append(r)
assert all(len(v) == 1 for v in inputs.values())
medians = []
for (species, group, tool), rows in sorted(replicates.items()):
    assert len(rows) == a.expected_rounds
    seconds = median(float(r['elapsed_seconds']) for r in rows)
    medians.append(dict(species=species, group=group, tool=tool, replicates=len(rows),
        seconds=seconds, minutes=seconds / 60,
        mbp_per_minute=int(rows[0]['bp']) * 60 / seconds / 1e6,
        peak_process_rss_gib=median(float(r['max_process_rss_kib']) / 2**20 for r in rows)))
summary = []
for group in groups:
    for tool in tools:
        rows = [r for r in medians if r['group'] == group and r['tool'] == tool]
        summary.append(dict(group=group, tool=tool, species=len(rows), replicates=a.expected_rounds,
            mean_seconds=mean(r['seconds'] for r in rows),
            sd_seconds=stdev(r['seconds'] for r in rows),
            mean_minutes=mean(r['minutes'] for r in rows),
            sd_minutes=stdev(r['minutes'] for r in rows),
            mean_mbp_per_minute=mean(r['mbp_per_minute'] for r in rows),
            sd_mbp_per_minute=stdev(r['mbp_per_minute'] for r in rows),
            mean_peak_process_rss_gib=mean(r['peak_process_rss_gib'] for r in rows)))
assembly = [dict(r) for r in summary if r['group'] == 'genome']
fastest = min(r['mean_minutes'] for r in assembly)
for r in assembly:
    r['relative_time'] = r['mean_minutes'] / fastest
a.output.mkdir(parents=True, exist_ok=True)
for name, rows in [('species_medians.tsv', medians), ('summary.tsv', summary), ('assembly_time.tsv', assembly)]:
    with (a.output/name).open('w') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

# Round presentation values only; keep machine-readable statistics at full precision.
labels = {'annevo': 'ANNEVO', 'tiberius': 'Tiberius', 'helixer': 'Helixer', 'carbon': 'Carbon-A'}
lines = ['| Tool | Time per original assembly (min, mean ± SD) | Relative elapsed time |',
         '| --- | ---: | ---: |']
for r in assembly:
    lines.append(f"| {labels[r['tool']]} | {r['mean_minutes']:.1f} ± {r['sd_minutes']:.1f} | {r['relative_time']:.1f}× |")
(a.output/'assembly_time.md').write_text('\n'.join(lines)+'\n')
