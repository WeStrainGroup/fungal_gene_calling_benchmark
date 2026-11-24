#!/usr/bin/env python3
"""Create one current input/output view, reusing all completed native predictions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

GROUPS = ('genome', '10-100kb', '5-10kb', '2-5kb', '1-2kb')
TOOLS = ('annevo', 'tiberius', 'carbon', 'helixer')

def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)

def link(source, destination):
    assert source.exists(), source
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        assert destination.resolve() == source.resolve(), destination
    else:
        assert not destination.exists(), destination
        destination.symlink_to(source)

def prepare(base, root, metadata):
    samples, sources = {}, {}
    for manifest, cohort in [('pilot_samples.tsv', 'three_predictor_full_genome'),
                             ('tiberius_samples.tsv', 'fungi_test_six'),
                             ('carbon_samples.tsv', 'carbon_fungi_set'),
                             ('helixer_samples.tsv', root.name)]:
        for row in csv.DictReader((metadata/manifest).open(), delimiter='\t'):
            species = row['species']
            if species not in samples:
                samples[species] = row
                sources[species] = base/cohort
    assert len(samples) == 17
    reused = []
    for species, source in sources.items():
        if source == root:
            continue
        link(source/'provenance'/species, root/'provenance'/species)
        link(source/'inputs'/species/'codon_table.json', root/'inputs'/species/'codon_table.json')
        for group in GROUPS:
            old = base/f'fragment_{group.replace("-", "_")}' if group in ('2-5kb', '10-100kb') else source
            for name in [f'{group}.fna'] + ([] if group == 'genome' else [f'{group}.fragments.tsv']):
                link(old/'inputs'/species/name, root/'inputs'/species/name)
            for tool in TOOLS[:-1]:
                native = old/'out'/tool/species/group
                gff = native/('prediction.gff' if tool == 'carbon' else 'prediction.gff3')
                assert gff.is_file() and gff.stat().st_size, gff
                link(native, root/'out'/tool/species/group)
                reused.append(dict(species=species, group=group, tool=tool,
                                   source=str(gff), sha256=hashlib.sha256(gff.read_bytes()).hexdigest()))
    assert len(reused) == 180
    table(root/'samples.tsv', list(samples.values()))
    tasks = [dict(species=s, group=g, tool=t) for s in samples for g in GROUPS for t in TOOLS
             if t == 'helixer' or sources[s] == root]
    assert len(tasks) == 160
    table(root/'tasks.tsv', tasks)
    (root/'provenance/reused_predictions.json').write_text(json.dumps(reused, indent=2)+'\n')
    if not (root/'INPUTS_READY').exists():
        print('Existing inputs linked; new RefSeq preparation is still running.')
        return
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (root/'inputs').rglob('*') if p.is_file()}
    (root/'input_sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    (root/'BENCHMARK_READY').write_text('17 species; 340 cases; 180 reused; 160 new predictions.\n')
    print((root/'BENCHMARK_READY').read_text())

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', type=Path, required=True)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    a = p.parse_args()
    prepare(a.base.resolve(), a.data.resolve(), a.metadata.resolve())
