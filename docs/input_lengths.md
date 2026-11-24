# Original assembly sequence lengths

Sequence lengths within each original assembly after retaining nuclear sequences ≥1 kb, before random fragmentation. Lengths are in kb, rounded to one decimal place.

| Species | Sequences | Min | Median | Max |
| --- | ---: | ---: | ---: | ---: |
| *Saccharomyces cerevisiae* | 16 | 230.2 | 765.0 | 1,531.9 |
| *Candida albicans* | 8 | 949.5 | 1,701.3 | 3,188.3 |
| *Agaricus bisporus* | 1,245 | 1.0 | 1.7 | 2,973.6 |
| *Cryphonectria parasitica* | 26 | 2.6 | 56.8 | 7,437.3 |
| *Parastagonospora nodorum* | 107 | 2.0 | 32.4 | 2,531.9 |
| *Puccinia striiformis* | 184 | 2.5 | 57.2 | 5,908.9 |
| *Punctularia strigosozonata* | 195 | 2.0 | 3.0 | 4,675.7 |
| *Tilletiopsis washingtonensis* | 38 | 1.2 | 523.7 | 1,093.0 |
| *Schizosaccharomyces pombe* | 3 | 2,452.9 | 4,539.8 | 5,579.1 |
| *Neurospora crassa* | 20 | 1.6 | 78.6 | 9,798.9 |
| *Coprinopsis cinerea* | 67 | 2.2 | 6.7 | 4,147.0 |
| *Rhizophagus irregularis* | 32 | 2,669.0 | 4,707.8 | 7,995.2 |
| *Aspergillus terreus* | 26 | 12.1 | 1,388.6 | 2,751.8 |
| *Coccidioides immitis* | 6 | 3,458.9 | 4,138.2 | 8,482.3 |
| *Debaryomyces hansenii* | 7 | 1,249.9 | 1,606.3 | 2,305.8 |
| *Fusarium graminearum* | 19 | 2.1 | 7.0 | 11,697.3 |
| *Laccaria bicolor* | 665 | 1.0 | 8.2 | 3,566.4 |


The [source table](../results/input_lengths.tsv) retains full-precision values. The displayed rows use `group=genome` and `scope=scored`, separately for each species.

## Reproduce

After extracting the frozen data release, run from the repository root:

```bash
python scripts/summarize_input_lengths.py --data data \
  --samples data/samples.tsv --output results/input_lengths.tsv
```

The script reads the actual FASTA files without modifying them and verifies that fragmentation preserves each genome's total bp. Counts and bp totals were also reconciled against all 85 rows of the published `results/inputs.tsv`.
All 180 five-number summaries were independently checked with NumPy; fragmented-record lengths were also derived from the frozen coordinate maps. Records below each target lower bound were confirmed to be intact original records.
