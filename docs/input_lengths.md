# Original assembly sequence lengths

Sequence lengths within each original assembly after retaining nuclear sequences ≥1 kb, before random fragmentation. Lengths are in kb; values below 10 kb retain one decimal place, and larger values are rounded to integers.

| Species | Sequences | Min | Median | Max |
| --- | ---: | ---: | ---: | ---: |
| *Agaricus bisporus* | 1,245 | 1.0 | 1.7 | 2,974 |
| *Punctularia strigosozonata* | 195 | 2.0 | 3.0 | 4,676 |
| *Coprinopsis cinerea* | 67 | 2.2 | 6.7 | 4,147 |
| *Fusarium graminearum* | 19 | 2.1 | 7.0 | 11,697 |
| *Laccaria bicolor* | 665 | 1.0 | 8.2 | 3,566 |
| *Parastagonospora nodorum* | 107 | 2.0 | 32 | 2,532 |
| *Cryphonectria parasitica* | 26 | 2.6 | 57 | 7,437 |
| *Puccinia striiformis* | 184 | 2.5 | 57 | 5,909 |
| *Neurospora crassa* | 20 | 1.6 | 79 | 9,799 |
| *Tilletiopsis washingtonensis* | 38 | 1.2 | 524 | 1,093 |
| *Saccharomyces cerevisiae* | 16 | 230 | 765 | 1,532 |
| *Aspergillus terreus* | 26 | 12 | 1,389 | 2,752 |
| *Debaryomyces hansenii* | 7 | 1,250 | 1,606 | 2,306 |
| *Candida albicans* | 8 | 950 | 1,701 | 3,188 |
| *Coccidioides immitis* | 6 | 3,459 | 4,138 | 8,482 |
| *Schizosaccharomyces pombe* | 3 | 2,453 | 4,540 | 5,579 |
| *Rhizophagus irregularis* | 32 | 2,669 | 4,708 | 7,995 |


The [source table](../results/input_lengths.tsv) retains full-precision values. The displayed rows use `group=genome` and `scope=scored`, separately for each species.

## Reproduce

After extracting the frozen data release, run from the repository root:

```bash
python scripts/summarize_input_lengths.py --data data \
  --samples data/samples.tsv --output results/input_lengths.tsv
```

The script reads the actual FASTA files without modifying them and verifies that fragmentation preserves each genome's total bp. Counts and bp totals were also reconciled against all 85 rows of the published `results/inputs.tsv`.
All 180 five-number summaries were independently checked with NumPy; fragmented-record lengths were also derived from the frozen coordinate maps. Records below each target lower bound were confirmed to be intact original records.
