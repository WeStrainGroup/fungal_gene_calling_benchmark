# Input sequence lengths

Statistics describe the frozen nuclear FASTA records actually used in the benchmark, which may be chromosomes, scaffolds or contigs. They are not CDS lengths or total genome sizes. Shared genomes are counted once.

The [full table](../results/input_lengths.tsv) provides each of the 17 species and the pooled `ALL` distribution for all five input conditions. `scope=all` includes every input record; `scope=scored` includes only records ≥1,000 bp. Quartiles use linear interpolation (R `quantile(type=7)`) with equal weight per record. N50 is the record length at which the cumulative length, ordered longest first, reaches half the total bp; it is not the median record length.

Across original inputs there are 3,435 records totaling 690,727,699 bp, with a minimum of 499 bp. Scoring retains 2,664 records totaling 690,167,379 bp, with a minimum of 1,000 bp. The excluded 771 records (560,320 bp) all belong to *Agaricus bisporus*. The scored original records have a pooled median of 3,236.5 bp but an N50 of 3,600,274 bp: short records are numerous while long records contain most bases.

## Original assemblies by species

All nuclear records, including those below 1 kb. Length quantiles and N50 are in kb; total length is in Mb. Exact bp values and the ≥1 kb subset are in the full table.

| Species | Records | Total Mb | Min | Q1 | Median | Q3 | Max | N50 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| *Saccharomyces cerevisiae* | 16 | 12.071 | 230.218 | 531.954 | 765.042 | 980.594 | 1,531.933 | 924.431 |
| *Candida albicans* | 8 | 14.283 | 949.511 | 1,151.457 | 1,701.279 | 2,245.472 | 3,188.341 | 2,231.883 |
| *Agaricus bisporus* | 2,016 | 32.614 | 0.499 | 0.794 | 1.191 | 1.884 | 2,973.556 | 1,225.131 |
| *Cryphonectria parasitica* | 26 | 43.870 | 2.579 | 7.304 | 56.798 | 3,359.665 | 7,437.273 | 5,118.729 |
| *Parastagonospora nodorum* | 107 | 37.164 | 2.005 | 4.029 | 32.376 | 474.105 | 2,531.949 | 1,045.363 |
| *Puccinia striiformis* | 184 | 89.491 | 2.546 | 33.927 | 57.189 | 84.298 | 5,908.853 | 4,700.190 |
| *Punctularia strigosozonata* | 195 | 34.172 | 2.001 | 2.430 | 3.005 | 4.263 | 4,675.683 | 2,122.796 |
| *Tilletiopsis washingtonensis* | 38 | 18.762 | 1.171 | 249.237 | 523.680 | 738.220 | 1,093.019 | 728.438 |
| *Schizosaccharomyces pombe* | 3 | 12.572 | 2,452.883 | 3,496.344 | 4,539.804 | 5,059.469 | 5,579.133 | 4,539.804 |
| *Neurospora crassa* | 20 | 41.038 | 1.646 | 9.294 | 78.550 | 4,311.148 | 9,798.893 | 6,000.761 |
| *Coprinopsis cinerea* | 67 | 36.150 | 2.208 | 3.452 | 6.705 | 9.392 | 4,146.986 | 3,468.139 |
| *Rhizophagus irregularis* | 32 | 146.773 | 2,669.031 | 3,614.262 | 4,707.757 | 5,270.495 | 7,995.222 | 5,085.394 |
| *Aspergillus terreus* | 26 | 29.331 | 12.095 | 30.166 | 1,388.580 | 1,860.592 | 2,751.824 | 1,912.493 |
| *Coccidioides immitis* | 6 | 28.948 | 3,458.857 | 3,590.154 | 4,138.235 | 5,026.670 | 8,482.323 | 4,323.945 |
| *Debaryomyces hansenii* | 7 | 12.152 | 1,249.940 | 1,465.962 | 1,606.296 | 2,029.283 | 2,305.761 | 2,007.515 |
| *Fusarium graminearum* | 19 | 36.458 | 2.080 | 4.498 | 7.037 | 14.095 | 11,697.295 | 8,914.601 |
| *Laccaria bicolor* | 665 | 64.877 | 1.015 | 3.475 | 8.175 | 17.798 | 3,566.367 | 784.691 |

## Reproduce

After extracting the frozen data release, run from the repository root:

```bash
python scripts/summarize_input_lengths.py --data data \
  --samples data/samples.tsv --output results/input_lengths.tsv
```

The script reads the actual FASTA files without modifying them and verifies that fragmentation preserves each genome's total bp. Counts and bp totals were also reconciled against all 85 rows of the published `results/inputs.tsv`.
All 180 five-number summaries were independently checked with NumPy; fragmented-record lengths were also derived from the frozen coordinate maps. Records below each target lower bound were confirmed to be intact original records.
