# Datasets and training exposure

All inputs are fungal nuclear DNA. Cohort labels refer to the authors' published species choices; metrics and fragmentation are computed here with a common protocol.

## ANNEVO-selected fungal species

| Species | Assembly | Code |
| --- | --- | ---: |
| Saccharomyces cerevisiae | GCF_000146045.2_R64 | 1 |
| Candida albicans | GCF_000182965.3_ASM18296v3 | 12 |

Source: [ANNEVO performance table](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/docs/performance.md) and [evaluation notes](https://github.com/xjtu-omics/ANNEVO/blob/42c920f70f52f96c7c06c665a639ec47d9304fff/docs/notes.md). ANNEVO states its test species were excluded from training/validation. Both are explicitly marked as Tiberius training species. The assembly versions above are those actually evaluated here, not an assertion that every upstream benchmark used identical annotation releases.

## Tiberius-selected primary fungal benchmark

| Species | Assembly | Code |
| --- | --- | ---: |
| Agaricus bisporus | GCF_000300555.1 | 1 |
| Cryphonectria parasitica | GCF_011745365.1 | 1 |
| Parastagonospora nodorum | GCF_000146915.1 | 1 |
| Puccinia striiformis | GCF_021901695.1 | 1 |
| Punctularia strigosozonata | GCF_000264995.1 | 1 |
| Tilletiopsis washingtonensis | GCF_003144115.1 | 1 |

Source: [frozen Tiberius fungi configuration](https://github.com/Gaius-Augustus/Tiberius/blob/b3ea48f263a075260dbf1e91a412e82ccdc96d3c/model_cfg/fungi.yaml). This is the six-species primary comparison, not the 228-accession extended table. Tiberius's documentation states that none of these six species were used to train Tiberius or ANNEVO. We attribute that statement to the benchmark authors; it is not an independent complete audit of the ANNEVO 2.3.3 training manifest, nor a claim that no close relative was used.

Input FASTA/GFF downloads use the exact RefSeq accession versions and NCBI MD5 files. Assemblies may contain short original scaffolds, retained unchanged in longer target bins. The extension reports their number and bp contribution.

## Carbon-A official fungal evaluation data

The [official data bucket](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-training-data) reports **5 seen fungi and 0 unseen fungi**. It also lists 644 fungal training assemblies. Thus this set can support within-assembly fragmentation testing, but cannot be presented as unseen-fungal-species evaluation. The five selected assemblies are:

| Species | Assembly | Handling |
| --- | --- | --- |
| Saccharomyces cerevisiae S288C | GCF_000146045.2 | Reuse the completed assembly and frozen references. |
| Schizosaccharomyces pombe | GCF_000002945.2 | Completed all four input conditions. |
| Neurospora crassa OR74A | GCF_000182925.2 | Completed all four input conditions. |
| Coprinopsis cinerea okayama7#130 | GCF_000182895.1 | Completed all four input conditions. |
| Rhizophagus irregularis DAOM-197198 | GCF_026210795.1 | Completed all four input conditions. |

Membership was checked against the official `eval/seen` directory and NCBI assembly metadata. Candida albicans is not in this five-species set. Our common protocol uses frozen RefSeq annotation files for the selected assembly versions; the original Carbon-A bucket may use a different annotation release. This is not an exact reproduction of Carbon-A's published scores.

Annotation quality varies: S. cerevisiae and S. pombe have specialist community resources; the N. crassa assembly metadata explicitly describes largely automatic annotation. Do not treat every reference model as experimentally verified. The set is suitable for paired fragmentation comparisons, with this limitation stated.

A newer assembly release or a different strain does not by itself establish an unseen species. No claim is made that the selected fungi are absent from all three models' training data.
