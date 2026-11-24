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

The same frozen Tiberius configuration identifies Agaricus bisporus, Parastagonospora nodorum, Punctularia strigosozonata and Tilletiopsis washingtonensis as used in Helixer training or validation. This overlap is retained and disclosed.

Input FASTA/GFF downloads use the exact RefSeq accession versions and NCBI MD5 files. Assemblies may contain short original scaffolds, retained unchanged in longer target bins. The extension reports their number and bp contribution.

## Carbon-A official fungal evaluation data

The [official data bucket](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-training-data) reports **5 seen fungi and 0 unseen fungi**. It also lists 644 fungal training assemblies. Thus this set can support within-assembly fragmentation testing, but cannot be presented as unseen-fungal-species evaluation. The five selected assemblies are:

| Species | Assembly | Handling |
| --- | --- | --- |
| Saccharomyces cerevisiae S288C | GCF_000146045.2 | Reuse the completed assembly and frozen references. |
| Schizosaccharomyces pombe | GCF_000002945.2 | Reuse the five completed input conditions. |
| Neurospora crassa OR74A | GCF_000182925.2 | Reuse the five completed input conditions. |
| Coprinopsis cinerea okayama7#130 | GCF_000182895.1 | Reuse the five completed input conditions. |
| Rhizophagus irregularis DAOM-197198 | GCF_026210795.1 | Reuse the five completed input conditions. |

Membership was checked against the official `eval/seen` directory and NCBI assembly metadata. Candida albicans is not in this five-species set. Our common protocol uses frozen RefSeq annotation files for the selected assembly versions; the original Carbon-A bucket may use a different annotation release. This is not an exact reproduction of Carbon-A's published scores.

Annotation quality varies: S. cerevisiae and S. pombe have specialist community resources; the N. crassa assembly metadata explicitly describes largely automatic annotation. Do not treat every reference model as experimentally verified. The set is suitable for paired fragmentation comparisons, with this limitation stated.

A newer assembly release or a different strain does not by itself establish an unseen species. No claim is made that the selected fungi are absent from all four models' training data.

## Helixer-selected fungal test species

Use the six fungal rows explicitly marked **test (used)** in [Nature Methods Supplementary Table S17](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-025-02939-1/MediaObjects/41592_2025_2939_MOESM1_ESM.pdf), rather than all reserved-test rows.

| Species | Assembly | Handling |
| --- | --- | --- |
| Aspergillus terreus | GCF_000149615.1 | New; genetic code 1. |
| Coccidioides immitis | GCF_000149335.2 | New; genetic code 1. |
| Debaryomyces hansenii | GCF_000006445.2 | New; genetic code 12, checked against nuclear CDS qualifiers. |
| Fusarium graminearum | GCF_000240135.3 | New; genetic code 1. |
| Laccaria bicolor | GCF_000143565.1 | New; genetic code 1. |
| Neurospora crassa | GCF_000182925.2 | Already present; reuse the exact frozen input and three-tool predictions. |

This labels Helixer's published test selection, not absence from ANNEVO, Tiberius or Carbon-A training. The frozen Tiberius fungi training list explicitly includes Coccidioides immitis, Debaryomyces hansenii and Neurospora crassa. Cross-tool training exposure for all additions is not yet fully audited; unknown exposure must not be called unseen. Latest Helixer software uses the published fungi checkpoint; changing postprocessing does not make the checkpoint newly trained.
