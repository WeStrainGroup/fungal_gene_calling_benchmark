# Third-party notices

The repository's MIT license covers our benchmark scripts and documentation. It does not relicense external predictors, model weights or source biological datasets.

## Carbon-A evaluator

`vendor/carbon_a/{evaluation.py,utils.py,constants.py}` are unmodified copies from [HuggingFaceBio/Carbon-A-1.2B](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B/tree/ee6359990109eef6adccfd19aba8ca64b48ae514/annotation_pipeline/src). Upstream licenses this project under Apache 2.0; the full license is included in `vendor/carbon_a/LICENSE`. Revision and individual SHA256 hashes are in `metadata/carbon_evaluator.json`. Carbon-A's original authors retain their rights. Our wrapper changes input normalization/reference selection around the evaluator; it does not modify these files.

## Other software and data

ANNEVO has its own non-commercial license; Tiberius and all dependencies retain their upstream licenses. No predictor weights or full ANNEVO/Tiberius source trees are redistributed here. Follow their official installation and licensing requirements.

The release includes public RefSeq references with original accession/version and checksum information, frozen fragment coordinates, and generated annotation results. Consult [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) and the original submissions for source-data conditions. Generated benchmark tables and figures may be reused with attribution to this repository and release; upstream dataset rights remain unchanged.
