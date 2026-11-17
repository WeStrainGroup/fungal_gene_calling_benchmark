"""Evaluate predicted CDS annotations against an NCBI reference assembly.

When ``sequence_report.jsonl`` is available, evaluation prefers assembled
chromosomes, then assembled linkage groups, then chromosome scaffolds; without
the report every record is evaluated. Six complementary metrics are reported:

* gene: one-to-one exact CDS-sequence matches on the same genomic record;
* protein: gene exact matches plus additional one-to-one exact matches after
  translating unmatched predicted and reference CDS with the same codon table;
* exon: exact merged CDS intervals from the longest CDS isoform per gene;
* 5' boundary and 3' boundary: individual CDS boundaries in transcript direction;
* nucleotide: base-pair overlap against those reference intervals.

Gene and protein metrics count the predicted genes whose ``confidence``
reaches the confidence threshold (default 0.1); genes without a confidence are
always counted, and a threshold of 0 counts every predicted gene. Exon,
boundary, and nucleotide metrics describe the whole predicted annotation.
Evaluation also reports AUROC, non-interpolated AUPRC, and point-biserial
correlation of the confidence against the one-to-one gene exact-match label,
together with exact and non-exact group means.

The 5' and 3' boundary metrics follow transcript direction, so their meanings are
reversed on the negative genomic strand. Interval arithmetic keeps memory
proportional to the number of CDS intervals rather than genome length.

Example:

    python src/evaluation.py --pred results/GCF_000146045.2/prediction \
        --genome_dir ncbi_dataset/data/GCF_000146045.2
"""

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import sys


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from constants import CONFIDENCE_THRESHOLD
from utils import load_codon_table, load_nuclear_evaluation_accessions


STRANDS = ("+", "-")
REFERENCE_CDS_FEATURES = frozenset(("mRNA", "transcript", "CDS"))
METRIC_NAMES = (
    "gene",
    "protein",
    "exon",
    "five_prime_boundary",
    "three_prime_boundary",
    "nucleotide",
)
INTERVAL_METRIC_NAMES = (
    "exon",
    "five_prime_boundary",
    "three_prime_boundary",
    "nucleotide",
)
REFERENCE_RECORD_PATTERN = re.compile(r"^lcl\|(.+?)_cds_")
PREDICTED_RECORD_PATTERN = re.compile(r"(?:^|\s)record_id=(\S+)")
PREDICTED_FRAME_PATTERN = re.compile(r"(?:^|\s)phase=([012])(?:\s|$)")
REFERENCE_FRAME_PATTERN = re.compile(r"\[frame=([123])\]")
GENE_ID_PATTERN = re.compile(
    r"\[db_xref=(?:[^,\]]+,)*GeneID:(\d+)(?:,|\])"
)
LOCUS_TAG_PATTERN = re.compile(r"\[locus_tag=([^\]]+)\]")
GENE_NAME_PATTERN = re.compile(r"\[gene=([^\]]+)\]")
CONFIDENCE_NAME = "confidence"


def iter_gff_rows(path):
    """Yield the nine tab-separated fields from each GFF3 feature row."""
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#") or line.isspace():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 GFF3 fields, "
                    f"found {len(fields)}"
                )
            yield fields


def gff_attribute_value(attribute_text, name):
    """Return one GFF3 attribute without parsing unrelated attributes."""
    prefix = name + "="
    if attribute_text.startswith(prefix):
        value_start = len(prefix)
    else:
        marker = ";" + prefix
        marker_start = attribute_text.find(marker)
        if marker_start == -1:
            return None
        value_start = marker_start + len(marker)

    value_end = attribute_text.find(";", value_start)
    if value_end == -1:
        return attribute_text[value_start:]
    return attribute_text[value_start:value_end]


def iter_selected_gff_features(path, feature_types, retained_records=None):
    """Yield selected GFF3 features, optionally restricted by record id.

    NCBI GFF files contain many gene, exon, and regulatory rows that CDS
    evaluation never uses. Reading the first three columns lets those rows be
    rejected before splitting their long attribute columns.
    """
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#") or line.isspace():
                continue

            leading_fields = line.split("\t", 3)
            if len(leading_fields) != 4:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 GFF3 fields, "
                    f"found {len(leading_fields)}"
                )
            record_id, _, feature_type, remaining_text = leading_fields
            if feature_type not in feature_types:
                continue
            if retained_records is not None and record_id not in retained_records:
                continue

            remaining_fields = remaining_text.rstrip("\n").split("\t", 5)
            if len(remaining_fields) != 6:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 GFF3 fields, "
                    f"found {3 + len(remaining_fields)}"
                )
            start, end, _, strand, _, attribute_text = remaining_fields
            yield record_id, feature_type, start, end, strand, attribute_text


def parse_predicted_cds_tracks(path):
    """Return predicted CDS intervals keyed by ``(record_id, strand)``."""
    tracks = defaultdict(list)
    for fields in iter_gff_rows(path):
        record_id, _, feature_type, start, end, _, strand, _, _ = fields
        if feature_type != "CDS":
            continue
        if strand not in STRANDS:
            raise ValueError(f"unsupported CDS strand {strand!r} in {path}")
        tracks[(record_id, strand)].append((int(start) - 1, int(end) - 1))
    return dict(tracks)


def parse_reference_cds_tracks(path, retained_records=None):
    """Return CDS intervals from the longest CDS isoform of each reference gene."""
    transcript_gene = {}
    transcript_cds = defaultdict(list)

    for feature in iter_selected_gff_features(
        path,
        REFERENCE_CDS_FEATURES,
        retained_records,
    ):
        record_id, feature_type, start, end, strand, attribute_text = feature
        if feature_type in ("mRNA", "transcript"):
            transcript_id = gff_attribute_value(attribute_text, "ID")
            parent = gff_attribute_value(attribute_text, "Parent")
            if transcript_id and parent:
                transcript_gene[(record_id, transcript_id)] = parent.split(",", 1)[0]
            continue

        if strand not in STRANDS:
            raise ValueError(f"unsupported CDS strand {strand!r} in {path}")
        parent = gff_attribute_value(attribute_text, "Parent")
        for transcript_id in (parent or "").split(","):
            if transcript_id:
                transcript_cds[(record_id, transcript_id)].append(
                    (int(start) - 1, int(end) - 1, strand)
                )

    transcripts_by_gene = defaultdict(dict)
    for (record_id, transcript_id), intervals in transcript_cds.items():
        gene_id = transcript_gene.get((record_id, transcript_id), transcript_id)
        transcripts_by_gene[(record_id, gene_id)][transcript_id] = intervals

    tracks = defaultdict(list)
    for (record_id, _), transcripts in transcripts_by_gene.items():
        longest_transcript = max(
            transcripts,
            key=lambda transcript_id: sum(
                end - start + 1
                for start, end, _ in transcripts[transcript_id]
            ),
        )
        intervals = transcripts[longest_transcript]
        # NCBI trans-spliced transcripts can contain CDS pieces on both strands.
        # Each piece remains on its annotated genomic strand for interval scoring.
        for start, end, strand in intervals:
            tracks[(record_id, strand)].append((start, end))
    return dict(tracks)


def iter_fasta_records(path):
    """Yield ``(header_without_marker, sequence)`` from a FASTA file."""
    header = None
    sequence_lines = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(line.strip())
    if header is not None:
        yield header, "".join(sequence_lines)


def parse_predicted_cds_fasta(path):
    """Return every predicted CDS as ``(record_id, sequence)``.

    A list is intentional: two identical predictions remain two predictions and
    therefore participate independently in precision accounting.
    """
    return [
        (record_id, sequence)
        for _, record_id, sequence in parse_predicted_cds_records(path)
    ]


def parse_predicted_cds_records(path):
    """Return predicted CDS records as ``(gene_id, record_id, sequence)``."""
    return [
        (gene_id, record_id, sequence)
        for gene_id, record_id, sequence, _ in (
            parse_predicted_cds_records_with_frames(path)
        )
    ]


def parse_predicted_cds_records_with_frames(path):
    """Return predicted CDS records with their zero-based translation frame."""
    predictions = []
    for header, sequence in iter_fasta_records(path):
        record_match = PREDICTED_RECORD_PATTERN.search(header)
        if record_match is None:
            raise ValueError(f"predicted FASTA header lacks record_id: {header}")
        frame_match = PREDICTED_FRAME_PATTERN.search(header)
        frame = int(frame_match.group(1)) if frame_match else 0
        gene_id = header.split(None, 1)[0]
        predictions.append((
            gene_id,
            record_match.group(1),
            sequence.upper(),
            frame,
        ))
    return predictions


def translate_cds_for_evaluation(sequence, frame, codon_to_amino_acid):
    """Translate one CDS for exact-match evaluation.

    Incomplete trailing codons are ignored, ambiguous codons become ``X``, and
    one terminal stop is removed. Internal stops remain in the sequence so an
    unusual reference CDS can still be compared instead of aborting evaluation.
    """
    sequence = sequence.upper()
    codon_end = len(sequence) - ((len(sequence) - frame) % 3)
    protein = "".join(
        codon_to_amino_acid.get(sequence[start:start + 3], "X")
        for start in range(frame, codon_end, 3)
    )
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


def parse_predicted_translated_proteins(path, codon_to_amino_acid):
    """Translate predicted CDS records as ``(gene_id, record_id, protein)``."""
    return [
        (
            gene_id,
            record_id,
            translate_cds_for_evaluation(sequence, frame, codon_to_amino_acid),
        )
        for gene_id, record_id, sequence, frame in (
            parse_predicted_cds_records_with_frames(path)
        )
    ]


def parse_predicted_gene_scores(path, score_name=CONFIDENCE_NAME):
    """Return one optional gene-level score keyed by predicted GFF3 gene ID."""
    scores = {}
    for fields in iter_gff_rows(path):
        _, _, feature_type, _, _, _, _, _, attribute_text = fields
        if feature_type != "gene":
            continue
        gene_id = gff_attribute_value(attribute_text, "ID")
        if gene_id is None:
            raise ValueError(
                f"predicted GFF3 gene lacks an ID in {path}"
            )
        if gene_id in scores:
            raise ValueError(f"duplicate predicted gene ID {gene_id!r} in {path}")
        score_text = gff_attribute_value(attribute_text, score_name)
        if score_text is None:
            scores[gene_id] = None
            continue
        if score_text in ("NA", "."):
            scores[gene_id] = None
            continue
        try:
            score = float(score_text)
        except ValueError as exc:
            raise ValueError(
                f"invalid {score_name}={score_text!r} for {gene_id!r} in {path}"
            ) from exc
        if not math.isfinite(score) or not 0.0 <= score <= 1.0:
            raise ValueError(
                f"{score_name} for {gene_id!r} must be between 0 and 1, "
                f"found {score_text!r}"
            )
        scores[gene_id] = score
    return scores


def reference_gene_id(header):
    """Return the most stable gene identifier available in an NCBI CDS header."""
    match = GENE_ID_PATTERN.search(header)
    if match:
        return "GeneID:" + match.group(1)
    match = LOCUS_TAG_PATTERN.search(header)
    if match:
        return "locus_tag:" + match.group(1)
    match = GENE_NAME_PATTERN.search(header)
    if match:
        return "gene:" + match.group(1)
    return "cds:" + header.split()[0]


def parse_reference_cds_fasta(path):
    """Return reference isoform sequences keyed by ``(record_id, gene_id)``."""
    return {
        gene_key: [sequence for sequence, _ in isoforms]
        for gene_key, isoforms in parse_reference_cds_fasta_with_frames(path).items()
    }


def parse_reference_cds_fasta_with_frames(path):
    """Return reference CDS and zero-based frame by genomic record and gene."""
    genes = defaultdict(list)
    for header, sequence in iter_fasta_records(path):
        record_match = REFERENCE_RECORD_PATTERN.search(header)
        if record_match is None:
            raise ValueError(f"unrecognized NCBI CDS FASTA header: {header}")
        frame_match = REFERENCE_FRAME_PATTERN.search(header)
        frame = int(frame_match.group(1)) - 1 if frame_match else 0
        key = (record_match.group(1), reference_gene_id(header))
        genes[key].append((sequence.upper(), frame))
    return dict(genes)


def parse_reference_translated_proteins(path, codon_to_amino_acid):
    """Translate reference CDS isoforms by ``(record_id, gene_id)``."""
    return {
        gene_key: [
            translate_cds_for_evaluation(sequence, frame, codon_to_amino_acid)
            for sequence, frame in isoforms
        ]
        for gene_key, isoforms in parse_reference_cds_fasta_with_frames(path).items()
    }


def merge_intervals(intervals):
    """Return the union of inclusive intervals in coordinate order."""
    if not intervals:
        return []
    ordered = sorted(intervals)
    if any(start > end for start, end in ordered):
        raise ValueError(f"invalid interval list: {ordered}")
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        previous = merged[-1]
        if start <= previous[1] + 1:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def merge_tracks(tracks):
    """Merge overlapping or adjacent CDS intervals independently per track."""
    return {track: merge_intervals(intervals) for track, intervals in tracks.items()}


def interval_overlap_length(left, right):
    """Return the number of bases covered by both merged interval lists."""
    overlap = 0
    left_index = 0
    right_index = 0
    while left_index < len(left) and right_index < len(right):
        left_start, left_end = left[left_index]
        right_start, right_end = right[right_index]
        overlap += max(
            0,
            min(left_end, right_end) - max(left_start, right_start) + 1,
        )
        if left_end < right_end:
            left_index += 1
        else:
            right_index += 1
    return overlap


def add_set_metric_counts(counts, metric_name, predicted_items, reference_items):
    """Accumulate TP/FP/FN counts for two sets."""
    true_positives = len(predicted_items & reference_items)
    counts[f"{metric_name}_tp"] += true_positives
    counts[f"{metric_name}_fp"] += len(predicted_items) - true_positives
    counts[f"{metric_name}_fn"] += len(reference_items) - true_positives


def interval_metric_counts(predicted_tracks, reference_tracks):
    """Count exon, boundary, and nucleotide matches."""
    predicted_tracks = merge_tracks(predicted_tracks)
    reference_tracks = merge_tracks(reference_tracks)
    counts = defaultdict(int)
    for metric_name in INTERVAL_METRIC_NAMES:
        for suffix in ("tp", "fp", "fn"):
            counts[f"{metric_name}_{suffix}"] = 0

    for track in set(predicted_tracks) | set(reference_tracks):
        _, strand = track
        predicted_exons = predicted_tracks.get(track, ())
        reference_exons = reference_tracks.get(track, ())

        predicted_bp = sum(end - start + 1 for start, end in predicted_exons)
        reference_bp = sum(end - start + 1 for start, end in reference_exons)
        overlap_bp = interval_overlap_length(predicted_exons, reference_exons)
        counts["nucleotide_tp"] += overlap_bp
        counts["nucleotide_fp"] += predicted_bp - overlap_bp
        counts["nucleotide_fn"] += reference_bp - overlap_bp

        predicted_left_boundaries = {start for start, _ in predicted_exons}
        predicted_right_boundaries = {end for _, end in predicted_exons}
        reference_left_boundaries = {start for start, _ in reference_exons}
        reference_right_boundaries = {end for _, end in reference_exons}

        if strand == "+":
            predicted_five_prime_boundaries = predicted_left_boundaries
            predicted_three_prime_boundaries = predicted_right_boundaries
            reference_five_prime_boundaries = reference_left_boundaries
            reference_three_prime_boundaries = reference_right_boundaries
        else:
            predicted_five_prime_boundaries = predicted_right_boundaries
            predicted_three_prime_boundaries = predicted_left_boundaries
            reference_five_prime_boundaries = reference_right_boundaries
            reference_three_prime_boundaries = reference_left_boundaries

        add_set_metric_counts(
            counts,
            "five_prime_boundary",
            predicted_five_prime_boundaries,
            reference_five_prime_boundaries,
        )
        add_set_metric_counts(
            counts,
            "three_prime_boundary",
            predicted_three_prime_boundaries,
            reference_three_prime_boundaries,
        )
        add_set_metric_counts(
            counts,
            "exon",
            set(predicted_exons),
            set(reference_exons),
        )

    return dict(counts)


def gene_exact_match_assignments(predicted_cds, reference_genes):
    """Return one-to-one exact matches keyed by prediction index.

    A reference gene may match any one of its annotated isoforms. The augmenting
    path prevents one prediction from recalling several genes with an identical
    CDS and prevents duplicate predictions from matching the same gene twice.
    """
    candidate_genes = defaultdict(list)
    for gene_key, isoforms in reference_genes.items():
        record_id, _ = gene_key
        for sequence in set(isoforms):
            candidate_genes[(record_id, sequence)].append(gene_key)

    matched_prediction_by_gene = {}
    matched_gene_by_prediction = {}

    def match_prediction(start_prediction):
        """Find and apply one augmenting path without recursive calls."""
        queued_predictions = [start_prediction]
        visited_predictions = {start_prediction}
        visited_genes = set()
        parent_prediction_by_gene = {}
        queue_index = 0

        while queue_index < len(queued_predictions):
            prediction_index = queued_predictions[queue_index]
            queue_index += 1
            prediction = predicted_cds[prediction_index]

            for gene_key in candidate_genes.get(prediction, ()):
                if gene_key in visited_genes:
                    continue
                visited_genes.add(gene_key)
                parent_prediction_by_gene[gene_key] = prediction_index

                previous_prediction = matched_prediction_by_gene.get(gene_key)
                if previous_prediction is None:
                    current_gene = gene_key
                    while current_gene is not None:
                        current_prediction = parent_prediction_by_gene[current_gene]
                        previous_gene = matched_gene_by_prediction.get(
                            current_prediction
                        )
                        matched_prediction_by_gene[current_gene] = current_prediction
                        matched_gene_by_prediction[current_prediction] = current_gene
                        current_gene = previous_gene
                    return True

                if previous_prediction not in visited_predictions:
                    visited_predictions.add(previous_prediction)
                    queued_predictions.append(previous_prediction)
        return False

    for prediction_index in range(len(predicted_cds)):
        match_prediction(prediction_index)
    return matched_gene_by_prediction


def gene_metric_counts(
    predicted_cds,
    reference_genes,
    matched_gene_by_prediction=None,
):
    """Count one-to-one exact CDS matches between predictions and genes."""
    if matched_gene_by_prediction is None:
        matched_gene_by_prediction = gene_exact_match_assignments(
            predicted_cds,
            reference_genes,
        )
    matched_genes = len(matched_gene_by_prediction)
    return {
        "gene_tp": matched_genes,
        "gene_fp": len(predicted_cds) - matched_genes,
        "gene_fn": len(reference_genes) - matched_genes,
        "predicted_genes": len(predicted_cds),
        "reference_genes": len(reference_genes),
    }


def protein_metric_counts(
    predicted_proteins,
    reference_proteins,
    inherited_gene_matches=None,
):
    """Count protein matches, retaining every exact CDS gene match.

    Translation frames can be uncertain for partial CDS. Exact nucleotide CDS
    matches supplied in ``inherited_gene_matches`` therefore remain matches at
    protein level even if their independently annotated frames differ. Protein
    matching is then applied only to the still-unmatched predictions and genes,
    allowing synonymous CDS differences to add matches without double counting.
    """
    if inherited_gene_matches is None:
        matched_gene_by_prediction = gene_exact_match_assignments(
            predicted_proteins,
            reference_proteins,
        )
    else:
        matched_gene_by_prediction = dict(inherited_gene_matches)
        matched_reference_genes = set(matched_gene_by_prediction.values())
        if len(matched_reference_genes) != len(matched_gene_by_prediction):
            raise ValueError("inherited gene matches must be one-to-one")
        if any(
            prediction_index < 0 or prediction_index >= len(predicted_proteins)
            for prediction_index in matched_gene_by_prediction
        ):
            raise ValueError("inherited gene match has an invalid prediction index")
        if not matched_reference_genes.issubset(reference_proteins):
            raise ValueError("inherited gene match is absent from protein references")

        unmatched_prediction_indices = [
            prediction_index
            for prediction_index in range(len(predicted_proteins))
            if prediction_index not in matched_gene_by_prediction
        ]
        unmatched_reference_proteins = {
            gene_key: isoforms
            for gene_key, isoforms in reference_proteins.items()
            if gene_key not in matched_reference_genes
        }
        additional_matches = gene_exact_match_assignments(
            [
                predicted_proteins[prediction_index]
                for prediction_index in unmatched_prediction_indices
            ],
            unmatched_reference_proteins,
        )
        matched_gene_by_prediction.update({
            unmatched_prediction_indices[local_prediction_index]: gene_key
            for local_prediction_index, gene_key in additional_matches.items()
        })

    matched_genes = len(matched_gene_by_prediction)
    return {
        "protein_tp": matched_genes,
        "protein_fp": len(predicted_proteins) - matched_genes,
        "protein_fn": len(reference_proteins) - matched_genes,
        "predicted_proteins": len(predicted_proteins),
        "reference_proteins": len(reference_proteins),
    }


def pearson_correlation(left_values, right_values):
    """Return Pearson correlation, or ``None`` when either side is constant."""
    if len(left_values) != len(right_values):
        raise ValueError("correlation vectors must have the same length")
    if len(left_values) < 2:
        return None

    left_mean = math.fsum(left_values) / len(left_values)
    right_mean = math.fsum(right_values) / len(right_values)
    left_deviations = [value - left_mean for value in left_values]
    right_deviations = [value - right_mean for value in right_values]
    left_sum_squares = math.fsum(value * value for value in left_deviations)
    right_sum_squares = math.fsum(value * value for value in right_deviations)
    if left_sum_squares == 0.0 or right_sum_squares == 0.0:
        return None

    cross_product = math.fsum(
        left * right
        for left, right in zip(left_deviations, right_deviations)
    )
    correlation = cross_product / math.sqrt(
        left_sum_squares * right_sum_squares
    )
    return max(-1.0, min(1.0, correlation))


def binary_ranking_auroc_auprc(scores, labels):
    """Return tie-aware AUROC and stepwise AUPRC for binary labels.

    Score ties are processed as a group. This gives tied positive/negative pairs
    half credit in AUROC and prevents input order from affecting either metric.
    AUPRC uses the average-precision (non-interpolated step) convention. Both
    metrics are undefined when either class is absent.
    """
    if len(scores) != len(labels):
        raise ValueError("ranking score and label vectors must have the same length")
    if any(label not in (0, 1, False, True) for label in labels):
        raise ValueError("ranking labels must be binary")

    positive_count = sum(bool(label) for label in labels)
    negative_count = len(labels) - positive_count
    if positive_count == 0 or negative_count == 0:
        return None, None

    ranked = sorted(
        zip(scores, (bool(label) for label in labels)),
        key=lambda item: item[0],
        reverse=True,
    )
    true_positives = 0
    false_positives = 0
    previous_recall = 0.0
    auroc_area = 0.0
    auprc_area = 0.0
    group_start = 0
    while group_start < len(ranked):
        score = ranked[group_start][0]
        group_end = group_start
        group_positives = 0
        while group_end < len(ranked) and ranked[group_end][0] == score:
            group_positives += int(ranked[group_end][1])
            group_end += 1
        group_negatives = group_end - group_start - group_positives

        next_true_positives = true_positives + group_positives
        next_false_positives = false_positives + group_negatives
        false_positive_rate_increment = group_negatives / negative_count
        auroc_area += false_positive_rate_increment * (
            true_positives + next_true_positives
        ) / (2.0 * positive_count)

        recall = next_true_positives / positive_count
        precision = next_true_positives / (
            next_true_positives + next_false_positives
        )
        auprc_area += (recall - previous_recall) * precision

        true_positives = next_true_positives
        false_positives = next_false_positives
        previous_recall = recall
        group_start = group_end

    return auroc_area, auprc_area


def confidence_exact_match_metrics(
    confidence_values,
    matched_gene_by_prediction,
):
    """Summarize association between confidence and exact gene matches."""
    scored_values = []
    exact_labels = []
    exact_scores = []
    nonexact_scores = []
    for prediction_index, score in enumerate(confidence_values):
        if score is None:
            continue
        is_exact = prediction_index in matched_gene_by_prediction
        scored_values.append(score)
        exact_labels.append(float(is_exact))
        if is_exact:
            exact_scores.append(score)
        else:
            nonexact_scores.append(score)

    auroc, auprc = binary_ranking_auroc_auprc(scored_values, exact_labels)
    return {
        "confidence_scored_genes": len(scored_values),
        "confidence_missing_genes": (
            len(confidence_values) - len(scored_values)
        ),
        "confidence_mean": (
            math.fsum(scored_values) / len(scored_values)
            if scored_values
            else None
        ),
        "confidence_exact_genes": len(exact_scores),
        "confidence_nonexact_genes": len(nonexact_scores),
        "confidence_exact_mean": (
            math.fsum(exact_scores) / len(exact_scores)
            if exact_scores
            else None
        ),
        "confidence_nonexact_mean": (
            math.fsum(nonexact_scores) / len(nonexact_scores)
            if nonexact_scores
            else None
        ),
        "confidence_exact_match_correlation": pearson_correlation(
            scored_values,
            exact_labels,
        ),
        "confidence_exact_match_auroc": auroc,
        "confidence_exact_match_auprc": auprc,
    }


def precision_recall_f1(true_positives, false_positives, false_negatives):
    """Return precision, recall, and F1 for integer match counts."""
    precision_denominator = true_positives + false_positives
    recall_denominator = true_positives + false_negatives
    precision = (
        true_positives / precision_denominator
        if precision_denominator
        else 0.0
    )
    recall = (
        true_positives / recall_denominator
        if recall_denominator
        else 0.0
    )
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return precision, recall, f1


def add_metric_rates(metrics):
    """Add precision, recall, and F1 values for every metric count group."""
    for metric_name in METRIC_NAMES:
        precision, recall, f1 = precision_recall_f1(
            metrics.get(f"{metric_name}_tp", 0),
            metrics.get(f"{metric_name}_fp", 0),
            metrics.get(f"{metric_name}_fn", 0),
        )
        metrics[f"{metric_name}_precision"] = precision
        metrics[f"{metric_name}_recall"] = recall
        metrics[f"{metric_name}_f1"] = f1


def filter_tracks(tracks, retained_records):
    """Restrict interval tracks to a record set, or return all tracks for None."""
    if retained_records is None:
        return tracks
    return {
        track: intervals
        for track, intervals in tracks.items()
        if track[0] in retained_records
    }


def filter_gene_inputs(predicted_cds, reference_genes, retained_records):
    """Restrict gene-level inputs to a record set."""
    if retained_records is None:
        return predicted_cds, reference_genes
    return (
        [
            prediction
            for prediction in predicted_cds
            if prediction[0] in retained_records
        ],
        {
            gene_key: isoforms
            for gene_key, isoforms in reference_genes.items()
            if gene_key[0] in retained_records
        },
    )


def required_file(directory, filename):
    """Return a required file path with a useful error when it is absent."""
    path = Path(directory) / filename
    if not path.is_file():
        raise FileNotFoundError(f"required evaluation file not found: {path}")
    return path


def default_codon_table_path(genome_dir):
    """Return an accession-specific codon table when one is bundled."""
    codon_table_directory = Path(__file__).resolve().parent.parent / "codon_tables"
    accession_path = codon_table_directory / f"{Path(genome_dir).name}.json"
    if accession_path.is_file():
        return accession_path
    return codon_table_directory / "standard.json"


def evaluate_prediction(
    prediction_prefix,
    genome_dir,
    confidence_threshold=CONFIDENCE_THRESHOLD,
    codon_table_path=None,
):
    """Load one prediction/reference pair and return count and rate metrics.

    Gene and protein metrics count the predicted genes whose confidence is at
    least ``confidence_threshold``; genes without a confidence
    (``confidence=NA``) are always counted, and ``0`` counts every predicted
    gene. Exon, boundary, and nucleotide metrics use every predicted CDS
    interval. Without ``codon_table_path`` the translation table is
    ``codon_tables/<accession>.json`` when the genome directory is named after
    an accession with a bundled table, otherwise the standard table.
    """
    prediction_prefix = Path(prediction_prefix)
    genome_dir = Path(genome_dir)
    if (
        not math.isfinite(confidence_threshold)
        or not 0.0 <= confidence_threshold <= 1.0
    ):
        raise ValueError("confidence_threshold must be between 0 and 1")
    codon_table_path = (
        Path(codon_table_path)
        if codon_table_path is not None
        else default_codon_table_path(genome_dir)
    )
    codon_to_amino_acid = load_codon_table(codon_table_path)
    retained_records = load_nuclear_evaluation_accessions(genome_dir)

    predicted_gff_path = str(prediction_prefix) + ".gff"
    predicted_gene_scores = parse_predicted_gene_scores(predicted_gff_path)
    predicted_tracks = parse_predicted_cds_tracks(predicted_gff_path)
    predicted_cds_frame_records = [
        prediction
        for prediction in parse_predicted_cds_records_with_frames(
            str(prediction_prefix) + ".fna"
        )
        if confidence_threshold <= 0.0
        or predicted_gene_scores.get(prediction[0]) is None
        or predicted_gene_scores[prediction[0]] >= confidence_threshold
    ]
    reference_tracks = parse_reference_cds_tracks(
        required_file(genome_dir, "genomic.gff"),
        retained_records,
    )
    reference_cds_path = required_file(genome_dir, "cds_from_genomic.fna")
    reference_genes = parse_reference_cds_fasta(reference_cds_path)
    reference_proteins = parse_reference_translated_proteins(
        reference_cds_path,
        codon_to_amino_acid,
    )

    predicted_tracks = filter_tracks(predicted_tracks, retained_records)
    reference_tracks = filter_tracks(reference_tracks, retained_records)
    if retained_records is not None:
        predicted_cds_frame_records = [
            prediction
            for prediction in predicted_cds_frame_records
            if prediction[1] in retained_records
        ]
    predicted_cds = [
        (record_id, sequence)
        for _, record_id, sequence, _ in predicted_cds_frame_records
    ]
    predicted_proteins = [
        (record_id, translate_cds_for_evaluation(sequence, frame, codon_to_amino_acid))
        for _, record_id, sequence, frame in predicted_cds_frame_records
    ]
    _, reference_genes = filter_gene_inputs([], reference_genes, retained_records)
    _, reference_proteins = filter_gene_inputs([], reference_proteins, retained_records)

    matched_gene_by_prediction = gene_exact_match_assignments(
        predicted_cds,
        reference_genes,
    )
    metrics = gene_metric_counts(
        predicted_cds,
        reference_genes,
        matched_gene_by_prediction,
    )
    metrics.update(protein_metric_counts(
        predicted_proteins,
        reference_proteins,
        matched_gene_by_prediction,
    ))
    metrics.update(confidence_exact_match_metrics(
        [
            predicted_gene_scores.get(gene_id)
            for gene_id, _, _, _ in predicted_cds_frame_records
        ],
        matched_gene_by_prediction,
    ))
    metrics.update(interval_metric_counts(predicted_tracks, reference_tracks))
    add_metric_rates(metrics)

    metrics["genome"] = genome_dir.name
    metrics["codon_table"] = codon_table_path.name
    metrics["confidence_threshold"] = confidence_threshold
    metrics["records"] = len(
        {record_id for record_id, _ in predicted_tracks}
        | {record_id for record_id, _ in reference_tracks}
    )
    metrics["record_filter"] = (
        "all (no sequence report)" if retained_records is None else "nuclear"
    )
    return metrics


def write_evaluation_metrics_json(path, metrics):
    """Atomically write the evaluation metrics of one accession."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "accession": metrics["genome"],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "metrics": metrics,
    }
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    with temporary_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    temporary_path.replace(output_path)


def print_evaluation_metrics(metrics):
    """Print a concise summary of one evaluation."""
    print(
        f"[evaluation] {metrics['genome']} filter={metrics['record_filter']} "
        f"codon_table={metrics['codon_table']} "
        f"confidence_threshold={metrics['confidence_threshold']:.2f} "
        f"({metrics['records']} records, "
        f"{metrics['reference_genes']} reference genes, "
        f"{metrics['predicted_genes']} predictions)"
    )
    for label, metric_name in (
        ("gene", "gene"),
        ("protein", "protein"),
        ("exon", "exon"),
        ("5' boundary", "five_prime_boundary"),
        ("3' boundary", "three_prime_boundary"),
        ("nucleotide", "nucleotide"),
    ):
        print(
            f"  {label:<11} P/R/F1: "
            f"{metrics[f'{metric_name}_precision']:.4f} / "
            f"{metrics[f'{metric_name}_recall']:.4f} / "
            f"{metrics[f'{metric_name}_f1']:.4f}"
        )
    scored_genes = metrics.get("confidence_scored_genes")
    if scored_genes:
        correlation = metrics["confidence_exact_match_correlation"]
        correlation_text = "NA" if correlation is None else f"{correlation:.4f}"
        auroc = metrics["confidence_exact_match_auroc"]
        auprc = metrics["confidence_exact_match_auprc"]
        auroc_text = "NA" if auroc is None else f"{auroc:.4f}"
        auprc_text = "NA" if auprc is None else f"{auprc:.4f}"
        exact_mean = metrics["confidence_exact_mean"]
        nonexact_mean = metrics["confidence_nonexact_mean"]
        exact_mean_text = "NA" if exact_mean is None else f"{exact_mean:.4f}"
        nonexact_mean_text = (
            "NA" if nonexact_mean is None else f"{nonexact_mean:.4f}"
        )
        print(
            "  confidence vs gene exact match: "
            f"AUROC/AUPRC={auroc_text}/{auprc_text}; "
            f"point-biserial r={correlation_text}; "
            f"mean exact/non-exact={exact_mean_text}/{nonexact_mean_text}; "
            f"exact/non-exact={metrics['confidence_exact_genes']}/"
            f"{metrics['confidence_nonexact_genes']}; "
            f"scored/missing={scored_genes}/"
            f"{metrics['confidence_missing_genes']}"
        )
    elif scored_genes == 0:
        print("  confidence vs gene exact match: unavailable (no scores)")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--pred",
        required=True,
        help="prediction prefix; reads <prefix>.gff and <prefix>.fna",
    )
    parser.add_argument(
        "--genome_dir",
        required=True,
        help=(
            "NCBI directory containing genomic.gff and cds_from_genomic.fna; "
            "optional sequence_report.jsonl enables nuclear-record selection"
        ),
    )
    parser.add_argument(
        "--codon-table",
        default=None,
        help=(
            "JSON translation table for predicted and reference CDS; defaults "
            "to codon_tables/<accession>.json when present, otherwise "
            "codon_tables/standard.json"
        ),
    )
    parser.add_argument(
        "--confidence-threshold",
        type=float,
        default=CONFIDENCE_THRESHOLD,
        help=(
            "gene-level confidence threshold: gene and protein metrics count "
            "predicted genes with confidence at or above this value "
            f"(default {CONFIDENCE_THRESHOLD}); genes without a confidence are "
            "always counted and 0 counts every gene"
        ),
    )
    parser.add_argument(
        "--output-json",
        help="atomically write the metrics to this JSON file",
    )
    return parser.parse_args()


def main():
    """Command-line entry point."""
    args = parse_args()
    if (
        not math.isfinite(args.confidence_threshold)
        or not 0.0 <= args.confidence_threshold <= 1.0
    ):
        raise SystemExit("--confidence-threshold must be between 0 and 1")

    metrics = evaluate_prediction(
        args.pred,
        args.genome_dir,
        confidence_threshold=args.confidence_threshold,
        codon_table_path=args.codon_table,
    )
    print_evaluation_metrics(metrics)
    if args.output_json:
        write_evaluation_metrics_json(args.output_json, metrics)
        print(f"[evaluation] wrote {args.output_json}")
    return metrics


if __name__ == "__main__":
    main()
