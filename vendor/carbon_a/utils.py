"""Utility functions shared by the Carbon-A annotation pipeline.

The functions in this module do not import numba, torch, or Viterbi backends.
They are intentionally plain NumPy/Python helpers used by post-processing,
evaluation, and the CPU/GPU Viterbi backends.
"""

import json
import math
import operator
from pathlib import Path

import numpy as np

from constants import (
    AMINO_ACID_SYMBOLS,
    BASE_A,
    BASE_C,
    BASE_G,
    BASE_N,
    BASE_T,
    BASE_VALUES,
    CODON_BASES,
    CODON_TABLE_WARNING_ALPHA,
    CONFIDENCE_THRESHOLD,
    COMPLEMENT_BASES,
    FIXED_POINT_SCORE_SCALE,
    MIN_INTERNAL_INTRON_LENGTH_BP,
    NEGATIVE_STRAND_PROBABILITY_COLUMN,
    NON_NUCLEAR_ASSEMBLY_UNIT,
    POSITIVE_STRAND_PROBABILITY_COLUMN,
    RECORD_ID_COLUMN,
    SEQUENCE_COLUMN,
    START_CODONS,
    STOP_CODONS,
    TRANSITION_COUNT,
)

_BASE_ENCODING_VALUES = bytearray([BASE_N] * 256)
_COMPLEMENT_ENCODING_VALUES = bytearray([BASE_N] * 256)
for _base, _base_value in BASE_VALUES.items():
    _BASE_ENCODING_VALUES[ord(_base)] = _base_value
    _BASE_ENCODING_VALUES[ord(_base.lower())] = _base_value

    _complement_value = BASE_VALUES[COMPLEMENT_BASES[_base]]
    _COMPLEMENT_ENCODING_VALUES[ord(_base)] = _complement_value
    _COMPLEMENT_ENCODING_VALUES[ord(_base.lower())] = _complement_value

BASE_ENCODING_TABLE = bytes(_BASE_ENCODING_VALUES)
COMPLEMENT_ENCODING_TABLE = bytes(_COMPLEMENT_ENCODING_VALUES)
REVERSE_COMPLEMENT_TRANSLATION = str.maketrans("ATCGNatcgn", "TAGCNtagcn")

def validate_max_boundary_intron_length(value):
    """Return a boundary-intron limit of -1 or a nonnegative integer.

    ``-1`` leaves boundary introns unlimited, ``0`` disables both boundary
    intron states, and a positive value sets their maximum visible length in
    bases. With a limit of ``0``, incomplete genes can only result from a record
    boundary truncating an exon.
    """
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise TypeError("max_boundary_intron_length must be an integer") from exc
    if value < -1:
        raise ValueError(
            "max_boundary_intron_length must be -1 or a nonnegative integer"
        )
    return value


def validate_reconstruction_batch_size(value):
    """Return a reconstruction batch size in bases, or ``None`` for auto sizing."""
    if value is None:
        return None
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise TypeError("reconstruction_batch_size must be an integer or None") from exc
    if value <= 0:
        raise ValueError("reconstruction_batch_size must be positive or None")
    return value


def validate_memory_budget_bytes(value):
    """Return a device memory budget in bytes, or ``None`` for the default."""
    if value is None:
        return None
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise TypeError("memory_budget_bytes must be an integer or None") from exc
    if value <= 0:
        raise ValueError("memory_budget_bytes must be positive or None")
    return value


# Codon tables are flat JSON objects: all 64 concrete codons map to one-letter
# amino acid symbols, and "*" marks the stop codons Viterbi recognizes.
ALL_CODONS = frozenset(
    first + second + third
    for first in CODON_BASES
    for second in CODON_BASES
    for third in CODON_BASES
)
def validate_codon_table(codon_to_amino_acid):
    """Validate and return a copy of a complete codon-to-amino-acid mapping."""
    if not isinstance(codon_to_amino_acid, dict):
        raise ValueError("codon table JSON must contain one object")

    supplied_codons = set(codon_to_amino_acid)
    missing_codons = sorted(ALL_CODONS - supplied_codons)
    unexpected_codons = sorted(supplied_codons - ALL_CODONS)
    if missing_codons or unexpected_codons:
        details = []
        if missing_codons:
            details.append("missing " + ",".join(missing_codons))
        if unexpected_codons:
            details.append("unexpected " + ",".join(unexpected_codons))
        raise ValueError(
            "codon table must contain all 64 ATCG codons: "
            + "; ".join(details)
        )

    invalid_amino_acids = {
        codon: amino_acid
        for codon, amino_acid in codon_to_amino_acid.items()
        if not isinstance(amino_acid, str)
        or len(amino_acid) != 1
        or amino_acid not in AMINO_ACID_SYMBOLS
    }
    if invalid_amino_acids:
        formatted_values = ", ".join(
            f"{codon}={amino_acid!r}"
            for codon, amino_acid in sorted(invalid_amino_acids.items())
        )
        raise ValueError(
            "amino acids must use the 20 standard one-letter symbols or '*': "
            + formatted_values
        )

    unsupported_stop_codons = sorted(
        codon
        for codon, amino_acid in codon_to_amino_acid.items()
        if amino_acid == "*" and codon not in STOP_CODONS
    )
    if unsupported_stop_codons:
        raise ValueError(
            "codon table stop codons must be a subset of TAA,TAG,TGA: "
            + ",".join(unsupported_stop_codons)
        )

    return dict(codon_to_amino_acid)


def load_codon_table(path):
    """Load a codon-to-amino-acid mapping from a JSON file."""
    path = Path(path)
    try:
        with path.open(encoding="utf-8") as handle:
            codon_to_amino_acid = json.load(handle)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"invalid codon table JSON at {path}:{exc.lineno}:{exc.colno}: {exc.msg}"
        ) from exc
    return validate_codon_table(codon_to_amino_acid)


def stop_codons_from_table(codon_to_amino_acid):
    """Return supported stop codons whose amino-acid symbol is ``*``."""
    codon_to_amino_acid = validate_codon_table(codon_to_amino_acid)
    return tuple(
        codon
        for codon in STOP_CODONS
        if codon_to_amino_acid[codon] == "*"
    )


# Translation uses a five-value alphabet: A/T/C/G plus one value shared by all
# ambiguous bases. Three values form a base-5 codon index into this table.
TRANSLATION_BASE_COUNT = 5
AMBIGUOUS_TRANSLATION_BASE = BASE_N
TRANSLATION_BASE_TO_VALUE = np.full(
    256,
    AMBIGUOUS_TRANSLATION_BASE,
    dtype=np.uint8,
)
for _base, _base_value in BASE_VALUES.items():
    TRANSLATION_BASE_TO_VALUE[ord(_base)] = _base_value
    TRANSLATION_BASE_TO_VALUE[ord(_base.lower())] = _base_value


def build_translation_lookup(codon_to_amino_acid):
    """Build the compact lookup array used by vectorized protein translation."""
    codon_to_amino_acid = validate_codon_table(codon_to_amino_acid)
    lookup = np.full(
        TRANSLATION_BASE_COUNT ** 3,
        ord("X"),
        dtype=np.uint8,
    )
    for codon, amino_acid in codon_to_amino_acid.items():
        first_base, second_base, third_base = (
            BASE_VALUES[base] for base in codon
        )
        codon_index = (
            first_base * TRANSLATION_BASE_COUNT * TRANSLATION_BASE_COUNT
            + second_base * TRANSLATION_BASE_COUNT
            + third_base
        )
        lookup[codon_index] = ord(amino_acid)
    return lookup


# The Viterbi states (State design in viterbi.py).
PREFIX_TYPE_COUNT = 6
STATE_COUNT = 30

# Prefix type ids carried by exon and intron states. Codon phase is the number
# of CDS bases already consumed in the current codon. Prefix types retain the
# additional bases needed to recognize a standard stop on the next emission.
PHASE0 = 0               # codon phase 0; codon boundary, no open prefix
PHASE1_PREFIX_T = 1      # codon phase 1; one-base prefix T
PHASE1_PREFIX_OTHER = 2  # codon phase 1; one-base prefix other
PHASE2_PREFIX_TA = 3     # codon phase 2; prefix TA, so A/G completes TAA/TAG
PHASE2_PREFIX_TG = 4     # codon phase 2; prefix TG, so A completes TGA
PHASE2_PREFIX_OTHER = 5  # codon phase 2; other prefix, so no stop can be completed

# State ids, grouped by biological role. The tables below map each prefix type
# to its exon and intron states.
INTERGENIC = 0
START_A = 1
START_AT = 2
EXON_PHASE0 = 3
EXON_PHASE1_PREFIX_T = 4
EXON_PHASE1_PREFIX_OTHER = 5
EXON_PHASE2_PREFIX_TA = 6
EXON_PHASE2_PREFIX_TG = 7
EXON_PHASE2_PREFIX_OTHER = 8
STOP = 9
GT_GC_AG_INTRON_FIRST_PHASE0 = 10
GT_GC_AG_INTRON_FIRST_PHASE1_PREFIX_T = 11
GT_GC_AG_INTRON_FIRST_PHASE1_PREFIX_OTHER = 12
GT_GC_AG_INTRON_FIRST_PHASE2_PREFIX_TA = 13
GT_GC_AG_INTRON_FIRST_PHASE2_PREFIX_TG = 14
GT_GC_AG_INTRON_FIRST_PHASE2_PREFIX_OTHER = 15
GT_GC_AG_INTRON_READY_PHASE0 = 16
GT_GC_AG_INTRON_READY_PHASE1_PREFIX_T = 17
GT_GC_AG_INTRON_READY_PHASE1_PREFIX_OTHER = 18
GT_GC_AG_INTRON_READY_PHASE2_PREFIX_TA = 19
GT_GC_AG_INTRON_READY_PHASE2_PREFIX_TG = 20
GT_GC_AG_INTRON_READY_PHASE2_PREFIX_OTHER = 21
AT_AC_INTRON_PHASE0 = 22
AT_AC_INTRON_PHASE1_PREFIX_T = 23
AT_AC_INTRON_PHASE1_PREFIX_OTHER = 24
AT_AC_INTRON_PHASE2_PREFIX_TA = 25
AT_AC_INTRON_PHASE2_PREFIX_TG = 26
AT_AC_INTRON_PHASE2_PREFIX_OTHER = 27
LEADING_BOUNDARY_INTRON = 28
TRAILING_BOUNDARY_INTRON = 29

EXON_STATE_BY_PREFIX_TYPE = np.array((
    EXON_PHASE0,
    EXON_PHASE1_PREFIX_T,
    EXON_PHASE1_PREFIX_OTHER,
    EXON_PHASE2_PREFIX_TA,
    EXON_PHASE2_PREFIX_TG,
    EXON_PHASE2_PREFIX_OTHER,
), dtype=np.int64)

GT_GC_AG_INTRON_FIRST_STATE_BY_PREFIX_TYPE = np.array((
    GT_GC_AG_INTRON_FIRST_PHASE0,
    GT_GC_AG_INTRON_FIRST_PHASE1_PREFIX_T,
    GT_GC_AG_INTRON_FIRST_PHASE1_PREFIX_OTHER,
    GT_GC_AG_INTRON_FIRST_PHASE2_PREFIX_TA,
    GT_GC_AG_INTRON_FIRST_PHASE2_PREFIX_TG,
    GT_GC_AG_INTRON_FIRST_PHASE2_PREFIX_OTHER,
), dtype=np.int64)

GT_GC_AG_INTRON_READY_STATE_BY_PREFIX_TYPE = np.array((
    GT_GC_AG_INTRON_READY_PHASE0,
    GT_GC_AG_INTRON_READY_PHASE1_PREFIX_T,
    GT_GC_AG_INTRON_READY_PHASE1_PREFIX_OTHER,
    GT_GC_AG_INTRON_READY_PHASE2_PREFIX_TA,
    GT_GC_AG_INTRON_READY_PHASE2_PREFIX_TG,
    GT_GC_AG_INTRON_READY_PHASE2_PREFIX_OTHER,
), dtype=np.int64)

AT_AC_INTRON_STATE_BY_PREFIX_TYPE = np.array((
    AT_AC_INTRON_PHASE0,
    AT_AC_INTRON_PHASE1_PREFIX_T,
    AT_AC_INTRON_PHASE1_PREFIX_OTHER,
    AT_AC_INTRON_PHASE2_PREFIX_TA,
    AT_AC_INTRON_PHASE2_PREFIX_TG,
    AT_AC_INTRON_PHASE2_PREFIX_OTHER,
), dtype=np.int64)


# Unlike codon phase above, GFF3 phase is the number of visible bases to skip at
# the beginning of a CDS feature before reaching the first complete codon. A
# boundary-truncated gene's first state retains enough information to recover it.
INITIAL_GFF3_PHASE_BY_STATE = {
    START_A: 0,
    START_AT: 2,
    EXON_PHASE0: 1,
    EXON_PHASE1_PREFIX_T: 0,
    EXON_PHASE1_PREFIX_OTHER: 0,
    EXON_PHASE2_PREFIX_TA: 2,
    EXON_PHASE2_PREFIX_TG: 2,
    EXON_PHASE2_PREFIX_OTHER: 2,
    STOP: 1,
}


def initial_gff3_phase_from_state(state):
    """Return the initial GFF3 phase encoded by the first visible CDS state."""
    try:
        return INITIAL_GFF3_PHASE_BY_STATE[int(state)]
    except KeyError as exc:
        raise ValueError(f"state {state} does not emit CDS") from exc


def reverse_complement(sequence):
    """Return the reverse complement string for a DNA sequence."""
    return sequence.translate(REVERSE_COMPLEMENT_TRANSLATION)[::-1]


def encode_both_strands(sequence):
    """Encode both strands as A/T/C/G/N = 0/1/2/3/4."""
    sequence_bytes = sequence.encode("ascii")
    forward_bases = np.frombuffer(
        sequence_bytes.translate(BASE_ENCODING_TABLE),
        dtype=np.int8,
    )
    reverse_bases = np.frombuffer(
        sequence_bytes.translate(COMPLEMENT_ENCODING_TABLE)[::-1],
        dtype=np.int8,
    )
    return forward_bases, reverse_bases


def probability_vector_from_arrow(value):
    """Convert one PyArrow list cell to a 1D float16 probability array."""
    if not hasattr(value, "values"):
        raise TypeError("probability column must contain PyArrow list values")
    probabilities = np.asarray(value.values, dtype=np.float16)

    if probabilities.ndim != 1:
        raise ValueError(f"probability vector must be 1D, got {probabilities.shape}")
    return probabilities


def prepare_decode_tasks(tasks):
    """Normalize Viterbi tasks and reject invalid shapes or lengths."""
    prepared_tasks = []
    for task_index, (probabilities, encoded_bases) in enumerate(tasks):
        probability_array = np.ascontiguousarray(probabilities, dtype=np.float16)
        base_array = np.ascontiguousarray(encoded_bases, dtype=np.int8)
        if probability_array.ndim != 1:
            raise ValueError(
                f"Viterbi task {task_index} probabilities must be 1D, "
                f"got {probability_array.shape}"
            )
        if base_array.ndim != 1:
            raise ValueError(
                f"Viterbi task {task_index} encoded bases must be 1D, "
                f"got {base_array.shape}"
            )
        if len(probability_array) != len(base_array):
            raise ValueError(
                f"Viterbi task {task_index} has {len(probability_array)} "
                f"probabilities but {len(base_array)} encoded bases"
            )
        if len(probability_array) == 0:
            raise ValueError(f"Viterbi task {task_index} is empty")
        prepared_tasks.append((probability_array, base_array))
    return prepared_tasks


def pack_decode_tasks(prepared_tasks, prepared_transition_log_odds):
    """Return prepared Viterbi tasks as one batch for a native backend.

    The batch is the P(CDS) bits (uint16), the encoded bases, the int32
    transition scores, and the record offsets, one more than the tasks. A
    single task is passed on without a copy.
    """
    record_offsets = np.zeros(len(prepared_tasks) + 1, np.int64)
    np.cumsum([len(encoded_bases) for _, encoded_bases in prepared_tasks], out=record_offsets[1:])
    if len(prepared_tasks) == 1:
        [(probabilities, encoded_bases)] = prepared_tasks
        [transition_scores] = prepared_transition_log_odds
    else:
        probabilities = np.concatenate([probabilities for probabilities, _ in prepared_tasks])
        encoded_bases = np.concatenate([encoded_bases for _, encoded_bases in prepared_tasks])
        transition_scores = np.concatenate(prepared_transition_log_odds)
    return probabilities.view(np.uint16), encoded_bases, transition_scores, record_offsets


def quantize_transition_log_odds(log_odds):
    """Return transition log-odds in nats as the Viterbi backends' int32 fixed-point scores.

    Every value is multiplied by ``FIXED_POINT_SCORE_SCALE`` and rounded to
    the nearest integer. float16 and float32 values are scaled in float32,
    other values in float64; the scale is a power of two, so the product is
    exact and rounds to the same integer in either precision.
    """
    log_odds = np.asarray(log_odds)
    exact_dtype = np.float32 if log_odds.dtype in (np.float16, np.float32) else np.float64
    scaled = np.rint(log_odds.astype(exact_dtype, copy=False) * exact_dtype(FIXED_POINT_SCORE_SCALE))
    limits = np.iinfo(np.int32)
    # As Python floats the extremes compare exactly with the int32 limits; NaN fails.
    if scaled.size and not limits.min <= float(scaled.min()) <= float(scaled.max()) <= limits.max:
        raise OverflowError("transition log-odds exceed the int32 fixed-point range")
    return scaled.astype(np.int32)


def prepare_transition_log_odds(transition_log_odds, prepared_tasks):
    """Check one ``(bases, 4)`` int32 transition log-odds array per task for the Viterbi backends.

    The columns follow ``TRANSITION_NAMES``: enter gene, leave gene, enter
    intron, leave intron. Each value is the transition network's log-odds
    against staying in ``FIXED_POINT_SCORE_SCALE`` units
    (``quantize_transition_log_odds``), aligned base by base with the task's
    probabilities. The arrays come back C-contiguous.
    """
    transition_log_odds = list(transition_log_odds)
    if len(transition_log_odds) != len(prepared_tasks):
        raise ValueError(
            f"{len(transition_log_odds)} transition log-odds arrays for "
            f"{len(prepared_tasks)} Viterbi tasks"
        )
    prepared_log_odds = []
    for task_index, (log_odds, (probabilities, _)) in enumerate(
        zip(transition_log_odds, prepared_tasks)
    ):
        log_odds = np.asarray(log_odds)
        if log_odds.shape != (len(probabilities), TRANSITION_COUNT):
            raise ValueError(
                f"Viterbi task {task_index} transition log-odds must have shape "
                f"{(len(probabilities), TRANSITION_COUNT)}, got {log_odds.shape}"
            )
        if log_odds.dtype != np.int32:
            raise TypeError(
                f"Viterbi task {task_index} transition log-odds must be int32 fixed-point scores, "
                f"got {log_odds.dtype}"
            )
        prepared_log_odds.append(np.ascontiguousarray(log_odds))
    return prepared_log_odds


def build_cds_sequence(genome_sequence, transcript_order_exons, strand):
    """Build a CDS string from genomic exon intervals in transcript order."""
    if not transcript_order_exons:
        return ""
    if strand == "+":
        return "".join(
            genome_sequence[start:end + 1] for start, end in transcript_order_exons
        )

    genomic_order_sequence = "".join(
        genome_sequence[start:end + 1]
        for start, end in reversed(transcript_order_exons)
    )
    return reverse_complement(genomic_order_sequence)


def compute_gff_cds_phases(transcript_order_exons, initial_gff3_phase=0):
    """Return GFF3 phases, preserving the first visible codon's phase."""
    if initial_gff3_phase not in (0, 1, 2):
        raise ValueError("initial GFF3 phase must be 0, 1, or 2")

    phases = []
    consumed_bases = 0
    for start, end in transcript_order_exons:
        phase = (initial_gff3_phase - consumed_bases) % 3
        phases.append(phase)
        consumed_bases += end - start + 1
    return phases


def has_internal_stop_codon(
    cds_sequence,
    translation_frame,
    stop_codons=STOP_CODONS,
):
    """Return whether a translated CDS contains a non-terminal stop codon.

    This is a diagnostic check: the Viterbi grammar should normally make the
    result false. A stop is terminal only when its third base is the final CDS
    base; a stop before an incomplete trailing codon is still internal.
    """
    if translation_frame not in (0, 1, 2):
        raise ValueError("translation frame must be 0, 1, or 2")

    codon_end = len(cds_sequence) - (
        (len(cds_sequence) - translation_frame) % 3
    )
    for codon_start in range(translation_frame, codon_end, 3):
        codon = cds_sequence[codon_start:codon_start + 3].upper()
        if codon in stop_codons and codon_start + 3 < len(cds_sequence):
            return True
    return False


def classify_cds_completion(
    cds_sequence,
    has_start_codon,
    has_stop_codon,
    has_short_intron=False,
    stop_codons=STOP_CODONS,
):
    """Return ``(complete, reason)`` for a reconstructed CDS sequence.

    Reasons are comma-separated so boundary-truncated predictions can state
    every visible problem inside one GFF3 attribute, for example
    ``no_start,no_stop``. ``short_intron``
    records a Viterbi path with fewer than four bases between adjacent exons.
    Start and stop presence comes from the decoded states, then the sequence
    is checked as a consistency guard.
    """
    reasons = []
    if len(cds_sequence) < 6:
        reasons.append("too_short")
    if len(cds_sequence) % 3 != 0:
        reasons.append("not_multiple_of_3")
    if has_short_intron:
        reasons.append("short_intron")
    if (
        not has_start_codon
        or len(cds_sequence) < 3
        or cds_sequence[:3] not in START_CODONS
    ):
        reasons.append("no_start")
    if (
        not has_stop_codon
        or len(cds_sequence) < 3
        or cds_sequence[-3:] not in stop_codons
    ):
        reasons.append("no_stop")

    if reasons:
        # Comma-separated: a semicolon would end the GFF3 attribute.
        return False, ",".join(reasons)
    return True, None


def binomial_cdf(k, n):
    """Return P(X <= k) for X ~ Binomial(n, 1/2).

    Below the median the terms shrink going down from ``k``, so the sum runs
    downward in log space and stops once a term no longer changes it. Above
    the median the symmetry of Binomial(n, 1/2) gives
    P(X <= k) = 1 - P(X <= n - k - 1).
    """
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if 2 * k >= n:
        return 1.0 - binomial_cdf(n - k - 1, n)

    def log_probability(count):
        return (
            math.lgamma(n + 1) - math.lgamma(count + 1)
            - math.lgamma(n - count + 1) - n * math.log(2)
        )

    largest = log_probability(k)
    relative_sum = 0.0
    for count in range(k, -1, -1):
        term = math.exp(log_probability(count) - largest)
        relative_sum += term
        if term < 1e-17 * relative_sum:
            break
    return math.exp(largest) * relative_sum


def codon_table_mismatch(
    genes, stop_codons, threshold=CONFIDENCE_THRESHOLD, alpha=CODON_TABLE_WARNING_ALPHA,
):
    """Return the codon groups of the complete genes that contradict the codon table.

    A codon table that does not match the genome splits genes at codons the
    genome reads through or fuses genes at stop codons it lacks, and the
    confidence network gives such genes almost no confidence. The scored
    complete genes (``confidence`` not None) form one group per stop codon of
    the table, the genes ending in it, and one per TAA, TAG or TGA the table
    reads as an amino acid, the genes containing it in frame. A group of ``n``
    genes, ``k`` of them reaching ``threshold``, contradicts the table when
    significantly fewer than half reach it: P(X <= k) < ``alpha`` for
    X ~ Binomial(n, 1/2), a one-sided binomial test. Returns
    ``(group, k, n)`` for every contradicting group, for example
    ``("ending in TAG", 65, 7877)``; the list is empty when no group
    contradicts the table.
    """
    read_through_codons = [codon for codon in STOP_CODONS if codon not in stop_codons]
    counts = {f"ending in {codon}": [0, 0] for codon in stop_codons}
    counts.update({f"with {codon} in frame": [0, 0] for codon in read_through_codons})
    for gene in genes:
        if gene["confidence"] is None:
            continue
        cds = gene["cds"]
        groups = [f"ending in {cds[-3:]}"]
        if read_through_codons:
            in_frame_codons = {cds[start:start + 3] for start in range(0, len(cds) - 3, 3)}
            groups += [f"with {codon} in frame" for codon in read_through_codons if codon in in_frame_codons]
        # A decoded complete gene always ends in a stop codon of the table; a
        # gene model that does not has no stop-codon group.
        for group in groups:
            if group in counts:
                counts[group][0] += int(gene["confidence"] >= threshold)
                counts[group][1] += 1
    return [
        (group, confident_count, gene_count)
        for group, (confident_count, gene_count) in counts.items()
        if gene_count and binomial_cdf(confident_count, gene_count) < alpha
    ]


def has_short_internal_intron(exons):
    """Return whether adjacent exons leave an internal intron under four bases."""
    for exon_index in range(1, len(exons)):
        previous_exon = exons[exon_index - 1]
        next_exon = exons[exon_index]
        if previous_exon[0] < next_exon[0]:
            intron_length = next_exon[0] - previous_exon[1] - 1
        else:
            intron_length = previous_exon[0] - next_exon[1] - 1
        if 0 < intron_length < MIN_INTERNAL_INTRON_LENGTH_BP:
            return True
    return False


def translate_cds_sequence(
    cds_sequence,
    frame,
    amino_acid_by_codon_index,
):
    """Translate one CDS sequence with a prepared codon-table lookup.

    A single stop at the final complete codon is removed from the protein.
    Any earlier stop indicates an invalid decoded CDS and raises an assertion
    instead of silently truncating the output. An incomplete trailing codon is
    ignored unless it follows a stop, in which case that stop is internal.
    """
    codon_end = len(cds_sequence) - ((len(cds_sequence) - frame) % 3)
    if codon_end <= frame:
        return ""

    sequence_bytes = np.frombuffer(cds_sequence.encode("ascii"), dtype=np.uint8)
    base_values = TRANSLATION_BASE_TO_VALUE[
        sequence_bytes[frame:codon_end]
    ].reshape(-1, 3)
    codon_indices = (
        base_values[:, 0] * TRANSLATION_BASE_COUNT * TRANSLATION_BASE_COUNT
        + base_values[:, 1] * TRANSLATION_BASE_COUNT
        + base_values[:, 2]
    )
    amino_acids = amino_acid_by_codon_index[codon_indices]

    stop_positions = np.flatnonzero(amino_acids == ord("*"))
    if stop_positions.size:
        only_terminal_stop = (
            stop_positions.size == 1
            and stop_positions[0] == len(amino_acids) - 1
            and codon_end == len(cds_sequence)
        )
        if not only_terminal_stop:
            stop_index = int(stop_positions[0])
            codon_start = frame + stop_index * 3
            stop_codon = cds_sequence[codon_start:codon_start + 3].upper()
            raise AssertionError(
                "protein translation encountered an internal stop codon: "
                f"codon={stop_codon}, cds_offset={codon_start}"
            )
        amino_acids = amino_acids[:-1]
    return amino_acids.tobytes().decode("ascii")


def assemble_gene_model(
    record_id,
    genome_sequence,
    oriented_exons,
    strand,
    record_length,
    first_cds_state,
    last_cds_state,
    stop_codons=STOP_CODONS,
):
    """Build one output gene from its exons and boundary CDS states."""
    initial_gff3_phase = initial_gff3_phase_from_state(first_cds_state)
    if strand == "+":
        transcript_order_exons = [
            (int(start), int(end)) for start, end in oriented_exons
        ]
    else:
        transcript_order_exons = [
            (record_length - 1 - int(end), record_length - 1 - int(start))
            for start, end in oriented_exons
        ]

    cds_sequence = build_cds_sequence(genome_sequence, transcript_order_exons, strand)
    if not cds_sequence:
        return None

    transcript_order_phases = compute_gff_cds_phases(
        transcript_order_exons,
        initial_gff3_phase,
    )
    genomic_order = sorted(
        zip(transcript_order_exons, transcript_order_phases),
        key=lambda exon_and_phase: exon_and_phase[0][0],
    )

    assert not has_internal_stop_codon(
        cds_sequence,
        initial_gff3_phase,
        stop_codons,
    ), (
        "Viterbi produced an internal stop codon: "
        f"record_id={record_id!r}, strand={strand}, "
        f"phase={initial_gff3_phase}, exons={transcript_order_exons}"
    )
    complete, incomplete_reason = classify_cds_completion(
        cds_sequence,
        first_cds_state == START_A,
        last_cds_state == STOP,
        has_short_intron=has_short_internal_intron(transcript_order_exons),
        stop_codons=stop_codons,
    )
    return {
        "record_id": record_id,
        "strand": strand,
        "cds": cds_sequence,
        "translation_frame": initial_gff3_phase,
        "complete": complete,
        "incomplete_reason": incomplete_reason,
        "exons": [exon for exon, phase in genomic_order],
        "phases": [phase for exon, phase in genomic_order],
    }


def format_gene_status(gene, separator):
    """Format shared complete/incomplete metadata for GFF3 or FASTA."""
    status = "complete" if gene["complete"] else "incomplete"
    fields = [f"status={status}"]
    if not gene["complete"]:
        reason = gene.get("incomplete_reason") or "unknown"
        fields.append(f"incomplete_reason={reason}")
    return separator.join(fields)


def format_gene_confidence(gene):
    """Format the per-gene confidence for GFF3 attributes or FASTA headers.

    Complete genes carry the confidence model's probability that the decoded
    gene is exactly right. Incomplete genes are not scored and are written as
    ``NA`` so that a confidence filter keeps them.
    """
    value = gene.get("confidence")
    formatted_value = "NA" if value is None else f"{value:.6f}"
    return f"confidence={formatted_value}"


def write_gff3(genes, output_path, source="Carbon-A"):
    """Write gene, mRNA, and CDS features in GFF3 format."""
    with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        for gene_number, gene in enumerate(genes, 1):
            gene_id = f"gene{gene_number:06d}"
            transcript_id = f"mrna{gene_number:06d}"
            gene_start = gene["exons"][0][0] + 1
            gene_end = gene["exons"][-1][1] + 1
            status_attributes = format_gene_status(gene, separator=";")
            status_attributes += ";" + format_gene_confidence(gene)

            handle.write(
                f"{gene['record_id']}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t"
                f"{gene['strand']}\t.\tID={gene_id};{status_attributes}\n"
            )
            handle.write(
                f"{gene['record_id']}\t{source}\tmRNA\t{gene_start}\t{gene_end}\t.\t"
                f"{gene['strand']}\t.\tID={transcript_id};Parent={gene_id}\n"
            )

            for (start, end), phase in zip(gene["exons"], gene["phases"]):
                handle.write(
                    f"{gene['record_id']}\t{source}\tCDS\t{start + 1}\t{end + 1}\t.\t"
                    f"{gene['strand']}\t{phase}\tID={transcript_id}.cds;Parent={transcript_id}\n"
                )


def write_cds_fasta(genes, output_path):
    """Write predicted CDS sequences in FASTA format."""
    with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
        for gene_number, gene in enumerate(genes, 1):
            gene_start = gene["exons"][0][0] + 1
            gene_end = gene["exons"][-1][1] + 1
            header = (
                f">gene{gene_number:06d} record_id={gene['record_id']} "
                f"strand={gene['strand']} loc={gene_start}-{gene_end} "
                f"phase={gene['translation_frame']} "
                f"{format_gene_status(gene, separator=' ')}"
            )
            header += " " + format_gene_confidence(gene)

            handle.write(header + "\n")
            cds_sequence = gene["cds"]
            for line_start in range(0, len(cds_sequence), 80):
                handle.write(cds_sequence[line_start:line_start + 80] + "\n")


def write_protein_fasta(
    genes,
    output_path,
    codon_to_amino_acid,
):
    """Write translated protein sequences in FASTA format.

    Returns the number of non-empty protein records written. Incomplete CDS
    records keep their reason in the header so downstream filters can decide
    whether to use them.
    """
    amino_acid_by_codon_index = build_translation_lookup(codon_to_amino_acid)
    written_count = 0
    with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
        for gene_number, gene in enumerate(genes, 1):
            cds_sequence = gene["cds"]
            protein_sequence = translate_cds_sequence(
                cds_sequence,
                gene["translation_frame"],
                amino_acid_by_codon_index,
            )
            if not protein_sequence:
                continue

            gene_start = gene["exons"][0][0] + 1
            gene_end = gene["exons"][-1][1] + 1
            header = (
                f">gene{gene_number:06d} record_id={gene['record_id']} "
                f"strand={gene['strand']} loc={gene_start}-{gene_end} "
                f"cds_len={len(cds_sequence)} aa_len={len(protein_sequence)} "
                f"{format_gene_status(gene, separator=' ')}"
            )
            header += " " + format_gene_confidence(gene)

            handle.write(header + "\n")
            for line_start in range(0, len(protein_sequence), 80):
                handle.write(protein_sequence[line_start:line_start + 80] + "\n")
            written_count += 1

    return written_count


def load_nuclear_evaluation_accessions(genome_dir):
    """Return the record accessions used for nuclear-gene evaluation.

    Prefer assembled chromosomes. If none exist, use assembled linkage groups;
    if those are also absent, use chromosome scaffolds. ``None`` means that no
    sequence report is available for filtering.
    """
    sequence_report_path = Path(genome_dir) / "sequence_report.jsonl"
    if not sequence_report_path.exists():
        return None

    assembled_chromosome_accessions = set()
    assembled_linkage_group_accessions = set()
    chromosome_scaffold_accessions = set()
    with sequence_report_path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            row = json.loads(line)
            accession = row.get("refseqAccession") or row.get("genbankAccession")
            if not accession:
                continue

            assembly_unit = row.get("assemblyUnit")
            location_type = row.get("assignedMoleculeLocationType")
            if assembly_unit == NON_NUCLEAR_ASSEMBLY_UNIT:
                continue

            role = row.get("role")
            if role == "assembled-molecule":
                if location_type == "Chromosome":
                    assembled_chromosome_accessions.add(accession)
                elif location_type == "Linkage Group":
                    assembled_linkage_group_accessions.add(accession)
            elif location_type == "Chromosome" and (role or "").endswith("-scaffold"):
                chromosome_scaffold_accessions.add(accession)

    return (
        assembled_chromosome_accessions
        or assembled_linkage_group_accessions
        or chromosome_scaffold_accessions
    )
