"""Shared constants for the Carbon-A annotation pipeline.

Keep cross-module data contracts and tunable pipeline parameters here. Private
implementation details, such as Viterbi state ids and evaluation regexes, stay
with the modules that own them.
"""

import math
from pathlib import Path


PROJECT_DIRECTORY = Path(__file__).resolve().parent.parent


# Model inference and output schema. The model reads windows of
# CONTEXT_LENGTH tokens (six bases each); consecutive windows of a long record
# overlap by OVERLAP_LENGTH tokens and overlapping probabilities are averaged.
MODEL_NAME = "HuggingFaceBio/Carbon-A-1.2B"
CONTEXT_LENGTH = 16384
OVERLAP_LENGTH = 1024
HEAD_NAMES = ("positive_strand", "negative_strand")
RECORD_ID_COLUMN = "record_name"
SEQUENCE_COLUMN = "sequence"
POSITIVE_STRAND_PROBABILITY_COLUMN = "pred_prob_positive_strand_cds"
NEGATIVE_STRAND_PROBABILITY_COLUMN = "pred_prob_negative_strand_cds"


# Sequence alphabet and biological constraints.
BASE_A = 0
BASE_T = 1
BASE_C = 2
BASE_G = 3
BASE_N = 4
BASE_VALUES = {
    "A": BASE_A,
    "T": BASE_T,
    "C": BASE_C,
    "G": BASE_G,
    "N": BASE_N,
}
COMPLEMENT_BASES = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
}
START_CODONS = ("ATG",)
STOP_CODONS = ("TAA", "TAG", "TGA")
CODON_BASES = "ATCG"
AMINO_ACID_SYMBOLS = frozenset("ACDEFGHIKLMNPQRSTVWY*")
NON_NUCLEAR_ASSEMBLY_UNIT = "non-nuclear"


# Viterbi and post-processing parameters.
# Keep emission probabilities finite before taking logarithms.
MIN_PROB = 1.0e-6
MAX_PROB = 1.0 - MIN_PROB

# Output gene models must contain at least 32 codons. A complete donor and
# acceptor occupy four distinct intron bases.
MIN_CDS_LENGTH_BP = 96
MIN_INTERNAL_INTRON_LENGTH_BP = 4

# Codon-table check. A codon table that does not match the genome splits genes
# at codons the genome reads through or fuses genes at stop codons it lacks,
# and the confidence network gives such genes almost no confidence. The check
# groups the complete genes by codon: one group per stop codon of the table
# (the genes ending in it) and one per TAA, TAG or TGA the table reads as an
# amino acid (the genes containing it in frame). It warns when, in any group,
# significantly fewer than half of the genes reach CONFIDENCE_THRESHOLD: k of n
# genes with P(X <= k) < CODON_TABLE_WARNING_ALPHA for X ~ Binomial(n, 1/2), a
# one-sided binomial test.
CODON_TABLE_WARNING_ALPHA = 0.01
# A failed check is recorded where a batch run will see it: appended to this
# file in the parent of the output directory, and in the output directory's
# name (run_*.sh share the suffix).
WARNINGS_FILE_NAME = "WARNINGS.txt"
WARNING_DIRECTORY_SUFFIX = "_WARNING"

# Intron classes are scored relative to the 99% GT-AG baseline.
GT_AG_INTRON_LOG_PRIOR = math.log(0.990 / 0.99)
GC_AG_INTRON_LOG_PRIOR = math.log(0.009 / 0.99)
AT_AC_INTRON_LOG_PRIOR = math.log(0.001 / 0.99)

# The two networks of the annotation pipeline. The transition network reads
# the coding-probability track and scores the four permitted transitions at
# every base; the confidence network reads the track together with the
# decoded roles of the whole strand and scores every complete gene. Both
# read coding probabilities clipped to NETWORK_INPUT_MIN_PROB as log-odds.
WEIGHTS_DIRECTORY = PROJECT_DIRECTORY / "weights"
TRANSITION_NETWORK_PATH = WEIGHTS_DIRECTORY / "transition_network.safetensors"
CONFIDENCE_NETWORK_PATH = WEIGHTS_DIRECTORY / "confidence_network.safetensors"
# Per-base transition tracks handed to the Viterbi backends, in this column
# order. They are stored as log-odds against staying in the backends'
# fixed-point units; the native backends and the Numba reference share the
# layout.
TRANSITION_NAMES = ("enter_gene", "leave_gene", "enter_intron", "leave_intron")
TRANSITION_COUNT = len(TRANSITION_NAMES)
ENTER_GENE, LEAVE_GENE, ENTER_INTRON, LEAVE_INTRON = range(TRANSITION_COUNT)
NETWORK_INPUT_MIN_PROB = 1.0e-4
# Both networks read long tracks in chunks of this many bases, each evaluated
# with the network's full context on both sides (tcn.chunked_inference), so a
# chromosome never has to fit in GPU memory at once. (The post-processing
# driver separately groups whole strands into host-side chunks, chunk_mbp.)
NETWORK_CHUNK_BP = 262_144
# Device memory the shipped confidence network needs per exon base of the
# genes it pools at once. The peak is in the attention MLP: the base's 32
# float32 features, the 32 float32 hidden values before and the 32 after GELU,
# and the base's int64 gene index. The genes of a strand are pooled in groups
# whose exon bases fit GPU_MEMORY_FRACTION of the free device memory at this
# size.
CONFIDENCE_NETWORK_BYTES_PER_EXON_BASE = 392
# Native Viterbi implementations accumulate int32 fixed-point scores.
# Transition log-odds are exchanged in the same units: nats multiplied by
# this scale and rounded.
FIXED_POINT_SCORE_SCALE = 16384.0

# Boundary introns model genes truncated at a record edge. A finite limit
# prevents an incomplete path from hiding a long region as one intron.
MAX_BOUNDARY_INTRON_LENGTH_BP = 20_000

# Records shorter than this are skipped.
MIN_RECORD_LENGTH_BP = 100_000
# Host-side chunk size in megabases by Viterbi backend. Post-processing takes
# about 50 MB of host memory per megabase of its largest chunk. The CPU
# backend decodes the strands of a chunk in parallel, one thread each, so it
# is given large chunks; the other backends gain nothing from them. A strand
# is never split: one longer than the chunk size is a chunk of its own.
CHUNK_SIZE_MBP = {"gpu": 25, "cpu": 1000, "numba": 25}
PARQUET_BATCH_ROWS = 1
PARQUET_READ_BUFFER_BYTES = 8 << 20
CODON_TABLE_PATH = PROJECT_DIRECTORY / "codon_tables" / "standard.json"


# GPU Viterbi backend parameters.
# Bases per segment; a batch whose transition log-odds are too large for this
# length to keep segment scores within int32 gets the longest length that does.
GPU_SEGMENT_LENGTH = 300
# Segments per group in stage 3 of the GPU backend; 0 takes the square root of
# the longest record's segment count, rounded up.
GPU_SEGMENT_GROUP_LENGTH = 0
# The device memory budget of a GPU Viterbi call and of the confidence
# network's pooling, as a fraction of the free device memory.
GPU_MEMORY_FRACTION = 0.90


# Gene-level results count the predicted genes whose confidence reaches this
# threshold; genes without a confidence are always counted.
CONFIDENCE_THRESHOLD = 0.1
