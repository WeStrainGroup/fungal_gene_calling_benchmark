"""Shared FASTA and frozen interval-evaluator helpers."""
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

REPO = Path(__file__).resolve().parents[1]
SPECIES = ('Saccharomyces_cerevisiae', 'Candida_albicans')
TOOLS = ('carbon', 'annevo', 'tiberius')
GROUPS = ('genome', '50-100kb', '5-10kb', '1-2kb')
sys.path.insert(0, str(REPO/'vendor/carbon_a'))
import evaluation as official


def fasta(path):
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, 'fasta')}


def sequence(model, genome):
    assert all(1 <= a <= b <= len(genome[model['record']]) for a, b in model['blocks'])
    seq = ''.join(genome[model['record']][a-1:b] for a, b in model['blocks'])
    return str(Seq(seq).reverse_complement()) if model['strand'] == '-' else seq
