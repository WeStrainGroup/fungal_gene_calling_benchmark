"""Regression checks for restoring exact frozen input bytes."""
import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from reconstruct_fragments import reconstruct


class ReconstructTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.data = Path(self.temp.name)
        self.species = self.data / 'inputs' / 'fungus'
        self.species.mkdir(parents=True)
        self.genome = self.species / 'genome.fna'
        self.genome.write_text('>chr description\nACGT\nNNTA\n')
        self.mapping = self.species / '1-2kb.fragments.tsv'
        self.output = self.species / '1-2kb.fna'
        self.expected = b'>second\nNNTA\n>first\nACGT\n'
        with self.mapping.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['record', 'chrom', 'start', 'end', 'bp', 'sequence_sha256'])
            for name, start, end, seq in [('second', 5, 8, b'NNTA'), ('first', 1, 4, b'ACGT')]:
                writer.writerow([name, 'chr', start, end, len(seq), hashlib.sha256(seq).hexdigest()])
        self.manifest()

    def manifest(self, output_bytes=None):
        entries = {self.genome: self.genome.read_bytes(), self.mapping: self.mapping.read_bytes(),
                   self.output: self.expected if output_bytes is None else output_bytes}
        (self.data / 'SHA256SUMS').write_text(''.join(
            f'{hashlib.sha256(content).hexdigest()}  {path.relative_to(self.data)}\n'
            for path, content in entries.items()))

    def test_exact_bytes_boundaries_order_and_verified_skip(self):
        self.assertEqual(reconstruct(self.data), (1, 1))
        self.assertEqual(self.output.read_bytes(), self.expected)
        self.assertEqual(reconstruct(self.data), (1, 0))

    def test_changed_original_and_map_rejected(self):
        for path in (self.genome, self.mapping):
            original = path.read_bytes()
            path.write_bytes(original + b'\n')
            with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                reconstruct(self.data)
            path.write_bytes(original)
        self.assertFalse(self.output.exists())

    def test_wrong_fragment_sequence_rejected(self):
        self.mapping.write_text(self.mapping.read_text().replace('5\t8', '4\t7'))
        self.manifest()
        with self.assertRaisesRegex(ValueError, 'Sequence SHA256 mismatch'):
            reconstruct(self.data)
        self.assertFalse(self.output.exists())
        self.assertFalse(self.output.with_suffix('.fna.tmp').exists())

    def test_duplicate_or_out_of_bounds_record_rejected(self):
        original = self.mapping.read_text()
        for text in (original.replace('first\t', 'second\t'), original.replace('5\t8', '0\t3')):
            self.mapping.write_text(text)
            self.manifest()
            with self.assertRaisesRegex(ValueError, 'Invalid record or coordinates'):
                reconstruct(self.data)

    def test_bad_existing_file_not_overwritten(self):
        self.output.write_bytes(b'bad')
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            reconstruct(self.data)
        self.assertEqual(self.output.read_bytes(), b'bad')

    def test_whole_fasta_hash_must_match(self):
        self.manifest(self.expected.replace(b'\n', b'\r\n'))
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            reconstruct(self.data)
        self.assertFalse(self.output.exists())
        self.assertFalse(self.output.with_suffix('.fna.tmp').exists())


if __name__ == '__main__':
    unittest.main()
