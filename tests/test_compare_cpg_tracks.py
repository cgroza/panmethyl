import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "module/resources/usr/bin/compare_cpg_tracks.py"


class CompareCpgTracksTests(unittest.TestCase):
    def run_comparison(self, left, right, directory):
        left_file = directory / "left.bed"
        left_file.write_text(left)
        right_file = directory / "right.bed.gz"
        data = right.encode()
        right_file.write_bytes(gzip.compress(data[:20]) + gzip.compress(data[20:]))
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(left_file), str(right_file),
             "--out", str(directory / "comparison"), "--temp-dir", str(directory)],
            text=True, capture_output=True,
        )

    def test_coordinates_ignore_ids_and_preserve_ambiguous_crosswalks(self):
        with tempfile.TemporaryDirectory() as work:
            directory = Path(work)
            result = self.run_comparison(
                "sample#1#chr1\t30\t32\t3,0,+\n"
                "GRCh38#0#chr1\t10\t12\t1,0,+\n"
                "sample#1#chr1\t10\t12\t1,0,+\n"
                "sample#1#chr1\t10\t12\t1,0,+\n"
                "sample#1#chr1\t20\t22\t>1 >2\n"
                "sample#2#chr1\t10\t12\t1,0,+\n"
                "sample#1#chr1\t20\t22\t9,0,+\n",
                "sample#1#chr1\t20\t22\t100,0,+\n"
                "sample#1#chr1\t10\t12\t101,0,+\n"
                "sample#1#chr1\t40\t42\t102,0,+\n"
                "CHM13#0#chr1\t10\t12\t101,0,+\n"
                "sample#2#chr1\t11\t13\t103,0,+\n",
                directory,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            output = directory / "comparison"
            self.assertEqual((output / "shared.tsv").read_text(),
                             "path\tstart\tend\tleft_cpg\tright_cpg\n"
                             "sample#1#chr1\t10\t12\t1,0,+\t101,0,+\n"
                             "sample#1#chr1\t20\t22\t9,0,+\t100,0,+\n"
                             "sample#1#chr1\t20\t22\t>1 >2\t100,0,+\n")
            self.assertEqual((output / "left_only.bed").read_text(),
                             "sample#1#chr1\t30\t32\t3,0,+\n"
                             "sample#2#chr1\t10\t12\t1,0,+\n")
            self.assertEqual((output / "right_only.bed").read_text(),
                             "sample#1#chr1\t40\t42\t102,0,+\n"
                             "sample#2#chr1\t11\t13\t103,0,+\n")
            self.assertEqual((output / "left_uncompared.bed").read_text(),
                             "GRCh38#0#chr1\t10\t12\t1,0,+\n")
            self.assertEqual((output / "right_uncompared.bed").read_text(),
                             "CHM13#0#chr1\t10\t12\t101,0,+\n")

    def test_invalid_inputs_do_not_create_output(self):
        for left, right, message in [
            ("p\t1\t3\ta", "p\t1\t3\tb\n", "incomplete row"),
            ("p\t1\t4\ta\n", "p\t1\t3\tb\n", "two-base"),
            ("p\tbad\t3\ta\n", "p\t1\t3\tb\n", "invalid BED"),
            ("p\t1\t3\n", "p\t1\t3\tb\n", "four nonempty"),
            ("p\t1\t3\ta\n", "q\t1\t3\tb\n", "No common paths"),
        ]:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as work:
                directory = Path(work)
                result = self.run_comparison(left, right, directory)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertFalse((directory / "comparison").exists())

    def test_existing_outputs_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as work:
            directory = Path(work)
            output = directory / "comparison"
            output.mkdir()
            sentinel = output / "shared.tsv"
            sentinel.write_text("keep\n")
            result = self.run_comparison("p\t1\t3\ta\n", "p\t1\t3\tb\n", directory)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("already exists", result.stderr)
            self.assertEqual(sentinel.read_text(), "keep\n")

    def test_truncated_gzip_is_rejected(self):
        with tempfile.TemporaryDirectory() as work:
            directory = Path(work)
            left = directory / "left.bed.gz"
            right = directory / "right.bed"
            left.write_bytes(gzip.compress(b"p\t1\t3\ta\n")[:-4])
            right.write_text("p\t1\t3\tb\n")
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(left), str(right),
                 "--out", str(directory / "comparison")],
                text=True, capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Error:", result.stderr)
            self.assertFalse((directory / "comparison").exists())


if __name__ == "__main__":
    unittest.main()
