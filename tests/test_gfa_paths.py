import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class GfaPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(
            ["cargo", "build", "--quiet", "--manifest-path",
             str(ROOT / "tagtobed/Cargo.toml"),
             "--bin", "lift_offsets", "--bin", "lift_edges"],
            check=True,
        )

    def tracks(self, records):
        with tempfile.TemporaryDirectory() as directory:
            graph = Path(directory) / "graph.gfa"
            graph.write_text(
                "S\t1\tACGC\nS\t2\tGCGT\nL\t1\t+\t2\t+\t0M\n"
                + records
            )
            index = subprocess.run(
                [sys.executable,
                 str(ROOT / "module/resources/usr/bin/index_nucleotide.py"),
                 str(graph), "CG"],
                check=True, text=True, capture_output=True,
            ).stdout
            index_file = Path(directory) / "index.tsv.gz"
            with gzip.open(index_file, "wt") as out:
                out.write(index)
            tracks = []
            for binary in ["lift_offsets", "lift_edges"]:
                result = subprocess.run(
                    [str(ROOT / "tagtobed/target/debug" / binary),
                     str(index_file), str(graph)],
                    check=True, text=True, capture_output=True,
                )
                tracks.extend(result.stdout.splitlines())
            return sorted(tracks)

    def test_mixed_paths_and_walks_with_sequence_offsets(self):
        tracks = self.tracks(
            "P\tlegacy[100-108]\t1+,2+\t*\n"
            "W\tsample\t1\tchr1\t100\t108\t>1>2\tXX:Z:optional\n"
            "W\tsample\t2\tchr1\t200\t208\t>1>2\n"
        )
        expected = []
        for name, start in [("legacy", 100), ("sample#1#chr1", 100),
                            ("sample#2#chr1", 200)]:
            expected.extend([
                f"{name}\t{start + 1}\t{start + 3}\t1,1,+",
                f"{name}\t{start + 3}\t{start + 5}\t>1 >2",
                f"{name}\t{start + 5}\t{start + 7}\t2,1,+",
            ])
        self.assertEqual(tracks, sorted(expected))

    def test_reverse_walks_and_single_node_unknown_start(self):
        tracks = self.tracks(
            "W\tsample\t0\treverse\t10\t18\t<2<1\n"
            "W\tsample\t0\tsingle\t*\t*\t>1\n"
            "W\tsample\t0\tmixed\t0\t8\t>1<2\n"
        )
        self.assertEqual(tracks, sorted([
            "sample#0#reverse\t12\t14\t2,1,+",
            "sample#0#reverse\t13\t15\t>1 >2",
            "sample#0#reverse\t16\t18\t1,1,+",
            "sample#0#single\t1\t3\t1,1,+",
            "sample#0#mixed\t1\t3\t1,1,+",
            "sample#0#mixed\t6\t8\t2,1,+",
        ]))


if __name__ == "__main__":
    unittest.main()
