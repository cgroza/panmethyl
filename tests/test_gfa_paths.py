import gzip
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


def bgzf_block(data):
    member = gzip.compress(data, mtime=0)
    return (
        member[:3] + b"\x04" + member[4:10]
        + b"\x06\x00BC\x02\x00" + struct.pack("<H", len(member) + 7)
        + member[10:]
    )


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

    def test_compressed_graphs_match_plain_index_and_tracks(self):
        data = (
            b"S\t1\tACGC\nS\t2\tGCGT\nL\t1\t+\t2\t+\t0M\n"
            b"P\tlegacy[100-108]\t1+,2+\t*\n"
            b"W\tsample\t0\tchr1\t10\t18\t<2<1\n"
        )
        bgzf = b"".join(
            bgzf_block(data[i:i + 32]) for i in range(0, len(data), 32)
        ) + bgzf_block(b"")
        variants = {
            "plain.gfa": data,
            "graph.gfa.gz": gzip.compress(data),
            "concatenated.gfa.gz": gzip.compress(data[:40]) + gzip.compress(data[40:]),
            "graph.gfa.bgz": bgzf,
            "compressed-without-extension": gzip.compress(data),
        }
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            expected = None
            for name, content in variants.items():
                with self.subTest(format=name):
                    graph = directory / name
                    graph.write_bytes(content)
                    outputs = []
                    for motif in ["A", "CG"]:
                        outputs.append(subprocess.run(
                            [sys.executable,
                             str(ROOT / "module/resources/usr/bin/index_nucleotide.py"),
                             str(graph), motif],
                            check=True, text=True, capture_output=True,
                        ).stdout)
                    index = directory / "index.tsv.gz"
                    with gzip.open(index, "wt") as out:
                        out.write(outputs[-1])
                    for binary in ["lift_offsets", "lift_edges"]:
                        result = subprocess.run(
                            [str(ROOT / "tagtobed/target/debug" / binary),
                             str(index), str(graph)],
                            check=True, text=True, capture_output=True,
                        )
                        outputs.append(sorted(result.stdout.splitlines()))
                    if expected is None:
                        expected = outputs
                    else:
                        self.assertEqual(outputs, expected)

    def test_truncated_compressed_graph_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            graph = Path(directory) / "truncated.gfa.gz"
            graph.write_bytes(gzip.compress(b"S\t1\tACGC\n")[:-4])
            index = Path(directory) / "index.tsv.gz"
            with gzip.open(index, "wt") as out:
                out.write("1\t1\t+\t1,1,+\n")
            commands = [
                [sys.executable,
                 ROOT / "module/resources/usr/bin/index_nucleotide.py", graph, "CG"],
                *[
                    [ROOT / "tagtobed/target/debug" / binary, index, graph]
                    for binary in ["lift_offsets", "lift_edges"]
                ],
            ]
            for command in commands:
                with self.subTest(tool=str(command[0])):
                    result = subprocess.run(
                        [str(arg) for arg in command], text=True, capture_output=True,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertTrue(result.stderr)

    @unittest.skipUnless(shutil.which("nextflow"), "Nextflow is not installed")
    def test_pipeline_indexes_and_lifts_compressed_graphs(self):
        data = (
            b"S\t1\tACGC\nS\t2\tGCGT\nL\t1\t+\t2\t+\t0M\n"
            b"P\tlegacy\t1+,2+\t*\nW\tsample\t0\tchr1\t10\t18\t>1>2\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            wrapper = directory / "index_nucleotide.py"
            wrapper.write_text(
                "#!/bin/sh\nexec "
                + shlex.join([
                    sys.executable,
                    str(ROOT / "module/resources/usr/bin/index_nucleotide.py"),
                ]) + ' "$@"\n'
            )
            wrapper.chmod(0o755)
            config = directory / "test.config"
            config.write_text(
                "apptainer.enabled = false\n"
                "nextflow.enable.moduleBinaries = false\n"
                "process {\n"
                f"  beforeScript = 'export PATH=\"{directory}:"
                f"{ROOT}/tagtobed/target/debug:$PATH\"'\n"
                "  withName: index_graph { memory = '1 GB' }\n"
                "  withName: lift_nucleotides { memory = '1 GB' }\n"
                "}\n"
            )
            bgzf = bgzf_block(data[:40]) + bgzf_block(data[40:]) + bgzf_block(b"")
            expected = None
            for suffix, content in [
                (".gfa", data), (".gfa.gz", gzip.compress(data)),
                (".gfa.bgz", bgzf), (".gfa.bgzf", bgzf),
            ]:
                with self.subTest(format=suffix):
                    graph = directory / ("graph" + suffix)
                    graph.write_bytes(content)
                    output = directory / ("out" + suffix)
                    result = subprocess.run(
                        ["nextflow", "run", str(ROOT), "-c", str(config),
                         "-ansi-log", "false", "--graph", str(graph),
                         "--lift", "true", "--out", str(output)],
                        cwd=directory, text=True, capture_output=True, timeout=180,
                        env={**os.environ, "NXF_SYNTAX_PARSER": "v2",
                             "NXF_OFFLINE": "true"},
                    )
                    self.assertEqual(
                        result.returncode, 0, result.stdout + result.stderr,
                    )
                    with gzip.open(output / "index/index.csv.gz", "rt") as index:
                        index_data = index.read()
                    with gzip.open(output / "tracks/GFA_paths.bed.gz", "rt") as tracks:
                        track_data = sorted(tracks.read().splitlines())
                    actual = (
                        (output / "index/node_sizes.csv").read_text(),
                        (output / "index/nodes_list.csv").read_text(),
                        index_data, track_data,
                    )
                    if expected is None:
                        expected = actual
                    else:
                        self.assertEqual(actual, expected)

            failed_output = directory / "failed-index"
            result = subprocess.run(
                ["nextflow", "run", str(ROOT), "-c", str(config),
                 "-ansi-log", "false", "--graph", str(graph),
                 "--motif", "CGT", "--out", str(failed_output)],
                cwd=directory, text=True, capture_output=True, timeout=180,
                env={**os.environ, "NXF_SYNTAX_PARSER": "v2",
                     "NXF_OFFLINE": "true"},
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unsupported motif", result.stdout + result.stderr)
            self.assertFalse((failed_output / "index/index.csv.gz").exists())


if __name__ == "__main__":
    unittest.main()
