#!/usr/bin/env python3
import argparse
import gzip
import io
from pathlib import Path
import sqlite3
import sys
import tempfile


def load_tracks(db, table, filename):
    with open(filename, "rb") as raw:
        stream = gzip.GzipFile(fileobj=raw) if raw.peek(2).startswith(b"\x1f\x8b") else raw
        with io.TextIOWrapper(stream) as tracks:
            batch = []
            for number, line in enumerate(tracks, 1):
                if not line.endswith("\n"):
                    raise ValueError(f"{filename}:{number}: incomplete row (missing newline)")
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 4 or not fields[0] or not fields[3]:
                    raise ValueError(f"{filename}:{number}: expected four nonempty BED fields")
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as error:
                    raise ValueError(f"{filename}:{number}: invalid BED coordinates") from error
                if start < 0 or end - start != 2:
                    raise ValueError(f"{filename}:{number}: expected a two-base CpG interval")
                batch.append((fields[0], start, end, fields[3]))
                if len(batch) == 10000:
                    db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?, ?, ?, ?)", batch)
                    batch.clear()
            db.executemany(f"INSERT OR IGNORE INTO {table} VALUES (?, ?, ?, ?)", batch)
    db.commit()


def write_rows(filename, rows, header=None):
    count = 0
    with open(filename, "w") as out:
        if header:
            out.write(header + "\n")
        for row in rows:
            out.write("\t".join(str(value) for value in row) + "\n")
            count += 1
    return count


def compare(left, right, output, temp_dir=None):
    if output.exists():
        raise ValueError(f"Output directory already exists: {output}")
    with tempfile.TemporaryDirectory(dir=temp_dir, prefix="compare-cpg-") as work:
        with sqlite3.connect(str(Path(work) / "tracks.sqlite")) as db:
            for table in ["left_track", "right_track"]:
                db.execute(f"""
                    CREATE TABLE {table} (
                        path TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL,
                        motif TEXT NOT NULL, PRIMARY KEY (path, start, end, motif)
                    ) WITHOUT ROWID
                """)
            load_tracks(db, "left_track", left)
            load_tracks(db, "right_track", right)
            db.execute("""
                CREATE TABLE common_paths AS
                SELECT path FROM left_track INTERSECT SELECT path FROM right_track
            """)
            db.execute("CREATE UNIQUE INDEX common_path_index ON common_paths(path)")
            common_count = db.execute("SELECT COUNT(*) FROM common_paths").fetchone()[0]
            if common_count == 0:
                raise ValueError("No common paths with CpG records; check path naming and inputs")

            output.mkdir(parents=True)
            shared = write_rows(
                output / "shared.tsv",
                db.execute("""
                    SELECT l.path, l.start, l.end, l.motif, r.motif
                    FROM left_track l JOIN right_track r
                    ON l.path = r.path AND l.start = r.start AND l.end = r.end
                    ORDER BY l.path, l.start, l.end, l.motif, r.motif
                """),
                "path\tstart\tend\tleft_cpg\tright_cpg",
            )
            counts = {"shared coordinate/ID pairs": shared}
            for side, other in [("left", "right"), ("right", "left")]:
                counts[f"{side} only occurrences"] = write_rows(
                    output / f"{side}_only.bed",
                    db.execute(f"""
                        SELECT t.path, t.start, t.end, t.motif FROM {side}_track t
                        JOIN common_paths p ON t.path = p.path
                        WHERE NOT EXISTS (
                            SELECT 1 FROM {other}_track o
                            WHERE o.path = t.path AND o.start = t.start AND o.end = t.end
                        )
                        ORDER BY t.path, t.start, t.end, t.motif
                    """),
                )
                counts[f"{side} uncompared occurrences"] = write_rows(
                    output / f"{side}_uncompared.bed",
                    db.execute(f"""
                        SELECT t.path, t.start, t.end, t.motif FROM {side}_track t
                        WHERE NOT EXISTS (SELECT 1 FROM common_paths p WHERE p.path = t.path)
                        ORDER BY t.path, t.start, t.end, t.motif
                    """),
                )
            counts["common paths with CpG records"] = common_count
            for label, count in counts.items():
                print(f"{label}: {count}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Compare CpG occurrences in two panmethyl --lift BED tracks by path coordinates.",
    )
    parser.add_argument("left", type=Path, help="First GFA_paths.bed[.gz]")
    parser.add_argument("right", type=Path, help="Second GFA_paths.bed[.gz]")
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--temp-dir", type=Path, help="Directory for the disk-backed SQLite database")
    args = parser.parse_args()
    try:
        compare(args.left, args.right, args.out, args.temp_dir)
    except (ValueError, OSError, EOFError, sqlite3.Error) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
