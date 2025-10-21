#!/usr/bin/python3
import polars
import sys
import gzip

schema = {'node':polars.String, 'pos':polars.Int64, 'strand':polars.String, 'gaf_depth':polars.Float32, 'gaf_score':polars.Float32}
columns = ['node', 'pos', 'strand', 'depth', 'score']

chunk = polars.read_csv(sys.argv[2], separator=' ', has_header=False,
                               new_columns=columns,
                               schema = schema)

if len(sys.argv) > 3:
    for f in sys.argv[3:]:
        chunk = polars.concat([chunk, polars.read_csv(f, separator=' ', has_header=False,
                                                    new_columns=columns,
                                                    schema=schema)],
                              how = 'vertical') \
                    .group_by('node', 'pos', 'strand') \
                    .agg(score=polars.col("depth").dot("gaf_score")/polars.col("depth").sum(),
                         depth=polars.col("depth").sum()) \
                    .fill_nan(0)

with gzip.open(sys.argv[1], 'wb') as csv:
    chunk.sort('node', 'pos', 'strand').write_csv(csv, separator=' ')
