#!/usr/bin/python3
import polars
import sys

polars.concat([polars.read_csv(f, separator=' ', has_header=False,
                               new_columns=['node', 'pos', 'strand', 'gaf_depth', 'gaf_score'],
                               schema = [polars.String, polars.Int64, polars.String, polars.Float32, polars.Float32]) for f in sys.argv[2:]], how = 'vertical') \
      .group_by('node', 'pos', 'strand') \
      .agg(score=polars.col("gaf_depth").dot("gaf_score") / polars.col("gaf_depth").sum(),
           depth=polars.col("gaf_depth").sum()) \
      .fill_nan(0) \
      .sort('node', 'pos', 'strand') \
      .write_csv(sys.argv[1], separator=' ')

