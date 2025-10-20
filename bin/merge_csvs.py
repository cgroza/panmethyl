#!/usr/bin/python3
import polars
import sys

polars.concat([polars.read_csv(f, separator=' ', has_header=False,
                               new_columns=['node', 'pos', 'strand', 'gaf_depth', 'gaf_score'],
                               schema = {'node':polars.String, 'pos':polars.Int64, 'strand':polars.String, 'gaf_depth':polars.Float32, 'gaf_score':polars.Float32}) for f in sys.argv[2:]], how = 'vertical') \
      .group_by('node', 'pos', 'strand') \
      .agg(score=polars.col("gaf_depth").dot("gaf_score") / polars.col("gaf_depth").sum(),
           depth=polars.col("gaf_depth").sum()) \
      .fill_nan(0) \
      .sort('node', 'pos', 'strand') \
      .write_csv(sys.argv[1], separator=' ')

