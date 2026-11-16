include { decompress_graph; annotate_VCF; index_graph; align_graphaligner; align_minigraph; align_giraffe; bamtags_to_BED; merge_CSV; BED_to_graph; annotate_BED; merge_BED; lift_epigenome; lift_nucleotides } from './module'

workflow {
  channel.fromPath(params.graph).set{graph_input_ch}
  if (params.graph.toString().toLowerCase() ==~ /.*[.](gz|bgz|bgzf)/) {
    graph_ch = decompress_graph(graph_input_ch)
  } else {
    graph_ch = graph_input_ch
  }

  graph_index = channel.empty()
  if(params.index) {
    channel.fromPath(params.index).map{index_dir ->
      [file(index_dir / 'node_sizes.csv', checkIfExists: true),
       file(index_dir / 'nodes_list.csv', checkIfExists: true),
       file(index_dir / 'index.csv.gz', checkIfExists: true)]
    }.set{graph_index}
  } else {
    index_graph(graph_ch, channel.value(params.motif)).set{graph_index}
  }

  graph_methylation_ch = channel.empty()

  gafs_ch = channel.empty()
  bams_ch = channel.empty()
  bamtags_bed_ch = channel.empty()

  if(params.lift) {
    lift_nucleotides(graph_ch, graph_index.map{index -> index[2]})
  }

  if(params.bams) {
    channel.fromPath(params.bams).splitCsv(header : true).map{
      row -> [row.sample, file(row.path, checkIfExists: true), row.preset]
    }.set{bams_ch}

    if(params.aligner == "minigraph") {
      align_minigraph(bams_ch.combine(graph_ch)).set{gafs_ch}
    }
    else if(params.aligner == "GraphAligner") {
      align_graphaligner(bams_ch.combine(graph_ch)).set{gafs_ch}
    }
    else if(params.align_graphaligner == "giraffe") {
      align_giraffe(bams_ch.combine(channel.fromPath(params.giraffe_index)).combine(channel.value(params.giraffe_index_prefix)))
    }
  }

  if (params.gafs) {
    channel.fromPath(params.gafs).splitCsv(header : true)
      .filter{row -> !row.sample.startsWith("#")}
      .map{row -> [row.sample, file(row.bam, checkIfExists: true), file(row.gaf, checkIfExists: true)]}.set{gafs_ch}
  }

  if(params.mods) {
    channel.fromPath(params.mods).splitCsv(header : true)
      .filter{row -> !row.sample.startsWith("#")}
      .map{row -> [row.sample, file(row.mods, checkIfExists: true)]}.set{bamtags_bed_ch}

  } else {
    bamtags_to_BED(gafs_ch.map{gaf -> [gaf[0], gaf[1]]}, channel.value(params.code)).set{bamtags_bed_ch}
  }


  if(params.graph_mods) {
    channel.fromPath(params.graph_mods).splitCsv(header : true).map{
      row -> [row.sample, file(row.path, checkIfExists: true)]
    }.set{graph_methylation_ch}
  } else {
    lift_epigenome(bamtags_bed_ch.combine(gafs_ch.map{gaf -> [gaf[0], gaf[2]]}, by: 0).combine(graph_index)).set{graph_methylation_ch}
  }

  merge_CSV(graph_methylation_ch.groupTuple(by: 0).combine(graph_index)).set{merged_ch}

  if (params.vcfs) {
    channel.fromPath(params.vcfs).splitCsv(header : true).map{
      row -> [row.sample, file(row.path, checkIfExists: true)]}.set{vcf_ch}
    annotate_VCF(vcf_ch.combine(merged_ch, by: 0))
  }

  if (params.bed) {
    BED_to_graph(graph_ch.combine(channel.fromPath(params.bed))).set{bed_ch}
    annotate_BED(merged_ch.combine(bed_ch).combine(graph_index.map{index -> index[0]}))
  }
}
