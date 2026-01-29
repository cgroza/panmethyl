params.out = "out"
params.bams = false
params.gafs = false
params.vcfs = false
params.graph_mods = false
params.graph = "graph.gfa"
params.code = "C+m."
params.missing = 0
params.motif = "CG"
params.aligner = "GraphAligner"
params.cpus = 40
params.memory =  '180G'
params.time = '24h'

include { annotate_VCF; index_graph; align_graphaligner; align_minigraph; bamtags_to_bed; epigenome_to_CSV; merge_CSV, bed_to_graph, epiannotate_bed, epiannotate_bed, merge_epiannotation} from './module'

workflow {
  Channel.fromPath(params.graph).set{graph_ch}
  index_graph(graph_ch, channel.value(params.motif))

  bam_methylation_ch = channel.empty()
  graph_methylation_ch = channel.empty()

  gafs_ch = channel.empty()
  bams_ch = channel.empty()

  if(params.bams) {
    Channel.fromPath(params.bams).splitCsv(header : true).map{
      row -> [row.sample, file(row.path, checkIfExists: true)]
    }.set{bams_ch}

    if(params.aligner == "minigraph") {
      align_minigraph(bams_ch.combine(graph_ch)).set{gafs_ch}
    }
    else if(params.aligner == "GraphAligner") {
      align_graphaligner(bams_ch.combine(graph_ch)).set{gafs_ch}
    }
  }

  if (params.gafs) {
    Channel.fromPath(params.gafs).splitCsv(header : true).map{
      row -> [row.sample, file(row.bam, checkIfExists: true), file(row.gaf, checkIfExists: true)]}.set{gafs_ch}
  }

  bamtags_to_bed(gafs_ch.combine(index_graph.out.graph_index),
                         channel.value(params.code),
                         channel.value(params.missing)).set{bam_methylation_ch}

  if(params.graph_mods) {
    Channel.fromPath(params.graph_mods).splitCsv(header : true).map{
      row -> [row.sample, file(row.path, checkIfExists: true)]
    }.set{graph_methylation_ch}
  }

  epigenome_to_CSV(bam_methylation_ch.concat(graph_methylation_ch).combine(index_graph.out.graph_index)).set{csv_ch}
  merge_CSV(csv_ch.groupTuple(by: 0)).set{merged_ch}

  if (params.vcfs) {
    Channel.fromPath(params.vcfs).splitCsv(header : true).map{
      row -> [row.sample, file(row.path, checkIfExists: true)]}.set(vcf_ch)
    annotate_vcf(vcf_ch.combine(merged_ch, by: 0))
  }

  if (params.bed) {
    bed_to_graph(index_graph.out.graph_index.combine(Channel.fromPath(params.bed))).set{bed_ch}
    epiannotate_bed(merged_ch.combine(bed_ch))
  }
}
