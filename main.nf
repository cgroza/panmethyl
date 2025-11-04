params.out = "out"
params.bams = false
params.gafs = false
params.vcfs = false
params.graph_mods = false
params.graph = "graph.gfa"
params.tag = "C+m."
params.motif = "CG"
params.aligner = "GraphAligner"
params.cpus = 40
params.memory =  '180G'
params.time = '24h'

include { * } from './module'

workflow {
  Channel.fromPath(params.graph).set{graph_ch}
  index_graph(graph_ch)

  bam_methylation_ch = channel.empty()
  graph_methylation_ch = channel.empty()

  if(params.bams) {
    Channel.fromPath(params.bams).splitCsv(header : true).map{
      row -> [row.sample, file(row.path)]
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
      row -> [row.sample, file(row.bam), file(row.gaf)]}.set{gafs_ch}
  }

  bamtags_to_methylation(gafs_ch.combine(index_graph.out.graph_index)).set{bam_methylation_ch}

  if(params.graph_mods) {
    Channel.fromPath(params.graph_mods).splitCsv(header : true).map{
      row -> [row.sample, file(row.path)]
    }.set{graph_methylation_ch}
  }

  methylation_to_csv(bam_methylation_ch.concat(graph_methylation_ch).combine(index_graph.out.graph_index)).set{csv_ch}
  merge_csv(csv_ch.groupTuple(by: 0)).set{merged_ch}

  if (params.vcfs) {
    Channel.fromPath(params.vcfs).splitCsv(header : true).map{
      row -> [row.sample, file(row.path)]}.set(vcf_ch)
    annotate_vcf(vcf_ch.combine(merged_ch, by: 0))
  }
}
