params.out = "out"
params.bams = "bams.csv"
params.graph = "graph.gfa"
params.tag = "C+m."
params.aligner = "GraphAligner"

process index_graph
{
  cpus 1
  memory "20G"
  time '3h'

  publishDir "${params.out}/index/", mode: 'copy'

  input:
  path(graph_path)

  output:
  tuple path("node_sizes.csv"), path("nodes_list.csv"), path("cpg_index.csv.gz"), emit: graph_index

  script:
  """
  awk '\$1 == "S" {print \$2, length(\$3)}' ${graph_path} > node_sizes.csv
  awk '{print \$1}' node_sizes.csv > nodes_list.csv
  index_cpg.py ${graph_path} | gzip > cpg_index.csv.gz
  """
}
process align_graphaligner {
  cpus 40
  memory '180G'
  time '12h'

  publishDir "${params.out}/gafs/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(graph_path)

  output:
  tuple val(sample_name), path("${sample_name}.gaf.gz"), emit: gafs

  script:
  """
  samtools fasta --threads 40 ${bam_path} | pigz  > ${sample_name}.fa.gz
  GraphAligner -x vg -a ${sample_name}.gaf -g ${graph_path} -f ${sample_name}.fa.gz
  pigz ${sample_name}.gaf
  """

}
process align_minigraph {
  cpus 40
  memory '180G'
  time '12h'

  publishDir "${params.out}/gafs/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(graph_path)

  output:
  tuple val(sample_name), path("${sample_name}.gaf.gz"), emit: gafs

  script:
  """
  samtools fasta --threads 40 ${bam_path} | pigz  > ${sample_name}.fa.gz
  minigraph --vc -c -N 1 -t 40 ${graph_path} ${sample_name}.fa.gz | pigz > ${sample_name}.gaf.gz
  """
}

process bamtags_to_methylation {
  cpus 2
  time '6h'
  memory '20G'

  publishDir "${params.out}/methylation/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(gaf_path),
    path(node_sizes), path(nodes_list), path(cpg_index)

  output:
  path("${sample_name}.graphMethylation")
  path("${sample_name}.graph5mC")

  script:
  """
  samtools index ${bam_path}
  tagtobed -b ${bam_path} -B ${params.tag} | pigz > ${sample_name}.5mC.gz

  join -t \$'\\t' -1 1 -2 1 <(gunzip -c ${gaf_path} | sort ) \
    <(gunzip -c ${sample_name}.5mC.gz | sort ) | \
    lift_5mC.py ${node_sizes} ${sample_name}.graph5mC
  nodes_methylation.py ${nodes_list} ${sample_name}.graph5mC ${cpg_index} | sort -t' ' -k1,1 -k2,2n | pigz > ${sample_name}.graphMethylation
  """
}

workflow {
  Channel.fromPath(params.graph).set{graph_ch}
  Channel.fromPath(params.bams).splitCsv(header : true).map{
    row -> [row.sample, file(row.path)]
  }.set{bams_ch}

  index_graph(graph_ch)
  if(params.aligner == "minigraph") {
    align_minigraph(bams_ch.combine(graph_ch)).set{gafs_ch}
  }
  else if(params.aligner == "GraphAligner") {
    align_graphaligner(bams_ch.combine(graph_ch)).set{gafs_ch}
  }
  bamtags_to_methylation(bams_ch.combine(gafs_ch, by : 0).combine(index_graph.out.graph_index))
}
