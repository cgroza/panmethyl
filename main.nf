params.out = "out"
params.bams = false
params.graph5mc = false
params.graph = "graph.gfa"
params.tag = "C+m."
params.aligner = "GraphAligner"
params.cpus = 40
params.memory =  '180G'
params.time = '24h'

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
  cpus params.cpus
  memory params.memory
  time params.time

  publishDir "${params.out}/gafs/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(graph_path)

  output:
  tuple val(sample_name), path(bam_path), path("${sample_name}.gaf.gz"), emit: gafs

  script:
  """
  samtools fasta --threads ${params.cpus} ${bam_path} | pigz  > ${sample_name}.fa.gz
  GraphAligner -t ${params.cpus} -x vg -a ${sample_name}.gaf -g ${graph_path} -f ${sample_name}.fa.gz
  pigz ${sample_name}.gaf
  """

}
process align_minigraph {
  cpus params.cpus
  memory params.memory
  time params.time

  publishDir "${params.out}/gafs/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(graph_path)

  output:
  tuple val(sample_name), path(bam_path), path("${sample_name}.gaf.gz"), emit: gafs

  script:
  """
  samtools fasta --threads ${params.cpus} ${bam_path} | pigz  > ${sample_name}.fa.gz
  minigraph --vc -c -N 1 -t ${params.cpus} ${graph_path} ${sample_name}.fa.gz | pigz > ${sample_name}.gaf.gz
  """
}

process bamtags_to_methylation {
  cpus 2
  time '6h'
  memory '60G'

  publishDir "${params.out}/methylation/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(gaf_path),
    path(node_sizes), path(nodes_list), path(cpg_index)

  output:
  tuple val(sample_name), path("${sample_name}.graph5mC")

  script:
  """
  samtools index ${bam_path}
  tagtobed -b ${bam_path} -B ${params.tag} | pigz > ${sample_name}.5mC.gz

  join -t \$'\\t' -1 1 -2 1 <(gunzip -c ${gaf_path} | sort ) \
    <(gunzip -c ${sample_name}.5mC.gz | sort ) | \
    lift_5mC.py ${node_sizes} ${sample_name}.graph5mC
  """
}

process methylation_to_csv {
  cpus 1
  time '6h'
  memory '60G'

  publishDir "${params.out}/methylation/", mode: 'copy'

  input:
  tuple val(sample_name), path(graph_methylation), path(node_sizes), path(nodes_list), path(cpg_index)

  output:
  tuple val(sample_name), path("${sample_name}.graphMethylation")

  script:
  """
  nodes_methylation.py ${nodes_list} ${sample_name}.graph5mC ${cpg_index} | sort -t' ' -k1,1 -k2,2n | pigz > ${sample_name}.graphMethylation
  """
}

process merge_csv {
  cpus 1
  time '6h'
  memory '60G'

  input:
  tuple val(sample_name), path("graphMethylation*.csv")

  output:
  path("${sample_name}.csv.gz")

  script:
  """
  merge_csvs.py ${sample_name}.csv.gz graphMethylation*.csv
  """
}

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
    bamtags_to_methylation(gafs_ch.combine(index_graph.out.graph_index)).set{bam_methylation_ch}
  }

  if(params.graph5mc) {
    Channel.fromPath(params.graph5mc).splitCsv(header : true).map{
      row -> [row.sample, file(row.path)]
    }.set{graph_methylation_ch}
  }
  methylation_to_csv(bam_methylation_ch.concat(graph_methylation_ch).combine(index_graph.out.graph_index)).set{csv_ch}
  merge_csv(csv_ch.groupTuple(by: 0))
}
