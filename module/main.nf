process annotate_VCF {
  publishDir "${params.out}/annotation/", mode: 'copy'
  input:
  tuple val(sample), path(vcf), path(mods)

  output:
  tuple val(sample), path("${sample}.mods.vcf.gz*"), path("${sample}.mods.tsv")

  script:
  """
  annotate_vcf.py ${vcf} ${mods} ${sample}.mods.vcf > ${sample}.mods.tsv
  bgzip ${sample}.mods.vcf
  tabix ${sample}.mods.vcf.gz
  """
}

process BED_to_graph {
  input:
  tuple path(graph), path(bed)
  output:
  tuple path("annotation.gaf"), path("annotation.bed")

  script:
  """
  awk '{\$4=\$4"_"NR; print \$0}' ${bed} | sort -k4 > annotation.bed
  vg annotate -x ${graph} -b annotation.bed -F > annotation.gaf
  """
}

process annotate_BED {
  publishDir "${params.out}/annotation", mode: 'copy'
  input:
  tuple val(sample), path(mods), path(gaf), path(bed)
  output:
  tuple val(sample), path("${sample}.bed")

  script:
  """
  annotate_bed.py ${gaf} ${mods} ${sample} > ${sample}_lifted.bed
  head -n1 ${sample}_lifted.bed > ${sample}.bed
  join -1 4 -2 1 ${bed} <(awk 'NR > 1' ${sample}_lifted.bed) | awk -v OFS='\t' '{sub(/_[0-9]+/,"",\$1); print(\$0)}' >> ${sample}.bed
  """
}

process merge_BED {
  input:
  path(beds)
  output:
  path("merged_epiannoation.bed")
  script:
  """
  template=\$(ls ${beds} | head -n1)

  echo -e "QNAME\tCHROM\tSTART\tEND\tPATH\tFORMAT" > body
  awk 'NR>1' \${template} | cut -f1-6 >> body

  for b in ${beds}
  do
  awk 'NR==1' \${b} | cut -f4 > \$b.col
  awk 'NR>1' \${b} | cut -f7 >> \$b.col
  done

  paste body *.col > merged_epiannoation.bed
  """
}

process index_graph {
  publishDir "${params.out}/index/", mode: 'copy'

  input:
  path(graph_path)
  val(motif)

  output:
  tuple path("node_sizes.csv"), path("nodes_list.csv"), path("index.csv.gz"), emit: graph_index

  script:
  """
  awk '\$1 == "S" {print \$2, length(\$3)}' ${graph_path} > node_sizes.csv
  awk '{print \$1}' node_sizes.csv > nodes_list.csv
  index_nucleotide.py ${graph_path} ${motif} | gzip > index.csv.gz
  """
}

process align_graphaligner {
  publishDir "${params.out}/gafs/", mode: 'copy', pattern: '*.gaf.gz'

  input:
  tuple val(sample_name), path(bam_path), path(graph_path)

  output:
  tuple val(sample_name), path(bam_path), path("${sample_name}.gaf.gz")

  script:
  """
  samtools fasta --threads ${task.cpus} ${bam_path} | pigz  > ${sample_name}.fa.gz
  GraphAligner -t ${task.cpus} -x vg -a ${sample_name}.gaf -g ${graph_path} -f ${sample_name}.fa.gz
  cut -f1-12,17 ${sample_name}.gaf | pigz > ${sample_name}.gaf.gz
  rm ${sample_name}.gaf
  """
}

process align_minigraph {
  publishDir "${params.out}/gafs/", mode: 'copy', pattern: '*.gaf.gz'

  input:
  tuple val(sample_name), path(bam_path), path(graph_path)

  output:
  tuple val(sample_name), path(bam_path), path("${sample_name}.gaf.gz")

  script:
  """
  samtools fasta --threads ${task.cpus} ${bam_path} | pigz  > ${sample_name}.fa.gz
  minigraph --vc -c -N 1 -t ${task.cpus} ${graph_path} ${sample_name}.fa.gz | cut -f1-12,19 | pigz > ${sample_name}.gaf.gz
  """
}

process bamtags_to_BED {
  publishDir "${params.out}/lifted/", mode: 'copy'

  input:
  tuple val(sample_name), path(bam_path), path(gaf_path),
    path(node_sizes), path(nodes_list), path(index)
  val(code)
  val(missing)

  output:
  tuple val(sample_name), path("${sample_name}.graph_mods")

  script:
  """
  samtools index ${bam_path}
  tagtobed -t ${task.cpus}  -T ${code[0]} -b ${bam_path} -B ${code} -m ${missing} | pigz > ${sample_name}.mods.gz

  join -t \$'\\t' -1 1 -2 1 <(pigz -c ${gaf_path} | sort ) \
    <(pigz -c ${sample_name}.mods.gz | sort ) | \
    lift.py ${node_sizes} ${sample_name}.graph_mods
  """
}

process epigenome_to_CSV {
  publishDir "${params.out}/lifted/", mode: 'copy'

  input:
  tuple val(sample_name), path(graph_mods), path(node_sizes), path(nodes_list), path(index)

  output:
  tuple val(sample_name), path("${sample_name}.csv.gz")

  script:
  """
  nodes_levels.py ${nodes_list} ${graph_mods} ${index} | sort -t' ' -k1,1 -k2,2n | pigz > ${sample_name}.csv.gz
  """
}

process merge_CSV {
  publishDir "${params.out}/levels/", mode: 'copy'

  input:
  tuple val(sample_name), path("graph_levels*.csv.gz")

  output:
  tuple val(sample_name), path("${sample_name}.csv.gz")

  script:
  """
  merge_csvs.py ${sample_name}.csv.gz graph_levels*.csv.gz
  """
}
