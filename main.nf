#! /usr/bin/env nextflow

def helpMessage() {
  log.info """
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run main.nf --inputdir "inputDirName" --outdir "outDirName"

  Arguments:
  --inputdir                The direcoty containing the input fastq or fastq.gz files (default is the current direcoty)
  --outdir                  The directory to place processing results (default is ./outdir)
  --paired_end              the fastq filename pattern (default *_R{1,2}*.fastq.gz)
  --single_end              the fastq filename pattern (default *.fastq.gz)
  --adapter                 The adapter fasta file (default is ./ref/adapters.fa . it is from bbmap program sitting in the "resource directory")
  --arti                    The contamination database (default ./ref/artifacts.fa.gz . it is from bbmap program sitting in the "resource directory")
  --type                    Specifies that the input files are paired reads or single-end|paired-end(default is paired-end)
  --minlen                  Min length cutoff of the quality controlled reads
  --kraken2_db              Kraken program database (default is at $HOME/software/kraken2/db/minikraken_8GB_20200312)
  --help                    This usage statement.
  """.stripIndent()
}
//--host                     The host reference fasta file (defualt is ./ref/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz)
seq_status = "raw"

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

Channel.fromPath(params.inputdir, type: 'dir')
  .ifEmpty { exit 1, "Cannot find  directory: ${params.inputdir}"}
  .set{rawreads_ch}


if (params.type == 'single-end') {
    fastq_input_to_fastqc = Channel
                        .fromPath(params.inputdir + '/' + params.single_end)
                        .map { file -> tuple(file.simpleName, file) }
    fastq_input_to_qc = Channel
                        .fromPath(params.inputdir + '/' + params.single_end)
                        .map { file -> tuple(file.simpleName, file) } 
    fastq_input_to_rawstats = Channel
                        .fromPath(params.inputdir + '/' + params.single_end)
                        .map { file -> tuple(file.simpleName, file) }                   
} else {
  Channel
      .fromFilePairs(params.inputdir + '/' + params.paired_end, size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.paired_end}"}
      //.into { fastq_reads_qc; fastq_reads_trim}
      .into { fastq_input_to_rawstats; fastq_input_to_fastqc; fastq_input_to_qc}
}

/*
Step 1
*/
process stats_rawreads {
  
  publishDir "${params.outdir}/raw", mode: "copy" 

  input:
  file fastqFile from fastq_input_to_rawstats.collect()
  file dir from rawreads_ch

  output:
  file "raw_stats.txt" into stats_rawreads_results
 

  script:
  """
  readstats.pl $dir > raw_stats.txt
  """
  
}


/*
Step 1
*/
process fastqc_rawreads {
  
  tag "$name"

  publishDir "${params.outdir}/raw/fastqc", mode: 'copy',
    saveAs: {filename ->
          if (filename.indexOf("zip") > 0)     "zips/$filename"
          else if (filename.indexOf("html") > 0)    "fastqc/$filename"
          else if (filename.indexOf("txt") > 0)     "fastqc_stats/$filename"
          else null            
      }
  input:
  set val(name), file(fastqFile) from fastq_input_to_fastqc

  output:
  file ("*fastqc.{zip,html}") into raw_fastqc_results

  script:
  
  """
  fastqc -o . $fastqFile
  """
  
}


/*
QC was implemented using the procedure described in the link
https://www.protocols.io/view/illumina-fastq-filtering-gydbxs6?step=4
*/
process bbduk {
  tag "$name"
  label "small_mem"
  
  publishDir "${params.outdir}/qc/seqs", mode: 'copy' , pattern: "*.qc.fastq.gz"
  publishDir "${params.outdir}/qc/bbduk_stats", mode: 'copy' , pattern: "*stats.txt"

  input:
  set val(name), file(reads) from fastq_input_to_qc
  path adapter from params.adapter
  path arti from params.arti
  //path host from params.host

  output:
  set val(name), file("*qc.fastq.gz") into qcreads_to_fastqc
  set val(name), file("*qc.fastq.gz") into qcreads_to_kraken2
  set val(name), file("*qc.fastq.gz") into qcreads_to_stats
  file "*stats.txt" into qc_stats
  //val "test" into test
  script:
  
  prefix_pe = name
  prefix_se = name
  //println "***************= $prefix_pe"

  ram = "-Xmx${task.memory.toGiga()}g"
  println "********************ram=$ram"
  
  if (params.type == 'single-end') {
    """
      echo ${prefix_se}
      
      bbduk.sh $ram \
        in=${reads} \
        ref=${adapter} \
        out=${prefix_se}_trim.fastq.gz \
        stats=${prefix_se}.adapter_trimstats.txt \
        refstats=${prefix_se}.adapter_refstats.txt \
        threads=${task.cpus} \
        ktrim=r k=23 qtrim=rl trimq=10  mink=11 hdist=1 maxns=0 \
        ftm=5 minlen=${params.minlen} tbo tpe ow=t rcomp=f hdist2=1 ftm=5 zl=4 

    # quality trimming
    bbduk.sh $ram \
        in=${prefix_se}_trim.fastq.gz \
        ref=${arti} \
        out=${prefix_se}_qFiltered.fastq.gz \
        stats=${prefix_se}.qtrimstats.txt \
        maq=5 trimq=10 qtrim=f ordered ow=t maxns=1 \
        minlen=${params.minlen} k=25 hdist=1 zl=6 \
        threads=${task.cpus}

    #artifact filtering
    bbduk.sh $ram \
        in=${prefix_se}_trim.fastq.gz \
        ref=${arti} \
        out=${prefix_se}_artiFiltered.fastq.gz \
        outm=${prefix_se}_arctiMatch.fastq.gz \
        stats=${prefix_se}.arti_trimstats.txt \
        refstats=${prefix_se}.arti_refstats.txt \
        ordered ow=t k=20 hdist=1 zl=6 \
        minlen=${params.minlen} \
        threads=${task.cpus}

    ln -s ${prefix_se}_artiFiltered.fastq.gz ${prefix_se}.qc.fastq.gz
    """
  } 
  
  else {
    """
    echo ${prefix_pe}
    
    bbduk.sh $ram \
      in1=${reads[0]} \
      in2=${reads[1]} \
      ref=${adapter} \
      out=${prefix_pe}_trim_R1.fastq.gz \
      out2=${prefix_pe}_trim_R2.fastq.gz \
      stats=${prefix_pe}.adapter_trimstats.txt \
      refstats=${prefix_pe}.adapter_refstats.txt \
      threads=${task.cpus} \
      ktrim=r k=23 qtrim=rl trimq=10  mink=11 hdist=1 maxns=0 \
      ftm=5 minlen=${params.minlen} tbo tpe ow=t rcomp=f hdist2=1 ftm=5 zl=4 

  # quality trimming
  bbduk.sh $ram \
      in=${prefix_pe}_trim_R1.fastq.gz \
      in2=${prefix_pe}_trim_R2.fastq.gz \
      ref=${arti} \
      out=${prefix_pe}_qFiltered_R1.fastq.gz \
      out2=${prefix_pe}_qFiltered_R2.fastq.gz \
      stats=${prefix_pe}.qtrimstats.txt \
      maq=5 trimq=10 qtrim=f ordered ow=t maxns=1 \
      minlen=${params.minlen} k=25 hdist=1 zl=6 \
      threads=${task.cpus}

  #artifact filtering
  bbduk.sh $ram \
      in=${prefix_pe}_trim_R1.fastq.gz \
      in2=${prefix_pe}_trim_R2.fastq.gz \
      ref=${arti} \
      out=${prefix_pe}_artiFiltered_R1.fastq.gz \
      out2=${prefix_pe}_artiFiltered_R2.fastq.gz \
      outm=${prefix_pe}_arctiMatch.fastq.gz \
      stats=${prefix_pe}.arti_trimstats.txt \
      refstats=${prefix_pe}.arti_refstats.txt \
      ordered ow=t k=20 hdist=1 zl=6 \
      minlen=${params.minlen} \
      threads=${task.cpus}

  ln -s ${prefix_pe}_artiFiltered_R1.fastq.gz  ${prefix_pe}_R1.qc.fastq.gz
  ln -s ${prefix_pe}_artiFiltered_R2.fastq.gz  ${prefix_pe}_R2.qc.fastq.gz
  """

}
}

Channel.fromPath("${params.outdir}", type: 'dir')
  .ifEmpty { exit 1, "Cannot find  directory: ${params.outdir}"}

      //.into { fastq_reads_qc; fastq_reads_trim}
  //.view()
  .set{qcreads_ch}

process stats_qcreads {
  
  publishDir "${params.outdir}/qc", mode: "copy" 

  
  input:
  file qc_stats from qcreads_to_stats.collect()
  file dir from qcreads_ch
  
  output:
  file "qc_stats.txt" into stats_qcreads_results

  script:
  
  """
  readstats.pl $dir/qc/seqs > qc_stats.txt
  """
}

process kraken2 {
  
  label "small_mem"
  tag "$name"

  publishDir "${params.outdir}/classifier", mode: 'copy' , pattern: "*.kraken2*"
  publishDir "${params.outdir}/classifier", mode: 'copy' , pattern: "*.bracken*"

  input:
  set val(name), file(fastqFile) from qcreads_to_kraken2

  output:
  file "*kraken2*" into kraken2_results
  file "*bracken*" into bracken_results

  script:
  sname = name
  

  """
  kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --use-names --report ${sname}.kraken2.report.txt --output ${sname}.kraken2.output.txt --gzip-compressed $fastqFile
  bracken -d ${params.kraken2_db} -i ${sname}.kraken2.report.txt -o ${sname}.bracken.output.txt -w ${sname}.braken.outreport.txt 
  """

}




/*
* quality controlled fastqc
*/
process fastqc_qcreads {
  tag "$name"

  publishDir "${params.outdir}/qc/fastqc", mode: 'copy',
    saveAs: {filename ->
          if (filename.indexOf("zip") > 0)     "zips/$filename"
          else if (filename.indexOf("html") > 0)    "fastqc/$filename"
          else if (filename.indexOf("txt") > 0)     "fastqc_stats/$filename"
          else null            
      }
  input:
   set val(name), file(qc_fastqFile) from qcreads_to_fastqc

  output:
  file ("*fastqc.{zip,html}") into qc_fastqc_results

  script:
  
  """
  fastqc -o . $qc_fastqFile
  """
  
}

process multiqc {

    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    input:
    file raw_fastqc from raw_fastqc_results.collect()
    file qc_fastqc from qc_fastqc_results.collect()
    file qc_stat from qc_stats.collect()
    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

  script:
  """
  multiqc .
  """
}



workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
