#!/usr/bin/env nextflow



rawFastq = Channel.fromPath("Raw_fastq/*fq.gz")
/*
* FastQC
*/
process fastqc {
        publishDir "FastQC", mode: "copy"
                //saveAs:{filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

        input:
                file fastq from rawFastq
        output:
                file  "*_fastqc.{html,zip}" into fastqc_results
        script:
                "fastqc $fastq"
}


