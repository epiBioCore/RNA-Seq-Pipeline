#!/usr/bin/env nextflow

/*
* parameters
* raw_fastq path
* threads
* length
* annotation

*/

Channel
	.fromFilePairs("Raw_fastq/*{1,2}.fq.gz")
	.ifEmpty { exit 1, "Cannot find any fastq files"}
	.into { read_files_star1; read_files_star2; read_files_fastqc }

Channel
	.fromPath("genes.gtf")
	.ifEmpty { exit 1, "Cannot find gtf file" }
	.into { star1_gtf; star2_gtf}

Channel
	.fromPath("genome.fa")
	.ifEmpty{ exit 1, "Cannot find fasta file" }
	.into {star1_fasta; star2_fasta }



/*
*FastQC
*/

process fastqc {
        publishDir "FastQC", mode: "copy"
                //saveAs:{filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

        input:
                set val(name),file(fastq) from read_files_fastqc
        output:
                file  "*_fastqc.{html,zip}" into fastqc_results
        script:
                "fastqc $fastq"
}



/*
* Building STAR index for 1pass
*/

process buildSTARIndex1pass {

	input:	
		file fasta from star1_fasta	
		file gtf from star1_gtf

	output:
		file "STAR1passIndex" into star_index
	script:
		"""
                mkdir STAR1passIndex
   
		STAR --runThreadN 14 --runMode genomeGenerate \\
		--genomeDir STAR1passIndex \\
		--genomeFastaFiles $fasta \\
		--sjdbGTFfile $gtf \\
		--sjdbOverhang 149
		"""
}



/*
* Aligning with STAR for 1 pass
* I will want to buffer the input for the read_files_star1
*/

process TestSTAR1pass {

	input:
		set val(x),file(reads) from read_files_star1
		file index from star_index.collect()

	output:
		file "*SJ.out.tab" into star1pass_gtf
 		file "*.bam" into star1pass_aln

	script:
		"""
		
		STAR --runThreadN 14 \\
		--genomeDir $index \\
		--readFilesIn $reads \\
		--readFilesCommand gunzip -c \\
		--runMode alignReads \\
		--outSAMattributes All \\
		--alignSJoverhangMin 8 \\
		--alignSJDBoverhangMin 1 \\
		--outFilterMismatchNmax 999 \\
		--outFilterMismatchNoverLmax 0.04 \\
		--alignIntronMin 20 \\
		--alignIntronMax 1000000 \\
		--alignMatesGapMax 1000000 \\
		--outSAMstrandField intronMotif \\
		--outSAMtype BAM SortedByCoordinate \\
		--outFileNamePrefix $x
		"""
			
}

/*
* building STAR index for second pass
*/

process STAR2passIndex {

	publishDir "STAR2passIndex", mode : "copy"

	input:
		file SJ_gtfs from star1pass_gtf.collect()
		file fasta from star2_fasta
		file gtf from star2_gtf

	output:
		file "STAR2passIndex" into Star2pass_index


	script:

	"""
	mkdir STAR2passIndex

	 STAR --runThreadN 14 \\
                --runMode genomeGenerate \\
                --genomeDir STAR2passIndex \\
                --genomeFastaFiles $fasta \\
		--sjdbFileChrStartEnd $SJ_gtfs \\
                --sjdbGTFfile $gtf \\
                --sjdbOverhang 149
      """

}

/*
* align bam files again with updated genome index
*/

process STAR2pass {
	publishDir "STAR2pass", mode : "copy"


	input:
		set val(x), file(reads) from read_files_star2
		file index from Star2pass_index.collect()

	output:
		file "*.bam" into star2pass_aln
		file "*.final.out" into star2pass_logs

	
	script:
                """
                STAR --runThreadN 14 \\
                --genomeDir $index \\
                --readFilesIn $reads \\
                --readFilesCommand gunzip -c \\
                --runMode alignReads \\
                --outSAMattributes All \\
                --alignSJoverhangMin 8 \\
                --alignSJDBoverhangMin 1 \\
                --outFilterMismatchNmax 999 \\
                --outFilterMismatchNoverLmax 0.04 \\
                --alignIntronMin 20 \\
                --alignIntronMax 1000000 \\
                --alignMatesGapMax 1000000 \\
                --outSAMstrandField intronMotif \\
                --outSAMtype BAM SortedByCoordinate \\
                --outFileNamePrefix $x
		"""

}


