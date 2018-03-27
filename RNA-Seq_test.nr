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
	.into { read_files; read_files_fastqc }

Channel
	.fromPath("genes.gtf")
	.ifEmpty { exit 1, "Cannot find gtf file" }
	.into { star1_gtf; star2_gtf; stringtie_gtf; stringtieMerge_gtf; stringtieCompare_gtf; stringtieFPKM_gtf}

Channel
	.fromPath("genome.fa")
	.ifEmpty{ exit 1, "Cannot find fasta file" }
	.into {star1_fasta; star2_fasta }



/*
*FastQC
*/

process fastqc {
        publishDir "FastQC", mode: "copy",
                saveAs:{filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename"}

        input:
                set val(name),file(fastq) from read_files_fastqc

        output:
                file  "*_fastqc.{html,zip}" into fastqc_results
        	file ".command.sh" into fastqc_commands

	script:
		prefix = fastq[0].toString() - ~/_1.*/
                "fastqc $fastq"
}


/*
* Trimm adaptors using Trimnomatic
*/

process adapterTrimming {
	publishDir "Trimmed_fastq", mode: "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }

	input:
		set val(name),file(read_files) from read_files

	output:
		file "*at.fq.gz" into read_files_star1,read_files_star2,trimmed_reads_fastqc
		file "*unpaired*" into unpaired_fastq
		file ".command.err" into trim_summary
		file "*.log" into trim_logs
		file ".command.sh" into trimmomatic_commands


	script:
		prefix = read_files[0].toString() - ~/_1.*/
		"""
		trimmomatic-0.33.jar PE -threads 5 -phred33 -trimlog Trimmomatic_${name}.log \\
		$read_files ${name}_1_at.fq.gz ${name}_1_unpaired.fq.gz ${name}_2_at.fq.gz ${name}_2_unpaired.fq.gz \\
		ILLUMINACLIP:/usr/share/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 
		"""

}

/* 
*Fastqc on trimmed reads
*/

process trimedFastQC {
	publishDir "Trimmed_fastQC", mode: "copy",
			saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }
	input:
		file reads from trimmed_reads_fastqc

	output:
		file "*_fastqc.{html,zip}" into trimmed_fastqc_results
		file ".command.sh" into trim_fastqc_commands	
	
	script:
		prefix = reads[0].toString() - ~/_1.*/
		"""
		fastqc $reads
		"""
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
		file ".command.sh" into star1pass_index_commands
	
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

process STAR1pass {

	input:
		file reads from read_files_star1
		file index from star_index.collect()

	output:
		file "*SJ.out.tab" into star1pass_gtf
 		file "*.bam" into star1pass_aln
		file ".command.sh" into star1pass_commands
	script:
		
		prefix = reads[0].toString() - ~/_1.*/
		"""
		STAR \\
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
		--outFileNamePrefix $prefix
		"""
			
}

/*
* building STAR index for second pass
*/

process STAR2passIndex {

	publishDir "STAR2passIndex", mode : "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "star2passIndex_commands.sh" : "$filename" }

	input:
		file SJ_gtfs from star1pass_gtf.collect()
		file fasta from star2_fasta
		file gtf from star2_gtf

	output:
		file "STAR2passIndex" into Star2pass_index
		file ".command.sh" into star2pass_index_commands

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
	publishDir "STAR2pass", mode : "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }	
	input:
		file reads from read_files_star2
		file index from Star2pass_index.collect()

	output:
		file "*.bam" into star2pass_aln
		file "*.final.out" into star2pass_logs
		file ".command.sh" into star2pass_commands
	
	script:
		prefix = reads[0].toString() - ~/_1.*/
                """
                STAR \\
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
                --outFileNamePrefix $prefix
		"""

}

/*
* filtering bam files to contain only properly paired reads
*/

process filterBam {
	publishDir "STAR2pass", mode : "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }

	input:
		file bam from star2pass_aln

	output:
		file "${prefix}_filtered_sortedByCoord.out.bam" into filtered_bams_assembly,filtered_bams_quantify
		file "*.bai" into bam_indices
		file ".command.sh" into filtering_commands
	script:
		prefix = bam.toString() - ~/Aligned.sortedByCoord.out.bam/
		"""
		samtools view -b -f 0x2 $bam > ${prefix}_filtered_sortedByCoord.out.bam
		samtools index ${prefix}_filtered_sortedByCoord.out.bam
		"""
}

/*
*stringtie assembly
*/
process stringtieAssembly {
	publishDir "Assembly", mode: "copy",
 		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }
	input:
		file bam from filtered_bams_assembly
		file gtf from stringtie_gtf.collect()

	output:
		file "*gtf" into assembled_gtfs, assembled_transcripts_fn
		file ".command.sh" into stringtieAssembly_commands

	script:
		prefix = bam.toString() - ~/_filtered.*/
		"""
		stringtie -p 14 -G $gtf -o ${bam.baseName}.gtf $bam
		"""
}

/*
*create a file with the names of the gtf files
*/

assembled_transcripts_fn
	.collectFile() {file -> ['mergedlist.txt', file.name + '\n'] }
	.set {GTFfilenames}


process stringtieMerge{
	publishDir "Assembly", mode: "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "stringtie_merge_commands.sh" : "$filename" }

	input:
		file gtf_filenames from GTFfilenames
		file gtf from stringtieMerge_gtf
		file assembled from assembled_gtfs.toList()
		
	output:
		file "merged.gtf" into merged_gtf,merged_gtf_quant
		file ".command.sh" into stringtieMerge_commands
	script:
		"""
		stringtie --merge -p 14 -G $gtf -o merged.gtf $gtf_filenames		
		"""
}


process gffcompare{
	publishDir "Assembly", mode: "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "gffcompare_commands.sh" : "$filename" }

	input:
		file gtf from stringtieCompare_gtf
		file merged from merged_gtf

	output:
		file "merged.*" into gtf_comparison
		file ".command.sh" into gffcompare_commands
	
	script:
		"""
		gffcompare -r $gtf -G -o merged $merged
		"""

}


/*
*quantifying transcripts
*/

process stringtieFPKM {
	publishDir "Ballgown", mode: "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }

	input:
		file bam from filtered_bams_quantify
		file gtf from merged_gtf_quant.collect()		

	output:
		file "${bam.baseName}" into ballgown
		file "${bam.baseName}_gene_abund.txt" into gene_abund
		file ".command.sh" into stringtieFPKM_commands
	script:
		prefix = bam[0].toString() - ~/_filtered.*/
		"""
		stringtie -p 14 -A ${bam.baseName}_gene_abund.txt -b ${bam.baseName} -G $gtf $bam
		"""

}
/*
* multiqc
*/

/*
*Combine all the trimm logs into one file
*/

trim_summary
	.collectFile(name: "Trimmomatic_summaries.txt",newLine: true)
	.set{All_trim_summaries}

fastqc_commands
                .concat(trimmomatic_commands,trim_fastqc_commands,star1pass_index_commands,star1pass_commands,star2pass_index_commands,star2pass_commands,filtering_commands,
                stringtieAssembly_commands,stringtieMerge_commands,gffcompare_commands,stringtieFPKM_commands)
                .collectFile(name: "RNA-Seq_commands.txt",newLine: true)
                .set{commands}

process multiqc {
	publishDir "Reports", mode: "copy"
	
	input:
		file fastqc from fastqc_results.collect()
		file star_log from star2pass_logs.collect()
		file trim_summary from All_trim_summaries.collect()
		file commands 
	
	output:
		file "multiqc_report.html" into multiqc_report
		file "RNA-Seq_commands.txt"

	script:
		"""
		multiqc .
		"""
}
