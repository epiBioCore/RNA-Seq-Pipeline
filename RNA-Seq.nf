#!/usr/bin/env nextflow

/*
*
*======================================================================================
*
*	UND Bioinformatics Core RNA-Seq Analysis Pipeline
*
*=====================================================================================
*
*author: Danielle Perley 
*email: danielle.perley@med.und.edu
*
*------------------------------------------------------------------------------------
*/


def helpMessage() {
	log.info"""

	usage: nextflow run RNA-Seq.nr --reads 'Raw_fastq/*_{1,2}.fq' --genome genome.fa --gtf genes.gtf --samples samples.txt --comparisons comparisons.txt

	Mandatory arguments:
	--reads		path to fastq files. Must be surrounded in quotes. ex 'Raw_fastq/*{1,2}.fq.gz', 'data/*_{R1,R2}*.fastq.gz'
	--genome	path to genome fasta file
	--gtf		path to gtf file

	Optional arguments:
	--singleEnd		specifies the fastq files are single end
	--length		the average length of reads. Default 100
	--strand		strandness of experiment. Acceptable values are "unstranded","first","second". Default: "unstranded".
	--threads		number of threads to use. Default: 1	
	--feautureCounts	perform read counting using featureCounts. Default true 
	--stringtie		perform transcript assembly using stringtie. Default true
	--noDE			do not perform Differential Expression analysis. 
        --samples       	a tab delimited file with sample names and treatment groups. More information provided below
        --comparisons   	a tab delimited file with comparisons. More information provided below.

	More Information:
	If performing a differential expression analysis, two tab-delimited files must be specified. A samples file that contains the coreNumber, sampleName, and treatment group. ex. A 
	A comparisons file with the desired comparisons. One comparison per row. The first column is the treatment group and the second group is the control. DESeq2 will calculate log Fold changes
	relative the the control group. ex.
	"""
} 


// variables

params.reads = 'Raw_fastq/*{1,2}.fq.gz'
params.length = 100
params.gtf = false
params.genome = false
params.singleEnd = false
params.threads = 1
params.samples = false
params.comparisons = false
params.strand = "unstranded"
params.featureCounts = true
params.stringtie = true
params.noDE = false

params.help = false
if (params.help){
    helpMessage()
    exit 0
}

Channel
	.fromFilePairs(params.reads, size: params.singleEnd ? 1:2)
	.ifEmpty { exit 1, "Fastq files not found: ${params.reads}\nPath must be enclosed with quotes. \nPath must contain at least one wildcard \n if this is single end data, please specify --singleEnd"}
	.into { read_files; read_files_fastqc }

// Validate inputs

if (params.gtf) {
	gtf = file(params.gtf)
	if ( !gtf.exists() ) {
		 exit 1, "GTF annotation file not found: ${params.gtf}" 
        } else {
	Channel
		.fromPath(params.gtf)
		.into { star1_gtf; star2_gtf; stringtie_gtf; stringtieMerge_gtf; stringtieCompare_gtf; stringtieFPKM_gtf; featureCounts_gtf}
       }
} else {
	exit 1, "No GTF annotation file specified"
}

if (params.genome) { 
        genome= file(params.genome)
	if ( !genome.exists() ) {
	exit 1, "Genome Fasta file not found: ${params.genome}."
	} else {
	Channel
		.fromPath(params.genome)
		.into {star1_fasta; star2_fasta; fasta_check }
        }
	fasta_check.println()
} else {
	exit 1, "No fasta file specified"
}


if (params.samples) {
	samples = file(params.samples)
	if ( !samples.exists() ) {
		exit 1, "Sample sheet not found: ${params.samples}"
	} else {
	Channel
		.fromPath(params.samples)
		.into { mq_samples; samples }
	}
} else if (!paramss.samples) {
	exit 1, "No sample sheet specified!"
}

if (!params.noDE && params.comparisons) {
	comparisons = file(params.comparisons)
	if ( !comparisons.exists() ) {
		exit 1, "List of Comparisons not found: ${params.comparisons}"
	} else {
	Channel
		.fromPath(params.comparisons)
		.set { comparisons}
	}
} else if (!params.noDE && !params.comparison) {
	exit 1, "No comparisons specified"
}		


pca_header = file("bin/pca_header.txt")
heatmap_header = file("bin/heatmap_header.txt")
multiqc_config = file("bin/multiqc_config.yaml")




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
		prefix = fastq[0].toString() - ~/_.*/
                "fastqc $fastq"
}


/*
* Trimm adaptors using Trimnomatic
*/

process adapterTrimming {
	publishDir "Trimmed_fastq", mode: "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${name}_commands.sh" : "$filename" }

	input:
		set val(name),file(read_files) from read_files

	output:
		file "*at.fq.gz" into read_files_star1,read_files_star2,trimmed_reads_fastqc
		file ".command.err" into trim_summary
		file "*.log" into trim_logs
		file ".command.sh" into trimmomatic_commands


	script:
		prefix1 = read_files[0].toString() - ~/\.f*q.gz/
		prefix2 = read_files[1].toString() - ~/\.f*q.gz/

		if (params.singleEnd) 
		"""
                trimmomatic-0.33.jar SE -threads 5 -phred33 -trimlog Trimmomatic_${name}.log \\
                $read_files ${prefix1}_at.fq.gz \\
                ILLUMINACLIP:/usr/share/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10
                """
		else 
		"""
		trimmomatic-0.33.jar PE -threads 5 -phred33 -trimlog Trimmomatic_${name}.log \\
		$read_files ${prefix1}_at.fq.gz ${prefix1}_unpaired.fq.gz ${prefix2}_at.fq.gz ${prefix2}_unpaired.fq.gz \\
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
		prefix = reads[0].toString() - ~/_.*/
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
   
		STAR --runThreadN ${params.threads} \\
		 --runMode genomeGenerate \\
		--genomeDir STAR1passIndex \\
		--genomeFastaFiles $fasta \\
		--sjdbGTFfile $gtf \\
		--sjdbOverhang ${params.length}
		"""
}



/*
* Aligning with STAR for 1 pass
* I will want to buffer the input for the read_files_star1
*/

process STAR1pass {
	maxForks 1

	input:
		file reads from read_files_star1
		file index from star_index.collect()

	output:
		file "*SJ.out.tab" into star1pass_gtf
 		file "*.bam" into star1pass_aln
		file ".command.sh" into star1pass_commands
	script:
		
		prefix = reads[0].toString() - ~/_.*/
		"""
		STAR --runThreadN ${params.threads} \\
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

	 STAR --runThreadN ${params.threads} \\
                --runMode genomeGenerate \\
                --genomeDir STAR2passIndex \\
                --genomeFastaFiles $fasta \\
		--sjdbFileChrStartEnd $SJ_gtfs \\
                --sjdbGTFfile $gtf \\
                --sjdbOverhang ${params.length}
      """

}

/*
* align bam files again with updated genome index
*/

process STAR2pass {
	publishDir "STAR2pass", mode : "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }	
	maxForks 1
	input:
		file reads from read_files_star2
		file index from Star2pass_index.collect()

	output:
		file "*.bam" into star2pass_aln
		file "*.final.out" into star2pass_logs
		file ".command.sh" into star2pass_commands
	
	script:
		prefix = reads[0].toString() - ~/_.*/
                """
                STAR --runThreadN ${params.threads} \\
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
* only if data is paired end
*/

if (!params.singleEnd) {
process filterBam {
	publishDir "STAR2pass", mode : "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }

	input:
		file bam from star2pass_aln

	output:
		file "${prefix}_filtered_sortedByCoord.out.bam" into filtered_bams
		file "*.bai" into bam_indices
		file ".command.sh" into filtering_commands
	script:
		prefix = bam.toString() - ~/Aligned.sortedByCoord.out.bam/
		"""
		samtools view -b -f 0x2 $bam > ${prefix}_filtered_sortedByCoord.out.bam
		samtools index ${prefix}_filtered_sortedByCoord.out.bam
		"""
}
}


/*
* if fastq files are paired-end, then the filtered bams will be input into the stringtie assembly,
* otherwise the bam files from star
*/

if (params.singleEnd) {
	star2pass_aln
		.into{bams_assembly; bams_quantify; bams_featureCounts }
} else {
	filtered_bams
		.into {bams_assembly; bams_quantify; bams_featureCounts }	
	
}

/*
*stringtie assembly
*/

if (params.stringtie) {
process stringtieAssembly {
	publishDir "Assembly", mode: "copy",
 		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }
	input:
		file bam from bams_assembly
		file gtf from stringtie_gtf.collect()

	output:
		file "*gtf" into assembled_gtfs, assembled_transcripts_fn
		file ".command.sh" into stringtieAssembly_commands

	script:
		prefix = bam.toString() - ~/_filtered.*|Aligned.*/
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
  		file "*for_multiqc.*" into gtf_multiqc
		file ".command.sh" into gffcompare_commands
	
	script:
		"""
		gffcompare -r $gtf -G -o merged $merged
		gffCompare_parse.py merged.stats
		"""

}


/*
*quantifying transcripts
*/

process stringtieFPKM {
	publishDir "Ballgown", mode: "copy",
		saveAs: {filename -> filename.indexOf(".sh") > 0 ? "${prefix}_commands.sh" : "$filename" }

	input:
		file bam from bams_quantify
		file gtf from merged_gtf_quant.collect()		

	output:
		file "$prefix" into ballgown
		file "$prefix/${prefix}.gtf" into sample_gtfs
		file ".command.sh" into stringtieFPKM_commands
	script:
		prefix = bam[0].toString() - ~/_filtered.*|Aligned.*/
		def stringtie_dir = ""
		if (params.strand =~/first/) {
			stringtie_dir = "--fr"
		} else if (params.strand =~/second/) {
			stringtie_dir = "--rf"
		}
		"""
		stringtie -p ${params.threads} $stringtie_dir -e -B -o ./$prefix/${prefix}.gtf -G $gtf $bam
		"""

}

/*
*sample_gtfs
*        .collectFile() {file -> ['GTFlist.txt',file.name + '\n'] }
*        .set {sampleGTFs}
*/


sample_gtfs
        .collectFile() {file -> ['GTFlist.txt', file.name.replaceFirst(~/\.gtf/,'')+ '\t' + file.name.replaceFirst(~/\.gtf/,'') +'/' + file.name + '\n'] }
        .set {sampleGTFs}



/*
*prepDE
*/

process prepDE {
	publishDir "Counts", mode: "copy"

	input:
		file ballgown from ballgown.collect()
		file sampleGTFs
	output:
		file "stringtie_gene_counts.csv" into stringtie_gene_counts
		file "stringtie_transcript_counts.csv" into transcript_counts
		file ".command.sh" into prepDE_commands

	script:
		"""
 		prepDE.py -i $sampleGTFs -g stringtie_gene_counts.csv -t stringtie_transcript_counts.csv -l ${params.length}
		"""

}

}


/* 
*featureCounts
*/

if (params.featureCounts) {
process featureCounts {
	publishDir "Counts",mode: "copy"

	input:
		file gtf from featureCounts_gtf
		file bam from bams_featureCounts.collect()

	output:
		file "featureCounts_gene_counts.csv" into featureCounts_gene_counts
		file "*.summary" into featureCounts_stats

	script:
		def fc_strand = 0
		if (params.strand =~ /first/ ) {
		   fc_strand = 2
		} else if (params.strand =~ /second/) {
		   fc_strand = 1
		}
		
	        """ 
		featureCounts -a $gtf -p -s $fc_strand -o featureCounts_gene_counts.csv $bam 
		"""
}

}

if (params.featureCounts) {
	featureCounts_gene_counts
			.set { DESeq2_gene_counts}
	print("Setting DESeq2 gene counts")			
} else {
	stringtie_gene_counts
			.set { DESeq2_gene_counts }
	print("using stringtie counts")
}

/*
*DESeq2
*/

if (!params.noDE) {
process DESeq2 {
	publishDir "DESeq2", mode: "copy"

	input:
		file DESeq2_gene_counts
		file comparisons
		file samples
		file pca_header
		file heatmap_header

	output:
		file "sample_distance_matrix_mqc.csv" into heatmap_config
		file "pca_config_mqc.yaml" into pca_config
		file "DE_summary.txt" into DE_results
		file "*.{txt,rda,csv,pdf}" into DESeq2_results

	script:
		"""
		DESeq2.R --counts=$DESeq2_gene_counts --samples=$samples --comparisons=$comparisons
		cat $pca_header pca_for_mq.yaml > pca_config_mqc.yaml
		cat $heatmap_header sample_distance_matrix_for_mq.csv > sample_distance_matrix_mqc.csv
		"""

}
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

/*
*fastqc_commands
*                .concat(trimmomatic_commands,trim_fastqc_commands,star1pass_index_commands,star1pass_commands,star2pass_index_commands,star2pass_commands,filtering_commands,
*                stringtieAssembly_commands,stringtieMerge_commands,gffcompare_commands,stringtieFPKM_commands,prepDE_commands)
*                .collectFile(name: "RNA-Seq_commands.txt",newLine: true)
*                .set{commands}
*/
process multiqc {
	publishDir "Reports", mode: "copy"
	
	input:
		file fastqc from fastqc_results.collect()
		file star_log from star2pass_logs.collect()
		file All_trim_summaries
		file featureCounts_stats
		file heatmap_config
		file pca_config
                file gtf_multiqc
                file multiqc_config
		file DE_results 
		file mq_samples  
	
	output:
		file "multiqc_report.html" into multiqc_report
	

	script:
		"""
		cut -f1-3 $mq_samples > samples_for_mq.txt
		multiqc -f --sample-names samples_for_mq.txt .
		"""
}

/*
* get software versions
*/
process softwareVersions {
	publishDir "Reports", mode: "copy"
	
	output:
		file "softwareVersions.txt" into softwareVersions

	script:
		'''
		fastqc --version &> softwareVersions.txt
		which trimmomatic-0.33.jar | cut -f5 -d "/" | cut -f1,2 -d "." &>> softwareVersions.txt
		STAR --version &>> softwareVersions.txt
		samtools --version | grep "samtools" &>> softwareVersions.txt
		echo stringtie $(stringtie -v) &>> softwareVersions.txt
		gffcompare --version &>> softwareVersions.txt
		'''
}	
