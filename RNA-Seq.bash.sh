#!/bin/bash -ue

samples=samples.txt
threads=14
length=150
##length is length of reads
genome=genome.fa
gtf=genes.gtf
deg=true
splicing=false
fastq=Complete_Set		
comparisons=comparisons.txt
strand=""
##strandness:
## "" -unstranded
##rf -reversely stranded - use for illumina directional libraries
##fr stranded
se=false
##if data is paired end, set se to false
outDir=RNA-Seq
bigwigs=true
##select true to generate bigwigs
annotation=mm10
#####

#################
##testing if options worked
##########################################





if [ -z "${samples:-}" ]
then
	echo "please specify a samples file"
	exit
elif [ ! -e $samples ]
then
        echo "sample file $samples not found"
        exit 1

fi



if [ -z "${comparisons:-}" ] && [ "$deg" = "true" ]
then
	echo "please specify a comparisons file"
	exit 1

elif [ ! -e $comparisons ] && [ "$deg" =  "true" ]
then
        echo "comparisons file $comparisons not found"
        exit 1

fi



if [ -z "${genome:-}" ]
then
	echo "please specify a reference genome fasta file"
	exit 1

elif [ ! -e $genome ]
then
	echo "Genome fasta file not found"
	exit 1

fi


if [ -z "${gtf:-}" ]
then
	echo "please specify a gtf file"
elif [ ! -e $gtf ]
then

	echo "gtf file not found"
	exit 1

fi


if [ -z "${length:-}" ]
then
	echo "please specify the length of your reads"
	exit 1
fi


if [ -z "${threads:-}" ]
then

	threads=1
fi


if [ -z "${outDir:-}" ]
then
	now=$(date +"%m_%d_%Y")
	outDir=RNA-Seq_out_${now}
	echo "No output directory specified. Using output directory $outDir"

fi



if [ -z "${splicing:-}" ]
then

        splicing=false
        echo "Splicing Analysis will not be done"
fi


if [ -z "${deg:-}" ]
then

        deg=false
        echo "Differential gene expression analysis will not be done"
fi


if [ "$deg" != "true" ] && [ "$splicing" != "true" ]
then
        echo "Only alignment with STAR will be done"
fi



if [ -z "${strand:-}" ]
then
	echo "strandness not specified. setting strandness to unstranded"
	strand=""	
fi

if [ -z "${se:-}" ]
then

	echo "data is single-end: $se "
	se=false
fi 


if [ -z "${fastq:-}" ]
then
	echo "Please specify the location of your fastq files"
	exit 1
fi


#####################
##set up directories
#can be changed


fastqc=($outDir/FastQC)
trimmed_fastq=($outDir/Trimmed_fastq)
trimmed_fastqc=($outDir/Trimmed_fastQC)
star1=($outDir/STAR1pass)
star1index=($outDir/STAR1index)
star2=($outDir/STAR2pass)
star2index=($outDir/STAR2index)
cufflinks=($outDir/Cufflinks)
counts=($outDir/Counts)
assembly=($outDir/Assembly)
ballgown=($outDir/Ballgown)
deseq=($outDir/DESeq2)
bw_dir=($outDir/BigWigs)
multiqc=($outDir/Multiqc)
enrich=($outDir/Cluster_Analysis)

###############step1. FastQC


if [ ! -d $fastqc ]
then
	mkdir -p $fastqc
fi



fastqc -t 10 -o $fastqc $fastq/*gz


######## Step2. Remove adaptors


if [ ! -d $trimmed_fastq ]
then
	mkdir -p $trimmed_fastq
fi


while read -r sampleName coreNumber group r1 r2
do
	if [ $se != "true" ]
	then

	read1=($fastq/${r1})
	read2=($fastq/${r2})
	i=$coreNumber
	echo trimming adaptors from $i
	 java -jar /usr/share/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $threads -phred33 -trimlog $trimmed_fastq/Trimmomatic_${i}.log \
         $read1 $read2 $trimmed_fastq/${i}_1_at.fq.gz $trimmed_fastq/${i}_1_unpaired.fq.gz $trimmed_fastq/${i}_2_at.fq.gz $trimmed_fastq/${i}_2_unpaired.fq.gz \
         ILLUMINACLIP:/usr/share/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 2> $trimmed_fastq/${i}_trimStats.txt
	
	gzip $trimmed_fastq/Trimmomatic_${i}.log
	else
	
	read1=($fastq/${r1})
        i=$coreNumber
        echo trimming adaptors from $i
         java -jar /usr/share/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $threads -phred33 -trimlog $trimmed_fastq/Trimmomatic_${i}.log \
         $read1  $trimmed_fastq/${i}_1_at.fq.gz \
         ILLUMINACLIP:/usr/share/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 2> $trimmed_fastq/${i}_trimStats.txt
	
	gzip $trimmed_fastq/Trimmomatic_${i}.log
	fi

done < $samples




###### FastQC again



 fastqc -t 10 -o $trimmed_fastq $trimmed_fastq/*at.fq.gz


## STAR 1 pass

if [ ! -d $star1index ]
then
	mkdir -p $star1index
fi


if [ ! -d $star1 ]
then
	mkdir -p $star1

fi

echo Generate index for STAR first pass
STAR --runThreadN $threads \
                 --runMode genomeGenerate \
                --genomeDir $star1index \
                --genomeFastaFiles $genome \
                --sjdbGTFfile $gtf \
                --sjdbOverhang $length




for i in $(ls $trimmed_fastq/*[12]_at.fq.gz | cut -f3 -d "/" | sed 's/_[12]_at.fq.gz//g' | sort | uniq)
do
       
	reads=$(ls $trimmed_fastq/${i}*_at.fq.gz)
	echo aligning $i for the first time
	STAR --runThreadN $threads \
                --genomeDir $star1index \
                --readFilesIn $reads \
                --readFilesCommand gunzip -c \
                --runMode alignReads \
                --outSAMattributes All \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --outSAMstrandField intronMotif \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix $star1/$i


done


###### STAR 2 pass

if [ ! -d $star2index ]
then
	mkdir -p $star2index
fi

if [ ! -d $star2 ]
then
	mkdir -p $star2
fi

SJ_gtfs=$(ls $star1/*SJ.out.tab)

STAR --runThreadN $threads \
                --runMode genomeGenerate \
                --genomeDir $star2index \
                --genomeFastaFiles $genome \
                --sjdbFileChrStartEnd $SJ_gtfs \
                --sjdbGTFfile $gtf \
                --sjdbOverhang $length



for i in $(ls $trimmed_fastq/*[12]_at.fq.gz | cut -f3 -d "/" | sed 's/_[12]_at.fq.gz//g' | sort | uniq)
do

		reads=$(ls $trimmed_fastq/${i}*_at.fq.gz)
	        echo aligning $i for the second time
		
                STAR --runThreadN $threads \
                --genomeDir $star2index \
                --readFilesIn $reads \
                --readFilesCommand gunzip -c \
                --runMode alignReads \
                --outSAMattributes All \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --outSAMstrandField intronMotif \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix $star2/$i


		if [ $se != "true" ]
		then

		samtools view -b -f 0x2 $star2/${i}Aligned.sortedByCoord.out.bam > $star2/${i}_properly_paired_sorted.bam
                samtools index $star2/${i}_properly_paired_sorted.bam

		else
                
	        mv $star2/${i}Aligned.sortedByCoord.out.bam $star2/${i}_sorted.bam 
		samtools index $star2/${i}_sorted.bam 
	

		fi
done


rm -r $star1 $star1index

if [ "$bigwigs" = "true" ]
then
	if [ ! -d $bw_dir ]
	then
		mkdir -p "$bw_dir"
	fi


	for file in $star2/*_sorted.bam
	do

		sample=$(echo $file | cut -f3 -d "/" | cut -f1 -d "_")
		
		bamCoverage -b $file --normalizeUsing CPM -of bigwig -o $bw_dir/${sample}.bw
	done

fi



##get expression using cufflinks

if [ ! -d $cufflinks ]
then

	mkdir $cufflinks
fi

lib=""



if [ "$strand" = "rf" ]
then
lib="--library-type fr-firststrand"

elif [ "$strand" = "fr" ]
then

lib="--library-type fr-secondstrand"

fi



if [ "$strand" = "rf" ]
then
	lib="--library-type fr-firststrand"
elif [ "$strand" = "fr" ]
then
	lib="--library-type fr-secondstrand"
fi	

for file in $star2/*sorted.bam
do
	echo "running cufflinks on $file"

	sample=$(echo $file | cut -f3 -d "/" | cut -f1 -d "_")
	cufflinks -p $threads -o $cufflinks/$sample \
           -G $gtf \
           -b $genome \
           -u \
           $lib \
           $file 2> $cufflinks/${sample}_stderr.txt
done


#combine
combineCufflinksFpkms.R in=$cufflinks




if [ "$deg" = "true" ]
then


#feature Counts

	if [ ! -d $counts ]
	then

		mkdir -p $counts
	fi

	fc_strand=""
	if [ $strand = "rf" ]
	then
	fc_strand="-s 2"

	elif [ $strand = "fr" ]
	then

	fc_strand="-s 1"

	fi


	if [ "$se" != "true" ]
	then            
                
		featureCounts -a $gtf $fc_strand -p -o $counts/featureCounts_gene_counts.txt $star2/*sorted.bam
                
	else            
		featureCounts -a $gtf $fc_strand -o $counts/featureCounts_gene_counts.txt $star2/*sorted.bam
                
	fi              
                        
       ##DESeq2         
                        
       ##add a prepfeatureCounts
                        
       if [ ! -d $deseq ]
       then             
                        
       	mkdir -p $deseq 
       fi               
                        
	prepfeatureCounts.R --counts=${counts}/featureCounts_gene_counts.txt --out=${counts}
	                
	cut -f2,3 $samples > $deseq/sample_sheet_for_debrowser.txt
                        
	DESeq2.R --counts=${counts}/featureCounts_for_DESeq2.csv --annotation=${counts}/gene_lengths.csv --species=${annotation} --fpkms=${cufflinks}/cufflink_fpkms_all_samples.csv --samples=${samples} --comparisons=${comparisons} --out=${deseq}

	if [ ! -d $enrich ]
	then

		mkdir -p $enrich

	fi

        ClusterProfiler.R --DE=${deseq} --org=${annotation} --out=${enrich}
fi                      
                        
                        
if [ "$splicing" = "true" ]
then                    
                        
                        
	if [ ! -d $assembly ]
	then            
		mkdir $assembly
                        
	fi              
	                
                        
	if [ $se != "true" ]
	then


	for i in $star2/*filtered_sortedByCoord.out.bam
	do
		sample=$(basename $i _filtered_sortedByCoord.bam)
		
		stringtie -p 14 -G $gtf -o $assembly/${sample}.gtf $i
	done

	else
	for i in $star2/*Aligned_sortedByCoord.out.bam
        do
                sample=$(basename $i Aligned_sortedByCoord.bam)
                
                stringtie -p 14 -G $gtf -o $assembly/${sample}.gtf $i
        done
	fi

	

	ls $assembly/*gtf | tr " " "\n" > $assembly/gtf_files.txt

	gtf_filenames=$assembly/gtf_files.txt

	stringtie --merge -p 14 -G $gtf -o $assembly/merged.gtf $gtf_filenames

	gffcompare -r $gtf -G -o $assembly/merged $assembly/merged.gtf
        gffCompare_parse.py $assembly/merged.stats > $assembly/parsed_stats.txt


	if [ $se != "true" ]
	then


	for i in $star2/*filtered_sortedByCoord.out.bam
	do
		sample=$(basename $i _filtered_sortedByCoord.bam)
		
		
		stringtie_dir=""
		if [ $strand = "fr" ]
		then

			strintie_dir="--fr"
		else
			stringtie_dir "--rf"

		fi
		
		stringtie -p $threads $stringtie_dir -e -B -o ballgown/${sample}/${sample}.gtf -G $gtf $bam

	done
	
	else

	for i in $star2/*Aligned_sortedByCoord.out.bam
        do
                sample=$(basename $i Aligned_sortedByCoord.bam)
                

		if [ $strand = "fr" ]
                then

                        strintie_dir="--fr"
                else
                        stringtie_dir "--rf"

                fi
                
                stringtie -p $threads $stringtie_dir -e -B -o ballgown/${sample}/${sample}.gtf -G $gtf $bam

        done

	fi
	## add ballgown
fi



######################
if [ ! -d $multiqc ]
then
	mkdir $multiqc
fi

cd $multiqc
ln -s ../../$trimmed_fastq/*trimStats.txt .
ln -s ../../$star2/*final.out .
ln -s ../../$fastqc/*  .
ln -s ../../$trimmed_fastq/*fastqc* .
ln -s ../../$assembly/*tsv .
ln -s ../../$counts/*summary .
ln -s ../../$deseq/*mq* .
ln -s ../../$deseq/DE_summary.txt


#multiqc -c ../../multiqc_config.yaml .
cd -

#get software versions
fastqc --version &> $outDir/softwareVersions.txt
which trimmomatic-0.39.jar | cut -f5 -d "/" | cut -f1,2 -d "." &>> $outDir/softwareVersions.txt
STAR --version &>> $outDir/softwareVersions.txt
samtools --version | grep "samtools" &>> $outDir/softwareVersions.txt
echo stringtie $(stringtie --version) &>> $outDir/softwareVersions.txt
gffcompare --version &>> $outDir/softwareVersions.txt
