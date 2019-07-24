#!/usr/bin/env python

import os
import re 


infile = open("merged.stats",'r')


#sp_data,mn_data = parse_gffcomp(infile)
#
#for k in sp_data:
#        outfile1.write(k,"\t",sp_data[k],"\n")
#
#
#for k in mn_data:
#        outfile2.write(k,"\t",mn_data[k],"\n")
#


def parse_gffcomp(infile):        
    ##initalize dicts sensitivity and precison
    sp_data=dict()
    
    #and missing and novel data
    mn_data=dict()
     ##parse summary
    # ##all loci
    # for line in log:
    #     query=re.search(r'Query mRNAs\s+:\s+(\d+) in\s+(\d+) loci',line )
    #     if query:
    #         parsed_data['query_mRNAs'] = int(query.group(1))
    #         parsed_data['total_loci'] = int(query.group(2))
    # ##ref loci    
    # for line in log:
    #     query=re.search(r'Reference mRNAs\s+:\s+(\d+) in\s+(\d+) loci',line )
    #     if query:
    #         parsed_data['known_mRNAs'] = int(query.group(1))
    #         parsed_data['know_loci'] = int(query.group(2))
    #         
    ###parse sensitivity and precision        
     
    #####make regexs
    ##SandP=sensitivity and Precision table
    ##example regex
    ##r'Base level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)'
    ##to match:
    ##         Base level:   100.0     |    34.9    |
    ## in plain english:
    ## match "Base level: 
    ## one or more spaces \s+
    ##capture the floating point number (\d+\.\d)
    ##one or more spaces \s+
    ##capture the floating point number (\d+\.\d)
    ## the sensitivity is group(1) and precision is group(2)
    ##
    ##Missing_novel= the Missed and Novel percentages
    ##example regex     
    ##r'Missed exons:\s+(\d+/\d+)\s+\(\s+(\d+\.\d%)'
    ##to match:
    ##           Missed exons:       0/135775  (  0.0%)
    ## in plain english:
    ## match: "Missed exons:"
    ## one or more spaces: \s+
    ##capture the fraction (\d+/\d)
    ##one or more spaces \s+
    ##a parenthesis \(
    ##one or more spaces \s+
    ##capture the percentage (\d+\.\d%)   I may decide to not capture the "%"
    ## depending on multiqc formating     
    ##
    ## (I don't think I will capture matching stats)
    ##  Matching:
    #     Matching intron chains:   14476        
    regexes={'SandP':{
                    'base':r'Base level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)',
                    'exon':r'Exon level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)',
                    'intron':r'Intron level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)',
                    'intron_chain':r'Intron chain level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)',
                    'transcript':r'Transcript level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)',
                    'locus':r'Locus level:\s+(\d+\.\d)\s+\|\s+(\d+\.\d)'},
            'Missing_Novel':{
                    'missed_exons':r'Missed exons:.*\(\s+(\d+\.\d)%',
                    'novel_exons':r'Novel exons:.*\(\s+(\d+\.\d)%',
                    'missed_introns':r'Missed introns:.*\(\s+(\d+\.\d)%',
                    'novel_introns':r'Novel introns.*\(\s+(\d+\.\d)%',
                    'missed_loci':r'Missed loci:.*\(\s+(\d+\.\d)%',
                    'novel_loci':r'Novel loci:.*\(\s+(\d+\.\d)%'}
    #         'Matching':{
    #                 'matching_intron_chains':r'Matching intron chains:\s+(\d+)',
    #                 'matching_transcripts':r'Matching transcripts:\s+(\d+)',
    #                 'matching_loci':r'Matching loci:\s+(d+)'}
    }
    
   
    ##the sensitivity and precision table
    for line in infile:
        for k,r in regexes['SandP'].items():
	    query=re.search(r,line)
            if query:
                sp_data[k]=query.group(1) + "\t" + query.group(2)
                
            
    infile.seek(0)	    
    #the matching stats
    #for line in log:
    #    for k,r in regexes['Matching'].items():
    #        query=re.search(r,line) 
    #        if query:
    #            parsed_data[k]=query.group(1)            
    
    ##the missing and novel stats
    for line in infile:
        for k,r in regexes['Missing_Novel'].items(): 
		query=re.search(r,line)
                if query:
                    mn_data[k]=query.group(1)
    
    	
    ##final values
    ##total union
    #for line in log:
    #    query=re.search(r'Total union super-loci across all input datasets:\s+(\d+)',line)
    #    if query:
    #        parsed_data['total_union']=query.group(1)
    #        
    #        
    #for line in log:
    #    query=re.search(r'(\d+) out of (\d+) consensus transcripts written in .*gtf \((\d+) discarded',line)
    #    if query:
    #        parsed_data['written_consensus_transcripts']=query.group(1)
    #        parsed_data['total_consensus_transcripts']=query.group(2)
    #        parsed_data['redundant_consensus_transcripts']=query.group(3)        
    # 
    
    
    return sp_data,mn_data
 


sp_data,mn_data = parse_gffcomp(infile)


#for k in sp_data:
#	print k

#### output data
##headers
outfile1 = open("sensitivity_and_precision_for_multiqc.tsv",'w')
outfile1.write("Context\tSensitivity\tPrecision\n")
outfile1.close()


outfile1=open("sensitivity_and_precision_for_multiqc.tsv",'a')

for k in sp_data:
        outfile1.write(k + "\t" + sp_data[k] +"\n")


##Missing and Novel
outfile2 = open("missing_novel_for_multiqc.tsv",'w')
outfile2.write("Context\tMissing\tNovel\n")
outfile2.close()

outfile2 = open("missing_novel_for_multiqc.tsv",'a')
outfile2.write("exons\t" + str(mn_data['missed_exons'])+ "\t" + str(mn_data['novel_exons']) +"\n" +
                "introns\t" + str(mn_data['missed_introns']) + "\t" + str(mn_data['novel_introns']) + "\n" +
	        "loci\t" + str(mn_data['missed_loci']) + "\t" + str(mn_data['novel_loci']) + "\n")

infile.close()
outfile1.close()
outfile2.close()
