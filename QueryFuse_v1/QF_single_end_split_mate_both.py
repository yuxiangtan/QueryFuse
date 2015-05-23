# -*- coding: cp936 -*-
"""
This script is a Gene_specific_fusion_query subfunction: single_end pocessing. Specially for step (6) split mate to both.
This script only process one end at a time. To process both, need to run for each end.
This script will finally generate an ordered and merged file that can be used to summarize and count spliting reads which the single aligned end is on query. 

=============================
Usage: python QF_single_end_split_mate_both.py
-h help

-i file_prefix that all the defualt can use                                                                         *[No default value]

-E end_info_in_English                                                                                              *[No default value, but should be first or second]

-q query_split_ID.txt                                                                                               [default value: file_prefixfolder, unmapped.bam_END_N_mate_on_query.psl_split_ID_subtract_ID.txt]

-Q to_both_ID(from the other mate) - location of the to_both_ID.txt(should use from the other mate)	            [default value: file_prefixfolder, singleton_sorted_to_to_both_query_and_other_ID.txt]

-w whole_gene_list.bed				                                                                    *[No default value]

-u Unmap_over_query_split_subtract_BED- Bed file which contain subtract information for unmap over query of the mate same as the query_split_ID         [default value file_prefixfolder, unmapped.bam_END_N_mate_on_query.psl_split_ID_subtract.bed]

-m mate_fa - fa file of the same mate as the query_split_ID                                                         [default value file_prefixfolder, unmapped.bam_END_N_mate.fa]

-U UNMAP_ON_QUERY_PSL - location of UNMAP_ON_QUERY_PSL                                                              [default value file_prefixfolder, unmapped.bam_END_N_mate_on_query.psl]
  
-t tophat_genome_reference - indexed genome_reference.		                                                    *[No default value]

-l Read_length - length of reads		                                                                    [default value 99]

-s SINGLE_ON_QUERY_BED		                                                                                    [default value file_prefixfolder, singleton_sorted_to_query.bed]

-g LOG_folder                                                                                                       *[No default value]

-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]

-F QueryFuse_path                                                                                                   *[No default value]

-a Align_percent: min percentage of alignment                                                                       [default value is 98]

============================

Python & Module requirement:
Versions: 2.7 or above
Module: No additional Python Module is required.

============================
Library file requirement:
Not Standalone version, few library file is required.
============================

"""

##By Yuxiang Tan
##Contact: yuxiang.tan@gmail.com
##Compatible Python Version:2.7 or above

###Code Framework

if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re, subprocess, time 
	import os
	import math
	from itertools import ifilter
	import QF_all_modules
	OUTPUT_SEP_CHAR='\t'
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 13:
		print __doc__
		sys.exit(3)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	file_prefix=None
	whole_gene_list=None
	tophat_genome_ref=None
        LOG_F=None
        END_N=None
	QF_path=None
        
        query_split=None
	single_to_both=None
	unmap_sub_bed=None
	unmap_query_psl=None
	read_len=99
	single_query_bed=None
        mate_fa=None
        resume_stat=0
        Align_percent=98
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:w:g:t:E:F:q:Q:u:o:U:s:m:l:r:R:P:Y:a:S:z')
        for opt in optlist:
            if opt[0] == '-h':
                print __doc__; sys.exit(2)
            elif opt[0] == '-i': file_prefix = opt[1]
            elif opt[0] == '-w': whole_gene_list = opt[1]
            elif opt[0] == '-g': LOG_F = opt[1]
            elif opt[0] == '-t': tophat_genome_ref =opt[1]
            elif opt[0] == '-E': END_N = opt[1]
            elif opt[0] == '-F': QF_path = opt[1]
            elif opt[0] == '-q': query_split =opt[1]
            elif opt[0] == '-Q': single_to_both= opt[1]
            elif opt[0] == '-u': unmap_sub_bed = opt[1]
            elif opt[0] == '-U': unmap_query_psl =opt[1]
	    elif opt[0] == '-s': single_query_bed = opt[1]
            elif opt[0] == '-m': mate_fa =opt[1]
            elif opt[0] == '-l': read_len = int(opt[1])
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-a': Align_percent =opt[1]
	    
        if LOG_F==None:
            print "Warning: LOG_F is not provided in QF_single_end_split_mate_both.py, exit."; sys.exit(2)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
        LOG_WHO=LOG_F+"/whole.log"
        log_error=open(LOG_ERR,"a")
        log_whole=open(LOG_WHO,"a")
        
        if file_prefix==None:
            print "Warning: file_prefix is not provided in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: file_prefix is not provided in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if whole_gene_list==None:
            print "Warning: whole_gene_list is not provided in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: whole_gene_list is not provided in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if tophat_genome_ref==None:
            print "Warning: tophat_genome_ref is not provided in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: tophat_genome_ref is not provided in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
                
        if END_N==None:
            print "Warning: End info is not provided in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: END_N is not provided in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
                
        if QF_path==None:
            print "Warning: QueryFuse_path info is not provided in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: QF_path is not provided in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)   
        
        if END_N=="first":
            END_O="second"
        elif END_N=="second":
            END_O="first"    
        else:
            print "Warning: End info is provided but not correct (neither first nor second) in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: END_N is provided but not correct (neither first nor second) in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if query_split==None:
            query_split=file_prefix+"unmapped.bam_"+END_N+"_mate_on_query.psl_split_ID_subtract_ID.txt"

	if single_to_both==None:
            single_to_both=file_prefix+"singleton_sorted_to_both_query_and_other_ID.txt"
	        	
        if unmap_sub_bed==None:
            unmap_sub_bed=file_prefix+"unmapped.bam_"+END_N+"_mate_on_query.psl_split_ID_subtract.bed"
	
        if unmap_query_psl==None:
            unmap_query_psl=file_prefix+"unmapped.bam_"+END_N+"_mate_on_query.psl"
	
        if mate_fa==None:
            mate_fa=file_prefix+"unmapped.bam_"+END_N+"_mate.fa"
        
        if single_query_bed==None:
            single_query_bed=file_prefix+"singleton_sorted_to_query.bed"
                
                
	#check whether the files provide is there.
	if not os.path.exists(whole_gene_list):
	    print "Warning: whole_gene_list:"+whole_gene_list+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: whole_gene_list:"+whole_gene_list+" is not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
	
        if not os.path.exists(tophat_genome_ref+".1.bt2"):
	    print "Warning: tophat_genome_ref:"+tophat_genome_ref+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: tophat_genome_ref:"+tophat_genome_ref+" is not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
	if not os.path.exists(file_prefix):
	    print "Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
                
        if not os.path.exists(query_split):
	    print "Warning: query_split:"+query_split+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: query_split ("+query_split+") not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(single_to_both):
            print "Warning: single_to_both:"+single_to_both+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: single_to_both("+single_to_both+") not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
                
        if not os.path.exists(unmap_sub_bed):
	    print "Warning: unmap_sub_bed:"+unmap_sub_bed+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: unmap_sub_bed("+unmap_sub_bed+") not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(unmap_query_psl):
	    print "Warning: unmap_query_psl:"+unmap_query_psl+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: unmap_query_psl("+unmap_query_psl+") not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(single_query_bed):
	    print "Warning: single_query_bed:"+single_query_bed+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: single_query_bed("+single_query_bed+") not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(mate_fa):
	    print "Warning: mate_fa:"+mate_fa+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: mate_fa ("+mate_fa+") is not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        
        #if not os.path.exists(out_folder):
	    #mkdir the folder
        #   os.makedirs(out_folder)
                
        if not os.path.exists(QF_path):
	    print "Warning: QueryFuse_path:"+QF_path+" is not found in QF_single_end_split_mate_both.py, exit."
            log_error.write("Warning: QueryFuse_path:"+QF_path+" is not found in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
	
        
        #Build the name for files
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID=query_split+"_split_mate_to_both_ID.txt"
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT=query_split+"_split_mate_to_both_ID_subtract.bed"
        UNMAP_MATE_FAI=mate_fa+".fai"
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FA=query_split+"_split_mate_to_both_ID_subtract.fa"
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FASTQ=query_split+"_split_mate_to_both_ID_subtract.fastq"
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM=query_split+"_split_mate_to_both_ID_subtract.sam"
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_BAM=query_split+"_split_mate_to_both_ID_subtract.bam"
        SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO=query_split+"_split_mate_to_both_ID_subtract_gene_anno.bed"


        #singleton(6) scenario
        step_name="get the split candidates which their mates are to both on "+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="Get mate_to_query_ID_subtract bed file of "+END_N+" end in QF_single_end_split_mate_both.py"
        singleton_6_cmd="Rscript "+QF_path+"/intersect_gene_list.R file.list1="+query_split+" file.list2="+single_to_both+" file.out="+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID
        if os.stat(query_split).st_size == 0:
            log_whole.write(query_split+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(query_split+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        elif os.stat(single_to_both).st_size == 0:
            log_whole.write(single_to_both+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(single_to_both+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
            singleton_6_cmd_status=QF_all_modules.resume_func(singleton_6_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=singleton_6_cmd_status[0]
	    QF_all_modules.key_step_check(singleton_6_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        
        #Use these IDs, get the subtract bed
        step_name="Get mate_to_query_ID_subtract bed file of "+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="get the subtract.fa of "+END_N+" end in QF_single_end_split_mate_both.py"
        Get_sub_bed_cmd=QF_path+"/grep_reads_from_Bed_by_ID.sh "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID+" "+unmap_sub_bed+" "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT+" Rscript "+QF_path+" "+LOG_ERR
        if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        elif os.stat(unmap_sub_bed).st_size == 0:
            log_whole.write(unmap_sub_bed+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(unmap_sub_bed+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    Get_sub_bed_cmd_status=QF_all_modules.resume_func(Get_sub_bed_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=Get_sub_bed_cmd_status[0]
	    QF_all_modules.key_step_check(Get_sub_bed_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
        
        #get the subtract.fa
        step_name="get the subtract.fa of "+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="convert .fa to fastq of"+END_N+" end in QF_single_end_split_mate_both.py"
        
	Get_sub_fa_cmd="fastaFromBed -fi "+mate_fa+" -bed "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT+" -fo "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FA+" -name"
            
	if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    Get_sub_fa_cmd_status=QF_all_modules.resume_func(Get_sub_fa_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=Get_sub_fa_cmd_status[0]
	    QF_all_modules.key_step_check(Get_sub_fa_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
              
        #convert .fa to fastq
        step_name="convert .fa to fastq of"+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="Bowtie2 with single end option of"+END_N+" end in QF_single_end_split_mate_both.py"
        convert_fa_fq_cmd="perl "+QF_path+"/fasta_to_fastq.pl "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FA+" > "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FASTQ
        if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FA).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FA+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FA+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    convert_fa_fq_cmd_status=QF_all_modules.resume_func(convert_fa_fq_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=convert_fa_fq_cmd_status[0]
	    QF_all_modules.key_step_check(convert_fa_fq_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
           
        #Bowtie2 single end option to the genome with most sensitive option
        step_name="Bowtie2 with single end option of"+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="transform sam into bam of"+END_N+" end in QF_single_end_split_mate_both.py"
        
	bowtie2_cmd="bowtie2 --very-sensitive-local -x "+tophat_genome_ref+" -U "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FASTQ+" > "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM
        
	if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FASTQ).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FASTQ+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_FASTQ+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    bowtie2_cmd_status=QF_all_modules.resume_func(bowtie2_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=bowtie2_cmd_status[0]
	    QF_all_modules.key_step_check(bowtie2_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
                
        #transform sam into bam
        step_name="transform sam into bam of"+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="intersectBam to see where the mates should align of"+END_N+" end in QF_single_end_split_mate_both.py"
        
	sam_to_bam_cmd="samtools view -S "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM+" -b > "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_BAM
        
	if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    sam_to_bam_cmd_status=QF_all_modules.resume_func(sam_to_bam_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=sam_to_bam_cmd_status[0]
	    QF_all_modules.key_step_check(sam_to_bam_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
           
        #intersectBam to see where these mates should align.
        step_name="intersectBam to see where the mates should align of"+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="merge into pairs of"+END_N+" end in QF_single_end_split_mate_both.py"
        intersectBam_cmd="intersectBed -abam "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_BAM+" -b "+whole_gene_list+" -bed -wo > "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO
        
	if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_BAM).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_BAM+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_BAM+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        elif os.stat(whole_gene_list).st_size == 0:
            log_whole.write(whole_gene_list+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(whole_gene_list+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    intersectBam_cmd_status=QF_all_modules.resume_func(intersectBam_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=intersectBam_cmd_status[0]
	    QF_all_modules.key_step_check(intersectBam_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
        
        #remove sam to save space
        if os.path.exists(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM):
            os.remove(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_SAM)
        
        #summary of the result by reading bed file
        #function is to figure out which has highest match, at the same time, figure out which end is the breakpoint by grouping (no need at all!!!)
        #GSFQ_single_end_split_mate_query_summary.sh $SINGLE_ON_QUERY_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO $8 $9

        #group into pairs and also filter our pair with too long or too much mismatch
        step_name="merge into pairs of"+END_N+" end in QF_single_end_split_mate_both.py"
        log_whole.write(step_name+'\n')
        next_step_name="This is the last step in this script"
        merger_cmd=QF_path+"/QF_single_end_split_mate_merger.sh "+SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO+" "+single_query_bed+" "+unmap_query_psl+" "+str(read_len)+" "+QF_path+" "+str(Align_percent)
        if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        elif os.stat(single_query_bed).st_size == 0:
            log_whole.write(single_query_bed+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(single_query_bed+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        elif os.stat(unmap_query_psl).st_size == 0:
            log_whole.write(unmap_query_psl+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n")
            log_error.write(unmap_query_psl+" is empty "+END_N+" end in QF_single_end_split_mate_both.py, exit.\n"); sys.exit(1)
        else:
	    merger_cmd_status=QF_all_modules.resume_func(merger_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=merger_cmd_status[0]
	    QF_all_modules.key_step_check(merger_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
        	
     


