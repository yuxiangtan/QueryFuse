# -*- coding: cp936 -*-
#   Copyright {2015} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
This script is a QueryFuse subfunction: pair_end pocessing

First, it will get query related paired aligned reads.
Second, after classifying reads into sub-categories, run analysis on each senario.

All the output files are in file_prefix(intermediate) folder and will be made use of in the summary step.

=============================
Usage: python QF_pair_end_process.py
-h help

-i file_prefix that all the defualt can use (for intermedia files)                                                  *[No default value]

-B BAM_file_folder_prefix that has all three needed input bam files                                                 *[No default value]

-o result_folder_prefix that all the final important result should go                                               *[No default value]

-w whole_gene_list.bed				                                                                    *[No default value]

-t tophat_genome_reference index folder - the folder path has indexed genome_reference.		                    *[No default value]

-T tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                    *[No default value]

-g LOG_folder                                                                                                       *[No default value]

-F QueryFuse_path                                                                                                   *[No default value]

-l Read_length - length of reads		                                                                    [default value 99]

-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]

-a Align_percent: min percentage of alignment                                                                       [default value is 98]

-q size_query - the value of blat -stepSize option value for blat to query					    [default value is 5]


============================

Python & Module requirement:
Versions: 2.7 or above
Module: No additional Python Module is required.

============================
Library file requirement:
Not Standalone version, few library file is required.
bowtie2 is need to reinstall and preloaded.
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
	##Liye own common function,class loading
	#from Constant_Library import *
	#from General_Library import *
	#from File_Class import *
	#from Sequencing_Library import *
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 8:
		print __doc__
		sys.exit(3)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	file_prefix=None
	whole_gene_list=None
	tophat_genome_ref=None
        LOG_F=None
        bam_fd=None
	QF_path=None
        outresult_fd=None
        genome_fa=None
        query_bed=None
        
        read_len="99"
	resume_stat=0
        Align_percent="98"
	size_query=5
     
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:B:o:w:g:t:T:F:l:r:R:P:Y:a:b:q:Q:O:S:z')
        for opt in optlist:
            if opt[0] == '-h':
                print __doc__; sys.exit(2)
            elif opt[0] == '-i': file_prefix = opt[1]
            elif opt[0] == '-B': bam_fd = opt[1]
            elif opt[0] == '-o': outresult_fd = opt[1]
            elif opt[0] == '-w': whole_gene_list = opt[1]
            elif opt[0] == '-g': LOG_F = opt[1]
            elif opt[0] == '-t': tophat_genome_ref =opt[1]
            elif opt[0] == '-T': genome_fa = opt[1]
            elif opt[0] == '-F': QF_path = opt[1]
            elif opt[0] == '-l': read_len = opt[1]
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-a': Align_percent =opt[1]
            elif opt[0] == '-Q': query_bed =opt[1]
            elif opt[0] == '-q': size_query =int(opt[1])
	    #because I check in pair end, only blat to query is used.
	    #elif opt[0] == '-O': blat_other =opt[1]
            
   
   
        if LOG_F==None:
            print "Warning: LOG_F is not provided in QF_pair_end_process.py, exit."; sys.exit(2)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
	LOG_WHO=LOG_F+"/whole.log"
        
        log_error=open(LOG_ERR,"a")
	log_whole=open(LOG_WHO,"a")
        
        #for parameter input needed.
        if file_prefix==None:
            print "Warning: file_prefix is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: file_prefix is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
        if bam_fd==None:
            print "Warning: bam_folder is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: bam_folder is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
        if outresult_fd==None:
            print "Warning: outresult_folder is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: outresult_folder is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)
                
        if whole_gene_list==None:
            print "Warning: whole_gene_list is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: whole_gene_list is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
        if tophat_genome_ref==None:
            print "Warning: tophat_genome_ref is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: tophat_genome_ref is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)
                
        if genome_fa==None:
            print "Warning: genome_fa path is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: genome_fa path is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)
                
        if QF_path==None:
            print "Warning: QueryFuse_path info is not provided in QF_pair_end_process.py, exit."
            log_error.write("Warning: QF_path is not provided in QF_pair_end_process.py, exit.\n"); sys.exit(1)   
        
        #for parameter with default
        if query_bed==None:
            query_bed=outresult_fd+"query_gene.bed"
            
	#check whether the files provide is there.
	if not os.path.exists(whole_gene_list):
	    print "Warning: whole_gene_list:"+whole_gene_list+" is not found in QF_pair_end_process.py, exit."
            log_error.write("whole_gene_list not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)
	
        if not os.path.exists(genome_fa):
	    print "Warning: genome_fa"+genome_fa+" is not found in QF_pair_end_process.py, exit."
            log_error.write("genome_fa not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(tophat_genome_ref+".1.bt2"):
	    print "Warning: tophat_genome_ref:"+tophat_genome_ref+" is not found in QF_pair_end_process.py, exit."
            log_error.write("tophat_genome_ref not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
	if not os.path.exists(file_prefix):
	    print "Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_pair_end_process.py, exit."
            log_error.write("intermediate file_prefix not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(bam_fd):
	    print "Warning: folder of bams:"+bam_fd+" is not found in QF_pair_end_process.py, exit."
            log_error.write("folder of bams not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(outresult_fd):
	    print "Warning: output result folder:"+outresult_fd+" is not found in QF_pair_end_process.py, exit."
            log_error.write("output result folder not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)
                 

        #if not os.path.exists(out_folder):
	    #mkdir the folder
        #   os.makedirs(out_folder)
                
        if not os.path.exists(QF_path):
	    print "QueryFuse_path:"+QF_path+" is not found in QF_pair_end_process.py, exit."
            log_error.write("QueryFuse_path not found in QF_pair_end_process.py, exit.\n"); sys.exit(1)

        resume_stat_loc=resume_stat
	
        #hardcoded input file names for preprocess
        PAIR_BAM=bam_fd+"paired.bam"
        PAIR_SORT=file_prefix+"paired_sorted"
        PAIR_SORT_BAM=bam_fd+"paired_sorted.bam"
        PAIR_TO_QUERY_BAM=PAIR_SORT+"_to_query.bam"
	PAIR_TO_QUERY_SAM=PAIR_SORT+"_to_query.sam"
        PAIR_TO_QUERY_BED=PAIR_SORT+"_to_query.bed"
        PAIR_TO_QUERY_FILTER_BED=PAIR_SORT+"_to_query_filter.bed"
    
        step_name="Get proper pairs have at least one end fell into query region in QF_pair_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="BAM to Bed for pair reads aligned to query in QF_pair_end_process.py"
        pair_to_bam_cmd="pairToBed -abam "+PAIR_SORT_BAM+" -b "+query_bed+" > "+PAIR_TO_QUERY_BAM+" ;samtools view "+PAIR_TO_QUERY_BAM+" > "+PAIR_TO_QUERY_SAM
	#print pair_to_bam_cmd
        pair_to_bam_status=QF_all_modules.resume_func(pair_to_bam_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=pair_to_bam_status[0]
	#this is the key step, f this fail, all the rest will fail. As a result, set system exit here.
	QF_all_modules.key_step_check(pair_to_bam_status, step_name, log_whole, log_error)
	log_whole.write("finished getting proper pairs have at least one end fell into query region"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')


        step_name="BAM to Bed for pair reads aligned to query in QF_pair_end_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="Get the IDs that are not fully within the range of query gene from the previous list in QF_pair_end_process.py"
        bamtobed_cmd="bamToBed -i "+PAIR_TO_QUERY_BAM+" > "+PAIR_TO_QUERY_BED
	bamtobed_cmd_status=QF_all_modules.resume_func(bamtobed_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=bamtobed_cmd_status[0]
        #this is the key step, f this fail, all the rest will fail. As a result, set system exit here.
	QF_all_modules.key_step_check(bamtobed_cmd_status, step_name, log_whole, log_error)
        log_whole.write("finished transforming BAM to Bed for pair reads aligned to query"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')

        step_name="Get the IDs that are not fully within the range of query gene from the previous list in QF_pair_end_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="Process Pair-align scenario a in QF_pair_end_process.py"
        substractBed_cmd="subtractBed -a "+PAIR_TO_QUERY_BED+" -b "+query_bed+" > "+PAIR_TO_QUERY_FILTER_BED
	substractBed_cmd_status=QF_all_modules.resume_func(substractBed_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=substractBed_cmd_status[0]
	#this is the key step, f this fail, all the rest will fail. As a result, set system exit here.
	QF_all_modules.key_step_check(bamtobed_cmd_status, step_name, log_whole, log_error)
        log_whole.write("finished getting the IDs that are not fully within the range of query gene"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')


        step_name="Process Pair-align scenario a in QF_pair_end_process.py"
        log_whole.write(step_name+'\n')
        #next_step_name="Process Pair-align scenario b in QF_pair_end_process.py"
        next_step_name="This is the end for pair_end_process"
	pairA_cmd=QF_path+"QF_pair_A.sh "+file_prefix+" "+bam_fd+" "+outresult_fd+" "+whole_gene_list+" "+genome_fa+" "+read_len+" "+tophat_genome_ref+" "+LOG_ERR+" "+QF_path+" "+Align_percent+" "+str(size_query)
        pairA_cmd_status=QF_all_modules.resume_func(pairA_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=pairA_cmd_status[0]
        QF_all_modules.optional_step_check(pairA_cmd_status, log_whole, log_error, "processing Pair-align scenario a in QF_pair_end_process.py")
	

        
#        step_name="Process Pair-align scenario b in QF_pair_end_process.py"
#        log_whole.write(step_name+'\n')
#        next_step_name="This is the end for pair_end_process"
#        pairB_cmd=QF_path+"QF_pair_B.sh "+file_prefix+" "+bam_fd+" "+outresult_fd+" "+whole_gene_list+" "+genome_fa+" "+read_len+" "+tophat_genome_ref+" "+LOG_ERR+" "+QF_path+" "+Align_percent+" "+str(size_query)
#	pairB_cmd_status=QF_all_modules.resume_func(pairB_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
#	resume_stat_loc=pairB_cmd_status[0]
#	QF_all_modules.optional_step_check(pairB_cmd_status, log_whole, log_error, "processing Pair-align scenario b in QF_pair_end_process.py")
        
