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
This script is a Gene_specific_fusion_query subfunction: single_end pocessing. Specially for (5) split mate to other gene candidate
This script only process one end at a time. To process both, need to run for each end.

This script will finally generate an ordered and merged file that can be used to summarize and count spliting reads which the single aligned end is on query. 
=============================

Usage: python QF_single_end_split_mate_other.py
-h help

-i file_prefix that all the defualt can use                                                                         *[No default value]

-E end_info_in_English                                                                                              *[No default value, but should be first or second]

-q query_split_ID.txt                                                                                               [default value: file_prefixfolder, unmapped.bam_END_N_mate_on_query.psl_split_ID_subtract_ID.txt]

-Q to_other_only_ID(from the other mate) - location of the to_other_only_ID.txt(should use from the other mate)	    [default value: file_prefixfolder, singleton_sorted_to_other_only_ID_END_O_mate.txt]

-u Unmap_over_query_split_subtract_BED- Bed file which contain subtract information for unmap over query of the mate same as the query_split_ID         [default value file_prefixfolder, unmapped.bam_END_N_mate_on_query.psl_split_ID_subtract.bed]

-m mate_fa - fa file of the same mate as the query_split_ID                                                         [default value file_prefixfolder, unmapped.bam_END_N_mate.fa]

-U UNMAP_ON_QUERY_PSL - location of UNMAP_ON_QUERY_PSL                                                              [default value file_prefixfolder, unmapped.bam_END_N_mate_on_query.psl]

-w whole_gene_list.bed				                                                                    *[No default value]

-t tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                    *[No default value]

-l Read_length - length of reads		                                                                    [default value 99]

-s SINGLE_ON_OTHER_BED		                                                                                    [default value file_prefixfolder, singleton_sorted_to_other.bed]

-g LOG_folder                                                                                                       *[No default value]

-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]

-F QueryFuse_path                                                                                                   *[No default value]

-O size_other - the value of blat -stepSize option value for blat to other					    [default value is 11]

-a Align_percent: min percentage of alignment                                                                       [default value is 98]

-m minscore parameter for blat (which will define the min alignment length allowed)				    [default value is 11]

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
	tophat_genome_fa=None
        LOG_F=None
        END_N=None
	QF_path=None
        query_split=None
	single_to_other=None
	unmap_sub_bed=None
	unmap_query_psl=None
	read_len=99
	single_other_bed=None
        mate_fa=None
        resume_stat=0
        Align_percent=98
	size_other=11
	MIN_SCORE="11"
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:w:g:t:E:F:q:Q:u:o:U:s:m:l:r:R:P:Y:a:b:O:S:M:z')
        for opt in optlist:
            if opt[0] == '-h':
                print __doc__; sys.exit(2)
            elif opt[0] == '-w': whole_gene_list = opt[1]
            elif opt[0] == '-i': file_prefix = opt[1]
            elif opt[0] == '-g': LOG_F = opt[1]
            elif opt[0] == '-t': tophat_genome_fa =opt[1]
            elif opt[0] == '-E': END_N = opt[1]
            elif opt[0] == '-F': QF_path = opt[1]
            elif opt[0] == '-q': query_split =opt[1]
            elif opt[0] == '-Q': single_to_other= opt[1]
            elif opt[0] == '-u': unmap_sub_bed = opt[1]
            elif opt[0] == '-U': unmap_query_psl =opt[1]
	    elif opt[0] == '-s': single_other_bed = opt[1]
            elif opt[0] == '-m': mate_fa =opt[1]
            elif opt[0] == '-l': read_len = int(opt[1])
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-a': Align_percent =opt[1]
            elif opt[0] == '-O': size_other =int(opt[1])
	    elif opt[0] == '-M': MIN_SCORE = opt[1]
	
        if LOG_F==None:
             log_whole.write("Warning: LOG_F is not provided in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(2)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
        LOG_WHO=LOG_F+"/whole.log"
        log_error=open(LOG_ERR,"a")
        log_whole=open(LOG_WHO,"a")
        
        if file_prefix==None:
            print "Warning: file_prefix is not provided in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: file_prefix is not provided in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if tophat_genome_fa==None:
            print "Warning: tophat_genome_fa is not provided in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: tophat_genome_fa is not provided in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
                
        if END_N==None:
            print "Warning: End info is not provided in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: END_N is not provided in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
                
        if QF_path==None:
            print "Warning: QueryFuse_path info is not provided in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: QF_path is not provided in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)   
        
	
        if END_N=="first":
            END_O="second"
        elif END_N=="second":
            END_O="first"    
        else:
            print "Warning: End info is provided but not correct (neither first nor second) in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: END_N is provided but not correct (neither first nor second) in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if query_split==None:
            query_split=file_prefix+"unmapped.bam_"+END_N+"_mate_on_query.psl_split_ID_subtract_ID.txt"

	if single_to_other==None:
            single_to_other=file_prefix+"singleton_sorted_to_other_only_ID_"+END_O+"_mate.txt"
	        	
        if unmap_sub_bed==None:
            unmap_sub_bed=file_prefix+"unmapped.bam_"+END_N+"_mate_on_query.psl_split_ID_subtract.bed"
	
        if unmap_query_psl==None:
            unmap_query_psl=file_prefix+"unmapped.bam_"+END_N+"_mate_on_query.psl"
	
        if mate_fa==None:
            mate_fa=file_prefix+"unmapped.bam_"+END_N+"_mate.fa"
        
        if single_other_bed==None:
            single_other_bed=file_prefix+"singleton_sorted_to_other.bed"
	        
	#check whether the files provide is there.
        if not os.path.exists(tophat_genome_fa):
	    print "Warning: tophat_genome_fa:"+tophat_genome_fa+" is not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: tophat_genome_fa:"+tophat_genome_fa+" is not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
                
	if not os.path.exists(file_prefix):
	    print "Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(query_split):
	    print "Warning: query_split:"+query_split+" is not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: query_split ("+query_split+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(single_to_other):
            print "Warning: single_to_other not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: single_to_other("+single_to_other+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
                
        if not os.path.exists(unmap_sub_bed):
	    print "Warning: unmap_sub_bed not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: unmap_sub_bed("+unmap_sub_bed+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(unmap_query_psl):
	    print "Warning: unmap_query_psl not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: unmap_query_psl("+unmap_query_psl+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(single_other_bed):
	    print "Warning: single_other_bed not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: single_other_bed("+single_other_bed+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(mate_fa):
	    print "Warning: mate_fa not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: mate_fa("+mate_fa+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        #if not os.path.exists(out_folder):
	    #mkdir the folder
        #   os.makedirs(out_folder)
                
        if not os.path.exists(QF_path):
	    print "Warning: QueryFuse_path("+QF_path+") not found in QF_single_end_split_mate_other.py, exit."
            log_error.write("Warning: QueryFuse_path("+QF_path+") not found in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        
        #Build the name for files
        SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID=query_split+"_split_mate_to_other_ID.txt"	
        SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL=query_split+"_split_mate_to_other_blat.psl"
        
        #singleton(5) scenario
        #get the split candidates which their mates are not query
        step_name="get the split candidates which their mates are not query on "+END_N+" end in QF_single_end_split_mate_other.py"
        log_whole.write(step_name+'\n')
        next_step_name="Blat the subtract to their mates of "+END_N+" end in QF_single_end_split_mate_other.py"
        singleton_5_cmd="Rscript "+QF_path+"/intersect_gene_list.R file.list1="+query_split+" file.list2="+single_to_other+" file.out="+SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID
        if os.stat(query_split).st_size == 0:
            log_whole.write(query_split+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(query_split+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        elif os.stat(single_to_other).st_size == 0:
            log_whole.write(single_to_other+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(single_to_other+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        else:
            singleton_5_cmd_status=QF_all_modules.resume_func(singleton_5_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=singleton_5_cmd_status[0]
	    QF_all_modules.key_step_check(singleton_5_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        
        #Blat the subtract to their mates 
        step_name="Blat the subtract to their mates of "+END_N+" end in QF_single_end_split_mate_other.py"
        log_whole.write(step_name+'\n')
        next_step_name="pre process the summary of "+END_N+" end in QF_single_end_split_mate_other.py"
        Blat_sub_cmd=QF_path+"/blat_to_mate_no_grouping.sh "+SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID+" "+single_other_bed+" "+whole_gene_list+" "+tophat_genome_fa+" "+unmap_sub_bed+" "+mate_fa+" "+SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL+" "+QF_path+" "+LOG_ERR+" "+str(size_other)+" "+Align_percent+" "+str(read_len)+" "+MIN_SCORE
        #Blat_sub_cmd="python "+QF_path+"/blat_to_mate_no_grouping.py -r "+SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID+" -M "+single_other_bed+" -B "+unmap_sub_bed+" -F "+tophat_genome_fa+" -f "+mate_fa+" -o "+SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL+" -Q "+QF_path+" -g "+LOG_ERR+" -s "+str(size_other)+" -a "+Align_percent+" -c FALSE -R "+str(read_len)
        if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_OTHER_ID+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        elif os.stat(single_other_bed).st_size == 0:
            log_whole.write(single_other_bed+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(single_other_bed+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
	elif os.stat(unmap_sub_bed).st_size == 0:
            log_whole.write(unmap_sub_bed+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(unmap_sub_bed+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        else:
	    Blat_sub_cmd_status=QF_all_modules.resume_func(Blat_sub_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=Blat_sub_cmd_status[0]
	    QF_all_modules.key_step_check(Blat_sub_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
        
	#pre process the summary, because the final summary need the information of both ends.
        step_name="pre process the summary of "+END_N+" end in QF_single_end_split_mate_other.py"
        log_whole.write(step_name+'\n')
        next_step_name="This is the last step in this script"
        summary_cmd=QF_path+"/QF_single_end_split_mate_other_summary.sh "+SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL+" "+single_other_bed+" "+unmap_query_psl+" "+str(read_len)+" "+QF_path+" "+LOG_ERR+" "+str(Align_percent)
        if os.stat(SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL).st_size == 0:
            log_whole.write(SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
        elif os.stat(single_other_bed).st_size == 0:
            log_whole.write(single_other_bed+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(single_other_bed+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
	elif os.stat(unmap_query_psl).st_size == 0:
            log_whole.write(unmap_query_psl+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n")
            log_error.write(unmap_query_psl+" is empty "+END_N+" end in QF_single_end_split_mate_other.py, exit.\n"); sys.exit(1)
	else:
	    summary_cmd_status=QF_all_modules.resume_func(summary_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=summary_cmd_status[0]
	    QF_all_modules.key_step_check(summary_cmd_status, step_name, log_whole, log_error)
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n') 
