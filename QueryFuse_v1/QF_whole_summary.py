# -*- coding: cp936 -*-
"""
This script is the fusion result final summary function.
It does mainly three jobs now:
merge all the split reads together and filter artifact duplicate 
extract spanning support reads based on the splitting reads.
use the leftover spanning reads to de novo fusion pairs.

This script will finally generate an ordered and merged file that can be used to summarize and count spliting reads which the single aligned end is on query. 

=============================
Usage: python QF_whole_summary.py
-h help

-i file_prefix that all the defualt can use (for intermedia files)                                                  *[No default value]

##-B BAM_file_folder_prefix that has all three needed input bam files                                                 *[No default value]

-o result_folder_prefix that all the final important result should go                                               *[No default value]

##-w whole_gene_list.bed				                                                                    *[No default value]

-g LOG_folder                                                                                                       *[No default value]

-F QueryFuse_path                                                                                                   *[No default value]

-l Read_length - length of reads		                                                                    [default value 99]

-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]

-a Align_percent: min percentage of alignment                                                                       [default value is 98]

-Q Query_bed - file have information of the query gene                                                              [default value is outresult_fd/query_gene.bed]

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


#group and score by following
#for each gene build one file and save the supporting reads.
#sort 5' and 3' separately and see which one has fewer groups (unique ID), the one has few is correct. 
#column 5 is the score from alignment tool.
#For each gene group, remove duplicate rows before reporting.(if a read aligned to two different location in a gene, then keep it, highly impossible)
#For the mate on query, for each gene group, grep the reads from query_psl_file
#From the extracted psl. get the coordinate aligned on the query. also sort 5' and 3' spearately and see which one has fewer groups.






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
	if len(sys.argv) < 6:
		print __doc__
		sys.exit(3)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	file_prefix=None
	##whole_gene_list=None
        LOG_F=None
        ##bam_fd=None
	QF_path=None
        outresult_fd=None
        query_bed=None
	
        read_len=99
	resume_stat=0
        Align_percent=98
	filter_num="2"
     
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:B:o:w:g:F:l:r:R:P:Y:a:b:Q:f:z')
        for opt in optlist:
            if opt[0] == '-h':
                print __doc__; sys.exit(2)
            elif opt[0] == '-i': file_prefix = opt[1]
            ##elif opt[0] == '-B': bam_fd = opt[1]
            elif opt[0] == '-o': outresult_fd = opt[1]
            ##elif opt[0] == '-w': whole_gene_list = opt[1]
            elif opt[0] == '-g': LOG_F = opt[1]
            elif opt[0] == '-F': QF_path = opt[1]
            elif opt[0] == '-l': read_len = int(opt[1])
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-a': Align_percent =opt[1]
            elif opt[0] == '-Q': query_bed =opt[1]
	    elif opt[0] == '-f': filter_num = opt[1]
   
   
        if LOG_F==None:
            print "LOG_F is not provided in QF_whole_summary.py, exit."; sys.exit(2)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
        LOG_WHO=LOG_F+"/whole.log"
        log_error=open(LOG_ERR,"a")
        log_whole=open(LOG_WHO,"a")
        
        #for parameter input needed.
        if file_prefix==None:
            print "file_prefix is not provided in QF_whole_summary.py, exit."
            log_error.write("file_prefix is not provided in QF_whole_summary.py, exit.\n"); sys.exit(1)
        
        if outresult_fd==None:
            print "outresult_folder is not provided in QF_whole_summary.py, exit."
            log_error.write("outresult_folder is not provided in QF_whole_summary.py, exit.\n"); sys.exit(1)
                
        if QF_path==None:
            print "QueryFuse_path info is not provided in QF_whole_summary.py, exit."
            log_error.write("QF_path is not provided in QF_whole_summary.py, exit.\n"); sys.exit(1)   
        
        #for parameter with default  
        if query_bed==None:
            query_bed=outresult_fd+"query_gene.bed"
                
                	
        if not os.path.exists(file_prefix):
	    print "intermediate file_prefix not found in QF_whole_summary.py, exit."
            log_error.write("intermediate file_prefix not found in QF_whole_summary.py, exit.\n"); sys.exit(1)

        if not os.path.exists(outresult_fd):
	    print "output result folder not found in QF_whole_summary.py, exit."
            log_error.write("output result folder not found in QF_whole_summary.py, exit.\n"); sys.exit(1)
                
        if not os.path.exists(query_bed):
	    print "query_bed not found in QF_whole_summary.py, exit."
            log_error.write("query_bed ("+query_bed+") not found in QF_whole_summary.py, exit.\n"); sys.exit(1)

        if not os.path.exists(QF_path):
	    log_error.write(QF_path+"\n")
	    print "QueryFuse_path not found in QF_whole_summary.py, exit."
            log_error.write("QueryFuse_path not found in QF_whole_summary.py, exit.\n"); sys.exit(1)
        
        #hardcoded input file names
        SINGLE_ON_QUERY_1_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM=file_prefix+"unmapped.bam_first_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_query_ID_subtract_gene_anno.bed_sort_merge_with_mate_summary.bed"
        SINGLE_ON_QUERY_2_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM=file_prefix+"unmapped.bam_second_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_query_ID_subtract_gene_anno.bed_sort_merge_with_mate_summary.bed"
        SINGLE_ON_QUERY_1_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM=file_prefix+"unmapped.bam_first_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_other_blat.psl_summary_sorted.bed_merge_with_mate_summary.bed_transform.bed"
        SINGLE_ON_QUERY_2_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM=file_prefix+"unmapped.bam_second_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_other_blat.psl_summary_sorted.bed_merge_with_mate_summary.bed_transform.bed"
        SINGLE_ON_QUERY_1_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM=file_prefix+"unmapped.bam_first_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_both_ID_subtract_gene_anno.bed_sort_merge_with_mate_summary.bed"
        SINGLE_ON_QUERY_2_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM=file_prefix+"unmapped.bam_second_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_both_ID_subtract_gene_anno.bed_sort_merge_with_mate_summary.bed"
        FUSION_BREAK_POINT_SUM=outresult_fd+"/fusion_break_point_summary.txt"
        SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED=file_prefix+"/split_summary_merged.bed"
        SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUP=file_prefix+"/split_summary_merged_dup_removed.bed"
        SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI=file_prefix+"/split_summary_merged_PCR_duplicate_removed.bed"
        SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN=outresult_fd+"/split_summary_merged_cleaned.bed"
        SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO=outresult_fd+"/split_summary_merged_cleaned_bp_anno.bed"
        FUSION_BREAK_POINT_SUM_COUNT=outresult_fd+"/fusion_break_point_summary_count.txt"
        FUSION_BREAK_POINT_SUM_READ=outresult_fd+"/fusion_break_point_summary_reads.txt"
        
        #check the existence of files and merge them.
        if not os.path.exists(query_bed):
            log_whole.write(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is already exist, remove it.")
            os.remove(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED)
        
        step_name="check the existence of files and merge them. in QF_whole_summary.py"
        log_whole.write(step_name+'\n')
        check_split_query_1_cmd="cat "+SINGLE_ON_QUERY_1_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM+" >> "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED
        check_split_query_2_cmd="cat "+SINGLE_ON_QUERY_2_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM+" >> "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED
        check_split_other_1_cmd="cat "+SINGLE_ON_QUERY_1_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM+" >> "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED
        check_split_other_2_cmd="cat "+SINGLE_ON_QUERY_2_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM+" >> "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED
        check_split_both_1_cmd="cat "+SINGLE_ON_QUERY_1_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM+" >> "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED
        check_split_both_2_cmd="cat "+SINGLE_ON_QUERY_2_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM+" >> "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED
        if not os.path.exists(SINGLE_ON_QUERY_1_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM):
            log_whole.write("Note: "+SINGLE_ON_QUERY_1_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM+" is not exist"+'\n')
        else:
            if os.stat(SINGLE_ON_QUERY_1_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM).st_size == 0:
                log_whole.write("Note: "+SINGLE_ON_QUERY_1_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM+" is empty"+'\n')
            else:
                subprocess.call(check_split_query_1_cmd,shell=True)
        
        if not os.path.exists(SINGLE_ON_QUERY_2_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM):
            log_whole.write("Note: "+SINGLE_ON_QUERY_2_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM+" is not exist"+'\n')
        else:
            if os.stat(SINGLE_ON_QUERY_2_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM).st_size == 0:
                log_whole.write("Note: "+SINGLE_ON_QUERY_2_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO_SUM+" is empty"+'\n')
            else:
                subprocess.call(check_split_query_2_cmd,shell=True)
        
        if not os.path.exists(SINGLE_ON_QUERY_1_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM):
            log_whole.write("Note: "+SINGLE_ON_QUERY_1_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM+" is not exist"+'\n')
        else:
            if os.stat(SINGLE_ON_QUERY_1_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM).st_size == 0:
                log_whole.write("Note: "+SINGLE_ON_QUERY_1_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM+" is empty"+'\n')
            else:
                subprocess.call(check_split_other_1_cmd,shell=True)
        
        if not os.path.exists(SINGLE_ON_QUERY_2_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM):
            log_whole.write("Note: "+SINGLE_ON_QUERY_2_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM+" is not exist"+'\n')
        else:
            if os.stat(SINGLE_ON_QUERY_2_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM).st_size == 0:
                log_whole.write("Note: "+SINGLE_ON_QUERY_2_SPLIT_MATE_OTHER_ID_SUBTRACT_ANNO_SUM+" is empty"+'\n')
            else:
                subprocess.call(check_split_other_2_cmd,shell=True)
        
        if not os.path.exists(SINGLE_ON_QUERY_1_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM):
            log_whole.write("Note: "+SINGLE_ON_QUERY_1_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM+" is not exist"+'\n')
        else:
            if os.stat(SINGLE_ON_QUERY_1_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM).st_size == 0:
                log_whole.write("Note: "+SINGLE_ON_QUERY_1_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM+" is empty"+'\n')
            else:
                subprocess.call(check_split_both_1_cmd,shell=True)
        
        if not os.path.exists(SINGLE_ON_QUERY_2_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM):
            log_whole.write("Note: "+SINGLE_ON_QUERY_2_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM+" is not exist"+'\n')
        else:
            if os.stat(SINGLE_ON_QUERY_2_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM).st_size == 0:
                log_whole.write("Note: "+SINGLE_ON_QUERY_2_SPLIT_MATE_BOTH_ID_SUBTRACT_ANNO_SUM+" is empty"+'\n')
            else:
                subprocess.call(check_split_both_2_cmd,shell=True)
           
        step_name="filter the artifact duplicates in QF_whole_summary.py"
        log_whole.write(step_name+'\n')
        next_step_name="further filter PCR artifact in QF_whole_summary.py"
        sort_merge_cmd="sort -k2 "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" | uniq > "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUP
        if not os.path.exists(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED):
            #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is not exist, exit."
            log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is not exist, which means no splitting reads are found. do span only grouping before exiting."+'\n')
            log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is not exist, which means no splitting reads are found. do span only grouping before exiting."+'\n')
            
	    PAIR_SORT=file_prefix+"paired_sorted"
	    PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED=PAIR_SORT+"_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed"
	    FUSION_SPAN_ONLY_GROUP_PREFIX=outresult_fd+"/fusion_span_only_grouped_summary"
	    span_only_grouping_cmd="python "+QF_path+"QF_span_only_grouping.py "+" -i "+PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED+" -o "+FUSION_SPAN_ONLY_GROUP_PREFIX+" -g "+LOG_F+" -l "+str(read_len)+" -f "+filter_num
	    span_only_grouping_cmd_status=QF_all_modules.resume_func(span_only_grouping_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
	    resume_stat=span_only_grouping_cmd_status[0]
	    QF_all_modules.optional_step_check(span_only_grouping_cmd_status, log_whole, log_error, "grouping the span_only reads in QF_summary_process.py")
	    sys.exit(4)
        else:
            if os.stat(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED).st_size == 0:
                #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is empty, exit."
                log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is empty, which means no splitting reads are found. exit."+'\n')
                log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED+" is empty, which means no splitting reads are found. exit."+'\n')
                PAIR_SORT=file_prefix+"paired_sorted"
		PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED=PAIR_SORT+"_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed"
		FUSION_SPAN_ONLY_GROUP_PREFIX=outresult_fd+"/fusion_span_only_grouped_summary"
		span_only_grouping_cmd="python "+QF_path+"QF_span_only_grouping.py "+" -i "+PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED+" -o "+FUSION_SPAN_ONLY_GROUP_PREFIX+" -g "+LOG_F+" -l "+str(read_len)+" -f "+filter_num
		span_only_grouping_cmd_status=QF_all_modules.resume_func(span_only_grouping_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
		resume_stat=span_only_grouping_cmd_status[0]
		QF_all_modules.optional_step_check(span_only_grouping_cmd_status, log_whole, log_error, "grouping the span_only reads in QF_summary_process.py")
		sys.exit(4)
            else:    
                sort_merge_cmd_status=QF_all_modules.resume_func(sort_merge_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
                resume_stat=sort_merge_cmd_status[0]
                QF_all_modules.key_step_check(sort_merge_cmd_status, step_name, log_whole, log_error)
	
        
        step_name="further filter PCR artifact in QF_whole_summary.py"
        log_whole.write(step_name+'\n')
        next_step_name="filter wrong match reads in QF_whole_summary.py"
        filter_PCR_cmd="python "+QF_path+"/eliminate_dup_for_summary_merged.py "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUP+" "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI
        filter_PCR_cmd_status=QF_all_modules.resume_func(filter_PCR_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
        resume_stat=filter_PCR_cmd_status[0]
        QF_all_modules.key_step_check(filter_PCR_cmd_status, step_name, log_whole, log_error)
	
        step_name="filter wrong match reads in QF_whole_summary.py"
        log_whole.write(step_name+'\n')
        next_step_name="annotate the bp for each read in QF_whole_summary.py"
        filter_wm_cmd="python "+QF_path+"/eliminate_wrong_split_pair_match_for_summary_merged.py "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" "+str(read_len)+" "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" "+str(Align_percent)
        if not os.path.exists(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI):
            #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" is not exist, exit."
            log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" is not exist, which means no meaningful splitting reads are found. exit."+'\n')
            log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" is not exist, which means no meaningful splitting reads are found. exit."+'\n')
            sys.exit(4)
        else:
            if os.stat(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI).st_size == 0:
                #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" is empty, exit."
                log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" is empty, which means no meaningful splitting reads are found. exit."+'\n')
                log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_DUPLI+" is empty, which means no meaningful splitting reads are found. exit."+'\n')
                sys.exit(4)
            else:
                filter_wm_cmd_status=QF_all_modules.resume_func(filter_wm_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
                resume_stat=filter_wm_cmd_status[0]
                QF_all_modules.key_step_check(filter_wm_cmd_status, step_name, log_whole, log_error)
	       
        step_name="annotate the bp for each read in QF_whole_summary.py"
        log_whole.write(step_name+'\n')
        next_step_name="group the same breakpoints and count and output in QF_whole_summary.py"
        anno_cmd="python "+QF_path+"/breakpoint_asign_summary_merged.py "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" "+str(read_len)+" "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" "+str(Align_percent)
        if not os.path.exists(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN):
            #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" is not exist, exit."
            log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" is not exist, which means no corretly aligned splitting reads are found. exit."+'\n')
            log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" is not exist, which means no corretly aligned splitting reads are found. exit."+'\n')
            sys.exit(4)
        else:
            if os.stat(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN).st_size == 0:
                #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" is empty, exit."
                log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" is empty, which means no corretly aligned splitting reads are found. exit."+'\n')
                log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN+" is empty, which means no corretly aligned splitting reads are found. exit."+'\n')
                sys.exit(4)
            else:
                anno_cmd_status=QF_all_modules.resume_func(anno_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
                resume_stat=anno_cmd_status[0]
                QF_all_modules.key_step_check(anno_cmd_status, step_name, log_whole, log_error)
	
        
        step_name="group the same breakpoints and count and output in QF_whole_summary.py"
        log_whole.write(step_name+'\n')
        next_step_name="This is the last step in this script"
        group_bp_cmd="python "+QF_path+"/fusion_report_summary_merged.py "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" "+query_bed+" "+FUSION_BREAK_POINT_SUM_COUNT+" "+FUSION_BREAK_POINT_SUM_READ
        if not os.path.exists(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO):
            #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" is not exist, exit."
            log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" is not exist, which means splitting reads can not be merged into an annotation file. exit."+'\n')
            log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" is not exist, which means splitting reads can not be merged into an annotation file. exit."+'\n')
            sys.exit(4)
        else:
            if os.stat(SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO).st_size == 0:
                #print "Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" is empty, exit."
                log_whole.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" is empty, which means splitting reads can not be merged into an annotation file. exit."+'\n')
                log_error.write("Warning: "+SINGLE_ON_QUERY_SPLIT_SUMMERY_MERGED_FILTER_CLEAN_ANNO+" is empty, which means splitting reads can not be merged into an annotation file. exit."+'\n')
                sys.exit(4)
            else:
                group_bp_cmd_status=QF_all_modules.resume_func(group_bp_cmd, resume_stat, step_name, next_step_name, LOG_OUT)
                resume_stat=group_bp_cmd_status[0]
                QF_all_modules.key_step_check(group_bp_cmd_status, step_name, log_whole, log_error)
	
        
    