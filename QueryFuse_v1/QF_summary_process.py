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
This script is a QueryFuse subfunction: summary process after finishing pair and single/unmapped process.

First, it will group the annotated bed from all steps and seperate to span and split groups
Second, find the fusion breakpoint from splitting reads.
Third, extract spanning reads for these breakpoints and final report them.
Fourth, use the rest of the spanning reads to de novo grouping and give range of potential fusion.
Fifth, filter and score the result and give the final output.

=============================
Usage: python QF_summary_process.py
-h help

-i file_prefix that all the defualt can use (for intermedia files)                                                  *[No default value]

-B BAM_file_folder_prefix that has all three needed input bam files                                                 *[No default value]

-o result_folder_prefix that all the final important result should go                                               *[No default value]

-w whole_gene_list.bed				                                                                    *[No default value]

##-t tophat_genome_reference index folder - the folder path has indexed genome_reference.		                    *[No default value]

-T tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                    *[No default value]

-g LOG_folder                                                                                                       *[No default value]

-F QueryFuse_path                                                                                                   *[No default value]

-l Read_length - length of reads		                                                                    [default value 99]

-r resume_status: check whether user want to skip finished step or start over                                       [default value 0, not resume]

-a Align_percent: min percentage of alignment                                                                       [default value is 98]

-s Standard deviation of fragment size                                                                              [default value is 100]

-t number of split_reads											    [default value 1]

-a number of span_reads												    [default value 1]

-u number of sum of supporting reads										    [default value 2]

-m minscore parameter for blat (which will define the min alignment length allowed)				    [default value is 11]

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


def rc(dna):
    import string
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq


if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re, subprocess, time 
	import os
	import math
	from itertools import ifilter
	import QF_all_modules
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 8:
		print __doc__
		sys.exit(3)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	file_prefix=None
	whole_gene_list=None
	##tophat_genome_ref=None
        LOG_F=None
        bam_fd=None
	QF_path=None
        outresult_fd=None
        genome_fa=None
        query_bed=None
        
        read_len=99
	resume_stat=0
        Align_percent=98
        read_std="100"
	split_n="1"
	span_n="1"
	sum_n="2"
	filter_num="2"
	size_query="5"
	size_other="11"
	MIN_SCORE="11"
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:B:o:w:g:t:T:F:f:l:r:R:P:Y:a:b:Q:S:s:L:N:u:f:U:O:m:z')
        for opt in optlist:
            if opt[0] == '-h':
                print __doc__; sys.exit(2)
            elif opt[0] == '-i': file_prefix = opt[1]
            elif opt[0] == '-B': bam_fd = opt[1]
            elif opt[0] == '-o': outresult_fd = opt[1]
            elif opt[0] == '-w': whole_gene_list = opt[1]
            elif opt[0] == '-g': LOG_F = opt[1]
            ##elif opt[0] == '-t': tophat_genome_ref =opt[1]
            elif opt[0] == '-T': genome_fa = opt[1]
            elif opt[0] == '-F': QF_path = opt[1]
            elif opt[0] == '-l': read_len = int(opt[1])
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-a': Align_percent =int(opt[1])
            elif opt[0] == '-Q': query_bed =opt[1]
            elif opt[0] == '-s': read_std =opt[1]
            elif opt[0] == '-L': split_n = opt[1]
	    elif opt[0] == '-N': span_n = opt[1]
	    elif opt[0] == '-u': sum_n =opt[1]
	    elif opt[0] == '-f': filter_num = opt[1]
	    elif opt[0] == '-U': size_query =int(opt[1])
	    elif opt[0] == '-O': size_other =int(opt[1])
	    elif opt[0] == '-m': MIN_SCORE = opt[1]  
   
        if LOG_F==None:
            print "Warning: LOG_F is not provided in QF_summary_process.py, exit."; sys.exit(2)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
        LOG_WHO=LOG_F+"/whole.log"
        log_error=open(LOG_ERR,"a")
        log_whole=open(LOG_WHO,"a")
	
        #for parameter input needed.
        if file_prefix==None:
            print "Warning: file_prefix is not provided in QF_summary_process.py, exit."
            log_error.write("file_prefix is not provided in QF_summary_process.py, exit.\n"); sys.exit(1)
        
        if bam_fd==None:
            print "bam_folder is not provided"
            log_error.write("bam_folder is not provided\n"); sys.exit(1)
        
        if outresult_fd==None:
            print "Warning: outresult_folder is not provided in QF_summary_process.py, exit."
            log_error.write("outresult_folder is not provided in QF_summary_process.py, exit.\n"); sys.exit(1)
                
        if whole_gene_list==None:
            print "whole_gene_list is not provided"
            log_error.write("whole_gene_list is not provided\n"); sys.exit(1)
        
        #if tophat_genome_ref==None:
        #    print "tophat_genome_ref is not provided"
        #    log_error.write("tophat_genome_ref is not provided\n"); sys.exit(1)
        #        
        if genome_fa==None:
            print "genome_fa path is not provided"
            log_error.write("genome_fa path is not provided\n"); sys.exit(1)
                
        if QF_path==None:
            print "Warning: QueryFuse_path info is not provided in QF_summary_process.py, exit."
            log_error.write("QF_path is not provided in QF_summary_process.py, exit.\n"); sys.exit(1)   
        
        #for parameter with default
        if query_bed==None:
            query_bed=outresult_fd+"query_gene.bed"
	
	QUERY_FA=query_bed[:-3]+"fa"
                
	# whether the files provide is there.
	if not os.path.exists(whole_gene_list):
	    print "whole_gene_list not found"
            log_error.write("whole_gene_list not found\n"); sys.exit(1)
#	
        if not os.path.exists(genome_fa):
	    print "genome_fa not found"
            log_error.write("genome_fa not found\n"); sys.exit(1)
        
#        if not os.path.exists(tophat_genome_ref+".1.bt2"):
#	    print "tophat_genome_ref not found"
#            log_error.write("tophat_genome_ref not found\n"); sys.exit(1)
#        
	if not os.path.exists(file_prefix):
	    print "Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_summary_process.py, exit."
            log_error.write("intermediate file_prefix not found in QF_summary_process.py, exit.\n"); sys.exit(1)
        
#        if not os.path.exists(bam_fd):
#	    print "folder of bams not found"
#            log_error.write("folder of bams not found\n"); sys.exit(1)
#        
        if not os.path.exists(outresult_fd):
	    print "Warning: output result folder:"+outresult_fd+" is not found in QF_summary_process.py, exit."
            log_error.write("output result folder not found in QF_summary_process.py, exit.\n"); sys.exit(1)
                 

        #if not os.path.exists(out_folder):
	    #mkdir the folder
        #   os.makedirs(out_folder)
                
        if not os.path.exists(QF_path):
	    print "Warning: QueryFuse_path:"+QF_path+" is not found in QF_summary_process.py, exit."
            log_error.write("QueryFuse_path not found in QF_summary_process.py, exit.\n"); sys.exit(1)

        resume_stat_loc=resume_stat
	
        #hardcoded input file names for preprocess
        PAIR_SORT=file_prefix+"paired_sorted"
        

        step_name="split read summary in QF_summary_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="merge reported events within a range together in QF_summary_process.py"
	finish_summary_cmd="python "+QF_path+"/QF_whole_summary.py -i "+file_prefix+" -o "+outresult_fd+" -Q "+query_bed+" -a "+str(Align_percent)+" -l "+str(read_len)+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -f "+filter_num
        finish_summary_cmd_status=QF_all_modules.resume_func(finish_summary_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=finish_summary_cmd_status[0]
        QF_all_modules.key_step_check(finish_summary_cmd_status, step_name, log_whole, log_error)
	print "finished spliting read summary in QF_summary_process.py"
        log_whole.write("finished spliting read summary in QF_summary_process.py"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')
	
	#to generate splitting fa first for the future subgrouping.
	FUSION_BREAK_POINT_SUM_COUNT=outresult_fd+"/fusion_break_point_summary_count.txt"
	FUSION_BREAK_POINT_SUM_READ=outresult_fd+"/fusion_break_point_summary_reads.txt"
	FUSION_BREAK_POINT_SUM_COUNT_MERGE=outresult_fd+"/fusion_break_point_summary_count_merged.txt"
	FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT=outresult_fd+"/fusion_break_point_summary_count_merged_ref_filter.txt"
	FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT_ADJ=outresult_fd+"/fusion_break_point_summary_count_merged_ref_filter_adjust.txt"
	SUBGROUP_FOLDER=outresult_fd+"/fusion_supporting_sub_group/"
	SUBGRAPH_FOLDER=outresult_fd+"/fusion_supporting_graphs/"
	ALL_REF_TEMP_FA=SUBGROUP_FOLDER+"/all_ref_template.fa"
	step_name="merge reported events within a range together in QF_summary_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="Generate splitting fa first for the future subgrouping in QF_summary_process.py"
	merge_summary_cmd="python "+QF_path+"/fusion_counts_merging.py -c "+FUSION_BREAK_POINT_SUM_COUNT+" -r "+FUSION_BREAK_POINT_SUM_READ+" -o "+FUSION_BREAK_POINT_SUM_COUNT_MERGE+" -O "+SUBGROUP_FOLDER+" -q "+query_bed+" -g 5"
        merge_summary_cmd_status=QF_all_modules.resume_func(merge_summary_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=merge_summary_cmd_status[0]
        QF_all_modules.key_step_check(merge_summary_cmd_status, step_name, log_whole, log_error)
	print "finished merging reported events within a range together in QF_summary_process.py"
        log_whole.write("finished merging reported events within a range together in QF_summary_process.py"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')
	
	
	if not os.path.exists(SUBGROUP_FOLDER):
		os.mkdir(SUBGROUP_FOLDER)
	FUSION_BREAK_POINT_SUM_READ_ID=SUBGROUP_FOLDER+"/fusion_break_point_summary_reads_ID.txt"
	FUSION_BREAK_POINT_SUM_READ_FA1=SUBGROUP_FOLDER+"/fusion_break_point_summary_reads_1.fa"
	FUSION_BREAK_POINT_SUM_READ_FA2=SUBGROUP_FOLDER+"/fusion_break_point_summary_reads_2.fa"
	UNMAP_1_FA=bam_fd+"unmapped.bam_first_mate.fa"
	UNMAP_2_FA=bam_fd+"unmapped.bam_second_mate.fa"
	h_UNMAP_1_FA=open(UNMAP_1_FA,"r")
	h_UNMAP_2_FA=open(UNMAP_2_FA,"r")
	FA1_line1=h_UNMAP_1_FA.readline()
	FA2_line1=h_UNMAP_2_FA.readline()
	if FA1_line1[0]!=">":
		print "Warning: The input of unmapped.bam_first_mate.fa in not in correct format with > at the first row in QF_summary_process.py, exit.\n, exit."
		log_error.write("The input of unmapped.bam_first_mate.fa in not in correct format with > at the first row in QF_summary_process.py, exit.\n"); sys.exit(1)
	if FA2_line1[0]!=">":
		print "Warning: The input of unmapped.bam_second_mate.fa in not in correct format with > at the first row in QF_summary_process.py, exit.\n, exit."
		log_error.write("The input of unmapped.bam_second_mate.fa in not in correct format with > at the first row in QF_summary_process.py, exit.\n"); sys.exit(1)
	row_int_fa1=0
	row_int_fa2=0
	for line in h_UNMAP_1_FA:
		if line[0]==">":
			break
		else:
			row_int_fa1+=1
	for line in h_UNMAP_2_FA:
		if line[0]==">":
			break
		else:
			row_int_fa2+=1
	
	h_UNMAP_1_FA.close()
	h_UNMAP_2_FA.close()
	if row_int_fa1!=row_int_fa1:
		print "Warning: The input of unmapped.bam_first_mate.fa and unmapped.bam_second_mate.fa are not in the same format in QF_summary_process.py, exit.\n, exit."
		log_error.write("The input of The input of unmapped.bam_first_mate.fa and unmapped.bam_second_mate.fa are not in the same format in QF_summary_process.py, exit.\n"); sys.exit(1)
	


	step_name="Generate splitting fa first for the future subgrouping in QF_summary_process.py"
        log_whole.write(step_name+'\n')
	#log_whole.write("row_int_fa1"+str(row_int_fa1)+'\n')
        next_step_name="Breakpoint correction for merged events and shift range  in QF_summary_process.py"
	get_split_fa_cmd1="cut -f13 "+FUSION_BREAK_POINT_SUM_READ+" | cut -f1 -d'/' > "+FUSION_BREAK_POINT_SUM_READ_ID
	get_split_fa_cmd2="perl "+QF_path+"/extract_fa_by_readID.pl "+UNMAP_1_FA+" "+FUSION_BREAK_POINT_SUM_READ_ID+" "+str(row_int_fa1)+" "+FUSION_BREAK_POINT_SUM_READ_FA1
	get_split_fa_cmd3="perl "+QF_path+"/extract_fa_by_readID.pl "+UNMAP_2_FA+" "+FUSION_BREAK_POINT_SUM_READ_ID+" "+str(row_int_fa2)+" "+FUSION_BREAK_POINT_SUM_READ_FA2
	get_split_fa_cmd_all=get_split_fa_cmd1+"; "+get_split_fa_cmd2+"; "+get_split_fa_cmd3
	get_split_fa_cmd_all_status=QF_all_modules.resume_func(get_split_fa_cmd_all, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=get_split_fa_cmd_all_status[0]
        QF_all_modules.key_step_check(get_split_fa_cmd_all_status, step_name, log_whole, log_error)
	print "finished generating splitting fa first for the future subgrouping"
        log_whole.write("finished generating splitting fa first for the future subgrouping"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')
	
	#summary the result
	PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED=PAIR_SORT+"_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed"
	WHOLE_FUSION_SUM=outresult_fd+"/whole_fusion_sum_all.txt"
	WHOLE_FUSION_SUM_FILTERED=outresult_fd+"/whole_fusion_sum_filtered.txt"
	FUSION_SPLIT_SPAN_SUPPORT=outresult_fd+"/fusion_split_span_support.txt" 
	FUSION_SPAN_ONLY=outresult_fd+"/fusion_span_only.bed"
	FUSION_SPAN_ONLY_GROUP_PREFIX=outresult_fd+"/fusion_span_only_grouped_summary"
	
	if os.path.exists(FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT):
		os.remove(FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT)
	
	if os.path.exists(ALL_REF_TEMP_FA):
		os.remove(ALL_REF_TEMP_FA)
	
	
	#this should move to before doing spanning because, the realign location can affect the spanning support reads.
	#based on the filtered result, read line by line to generate template and supporting read graph for each fusion events.
	h_fusion_sum_filter=open(FUSION_BREAK_POINT_SUM_COUNT_MERGE,"r")
	for filter_fusion_line in h_fusion_sum_filter:
		#to get the correct format
		filter_fusion=filter_fusion_line.strip().split("\t")
		sub_read=SUBGROUP_FOLDER+"/"+"_".join(filter_fusion[:11])+"grouped_reads.txt"
		sub_read_ID=SUBGROUP_FOLDER+"/"+"_".join(filter_fusion[:11])+"grouped_reads_ID.txt"
		sub_read_fa1=SUBGROUP_FOLDER+"/"+"_".join(filter_fusion[:11])+"grouped_seqs.fa1"
		sub_read_fa2=SUBGROUP_FOLDER+"/"+"_".join(filter_fusion[:11])+"grouped_seqs.fa2"
		sub_read_fa=SUBGROUP_FOLDER+"/"+"_".join(filter_fusion[:11])+"grouped_seqs.fa"
		sub_read_file=SUBGROUP_FOLDER+"/"+"_".join(filter_fusion[:11])+"grouped_seqs.txt"
		get_split_fa_cmd1="cut -f13 "+sub_read+" | cut -f1 -d'/' > "+sub_read_ID
		get_split_fa_cmd2="perl "+QF_path+"/extract_fa_by_readID.pl "+FUSION_BREAK_POINT_SUM_READ_FA1+" "+sub_read_ID+" "+str(row_int_fa1)+" "+sub_read_fa1
		get_split_fa_cmd3="perl "+QF_path+"/extract_fa_by_readID.pl "+FUSION_BREAK_POINT_SUM_READ_FA2+" "+sub_read_ID+" "+str(row_int_fa2)+" "+sub_read_fa2
		#get_split_fa_cmd4="cat "+sub_read_fa1+" > "+sub_read_fa+"; cat "+sub_read_fa2+" >> "+sub_read_fa
		get_split_fa_cmd_all=get_split_fa_cmd1+"; "+get_split_fa_cmd2+"; "+get_split_fa_cmd3#+"; "+get_split_fa_cmd4
		subprocess.call(get_split_fa_cmd_all, shell = True)
		os.remove(sub_read_ID)
		
		#convert the fa into the format I want. into a function.
		h_sub_read_fa1=open(sub_read_fa1,"r")
		h_sub_read_file=open(sub_read_file,"w")
		data_sub_read_fa1=[line.strip() for line in h_sub_read_fa1]
		i_data_sub_read_fa1=0
		while i_data_sub_read_fa1 < len(data_sub_read_fa1):
			if data_sub_read_fa1[i_data_sub_read_fa1][0]==">":
				ID_line=data_sub_read_fa1[i_data_sub_read_fa1][1:]+"/1\t"
				seq_line=""
				for j in range(row_int_fa1):
					seq_line=seq_line+data_sub_read_fa1[i_data_sub_read_fa1+j+1]
				h_sub_read_file.write(ID_line+seq_line+"\n")
			i_data_sub_read_fa1+=row_int_fa1+1
			
		h_sub_read_fa1.close()
		#h_sub_read_file.close()
		
		h_sub_read_fa2=open(sub_read_fa2,"r")
		#h_sub_read_file=open(sub_read_file,"w")
		data_sub_read_fa2=[line.strip() for line in h_sub_read_fa2]
		i_data_sub_read_fa2=0
		while i_data_sub_read_fa2 < len(data_sub_read_fa2):
			if data_sub_read_fa2[i_data_sub_read_fa2][0]==">":
				ID_line=data_sub_read_fa2[i_data_sub_read_fa2][1:]+"/2\t"
				seq_line=""
				for j in range(row_int_fa1):
					seq_line=seq_line+data_sub_read_fa2[i_data_sub_read_fa2+j+1]
				h_sub_read_file.write(ID_line+seq_line+"\n")
			i_data_sub_read_fa2+=row_int_fa1+1
			
		h_sub_read_fa2.close()
		h_sub_read_file.close()
		
		os.remove(sub_read_fa1)
		os.remove(sub_read_fa2)
		
		
		
		#do supportig reads merge bp, adjustment and graphing.
		subgroup_adjust_cmd="python "+QF_path+"/supporting_reads_merge_bp_adjustment.py "+" -d "+sub_read_file+" -k "+"^".join(filter_fusion_line.strip().split("\t"))+" -o "+SUBGROUP_FOLDER+" -O "+SUBGRAPH_FOLDER+" -q "+FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT+" -Q "+sub_read+" -F "+QF_path+" -m "+MIN_SCORE
		subprocess.call(subgroup_adjust_cmd, shell = True)
	
	
	step_name="Breakpoint correction for merged events and shift range  in QF_summary_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="Group the spanning reads by using the existing breakpoints defined by splitting reads in QF_summary_process.py"
        bp_correct_cmd="python "+QF_path+"/breakpoint_adjustment_shift_range_scorening.py -o "+FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT_ADJ+" -O "+SUBGRAPH_FOLDER+" -w "+whole_gene_list+" -g "+LOG_ERR+" -t "+ALL_REF_TEMP_FA+" -T "+genome_fa+" -F "+QF_path+" -a "+str(Align_percent)+" -Q "+query_bed+" -q "+QUERY_FA+" -i "+str(size_query)+" -I "+str(size_other)+" -l "+str(read_len)+" -f "+file_prefix+" -m "+MIN_SCORE
	bp_correct_cmd_status=QF_all_modules.resume_func(bp_correct_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=bp_correct_cmd_status[0]
        QF_all_modules.key_step_check(bp_correct_cmd_status, step_name, log_whole, log_error)
	print "finished breakpoint correction for merged events and shift range ."
        log_whole.write("finished breakpoint correction for merged events and shift range ."+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')

	
	step_name="Group the spanning reads by using the existing breakpoints defined by splitting reads in QF_summary_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="Filter the summary with split,span,sum limit in QF_summary_process.py"
        group_span_from_split_cmd="python "+QF_path+"fusion_report_span_from_split.py "+" "+FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT_ADJ+" "+PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED+" "+str(read_len)+" "+read_std+" "+WHOLE_FUSION_SUM+" "+FUSION_SPLIT_SPAN_SUPPORT+" "+FUSION_SPAN_ONLY
        #group_span_from_split_cmd="python "+QF_path+"fusion_report_span_from_split.py "+" "+FUSION_BREAK_POINT_SUM_COUNT_MERGE+" "+PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_PAIRED_ANNO_BED+" "+str(read_len)+" "+read_std+" "+WHOLE_FUSION_SUM+" "+FUSION_SPLIT_SPAN_SUPPORT+" "+FUSION_SPAN_ONLY
	group_span_from_split_cmd_status=QF_all_modules.resume_func(group_span_from_split_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=group_span_from_split_cmd_status[0]
        QF_all_modules.optional_step_check(group_span_from_split_cmd_status, log_whole, log_error, "Grouping the spanning reads by using the existing breakpoints defined by splitting reads")
	print "finished grouping the spanning reads by using the existing breakpoints defined by splitting reads"
        log_whole.write("finished grouping the spanning reads by using the existing breakpoints defined by splitting reads"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')

	
	step_name="Filter the summary with split,span,sum limit in QF_summary_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="Group the span_only reads in QF_summary_process.py"
        summary_filter_cmd="Rscript "+QF_path+"whole_summary_score_ranking.R infile="+WHOLE_FUSION_SUM+" file.out="+WHOLE_FUSION_SUM_FILTERED+" split_n="+split_n+" span_n="+span_n+" sum_n="+sum_n
	#summary_filter_cmd="python "+QF_path+"whole_summary_filter.py "+" -i "+WHOLE_FUSION_SUM+" -o "+WHOLE_FUSION_SUM_FILTERED+" -t "+split_n+" -a "+span_n+" -u "+sum_n
        summary_filter_cmd_status=QF_all_modules.resume_func(summary_filter_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=summary_filter_cmd_status[0]
        QF_all_modules.optional_step_check(summary_filter_cmd_status, log_whole, log_error, "filtering the summary with split,span,sum limit")
	
		
	step_name="Group the span_only reads in QF_summary_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="This is the last step in this script"
        span_only_grouping_cmd="python "+QF_path+"QF_span_only_grouping.py "+" -i "+FUSION_SPAN_ONLY+" -o "+FUSION_SPAN_ONLY_GROUP_PREFIX+" -g "+LOG_F+" -l "+str(read_len)+" -f "+filter_num
        span_only_grouping_cmd_status=QF_all_modules.resume_func(span_only_grouping_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
	resume_stat_loc=span_only_grouping_cmd_status[0]
	QF_all_modules.optional_step_check(span_only_grouping_cmd_status, log_whole, log_error, "grouping the span_only reads in QF_summary_process.py")
	

	
