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
This script is a Gene_specific_fusion_query subfunction: single_end pocessing

bam_folder
intermedia_folder=file_prefix
result_folder to put result in

tophat_ref
human_genome.fa should be different

3 bams as option

need to seperate the bam files with the intermediate files, as a result, can keep the names but not the real paths
also, because there are too many outputs, as a result, all use hardcoded name rather then optional names.


This script will finally generate an ordered and merged file that can be used to summarize and count spliting reads which the single aligned end is on query. 
=============================
Usage: python QF_single_end_process.py
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

-U size_query - the value of blat -stepSize option value for blat to query					    [default value is 5]

-O size_other - the value of blat -stepSize option value for blat to other					    [default value is 11]

-Q query_bed path												    *[No default value]

-q query_ID													    *[No default value]

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

#resume the run function
#def resume_func(cmd, resume_stat, step_name, next_step_name, LOG_OUT):
#    if resume_stat == 1:
#        h_LOG_OUT=open(LOG_OUT,"r")
#        count_res = 0
#        for line in h_LOG_OUT:
#	    if next_step_name in line:
#                count_res +=1
#	h_LOG_OUT.close()            
#        if count_res==0:
#            subprocess.call(cmd,shell=True)
#	    resume_stat=0
#        else:
#            print step_name+" pass"
#    else:
#        log_outf=open(LOG_OUT,"a")
#        log_outf.write(step_name+"\n")
#        subprocess.call(cmd,shell=True)
#        log_outf.close()
#    return resume_stat


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
	if len(sys.argv) < 11:
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
        query_ID=None
        
        read_len=99
	resume_stat=0
        Align_percent=98
        size_query=5
	size_other=11
     
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:B:o:w:g:t:T:F:l:r:R:P:Y:a:b:U:O:Q:q:S:z')
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
            elif opt[0] == '-l': read_len = int(opt[1])
            elif opt[0] == '-r': resume_stat = int(opt[1])
            elif opt[0] == '-a': Align_percent =int(opt[1])
            elif opt[0] == '-U': size_query =int(opt[1])
	    elif opt[0] == '-O': size_other =int(opt[1])
            elif opt[0] == '-Q': query_bed =opt[1]
            elif opt[0] == '-q': query_ID =opt[1]
            
   
   
        if LOG_F==None:
            print "Warning: LOG_F is not provided in QF_single_end_process.py, exit."; sys.exit(2)
        
        if not os.path.exists(LOG_F):
	    os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
        LOG_WHO=LOG_F+"/whole.log"
        log_error=open(LOG_ERR,"a")
        log_whole=open(LOG_WHO,"a")
	
        #for parameter input needed.
        if file_prefix==None:
            print "Warning: file_prefix is not provided in QF_single_end_process.py, exit."
            log_error.write("file_prefix is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
        if bam_fd==None:
            print "Warning: bam_folder is not provided in QF_single_end_process.py, exit."
            log_error.write("bam_folder is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
        if outresult_fd==None:
            print "Warning: outresult_folder is not provided in QF_single_end_process.py, exit."
            log_error.write("outresult_folder is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)
                
        if whole_gene_list==None:
            print "Warning: whole_gene_list is not provided in QF_single_end_process.py, exit."
            log_error.write("whole_gene_list is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
        if tophat_genome_ref==None:
            print "Warning: tophat_genome_ref is not provided in QF_single_end_process.py, exit."
            log_error.write("tophat_genome_ref is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)
                
        if genome_fa==None:
            print "Warning: genome_fa path is not provided in QF_single_end_process.py, exit."
            log_error.write("genome_fa path is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)
                
        if QF_path==None:
            print "Warning: QueryFuse_path info is not provided in QF_single_end_process.py, exit."
            log_error.write("QF_path is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)   
        
        if query_ID==None:
            print "Warning: The ID of query gene is not provided in QF_single_end_process.py, exit."
            log_error.write("The ID of query gene is not provided in QF_single_end_process.py, exit.\n"); sys.exit(1)   
        
        
        #for parameter with default
        if query_bed==None:
            query_bed=outresult_fd+"query_gene.bed"
            
	#check whether the files provide is there.
	if not os.path.exists(whole_gene_list):
	    print "Warning: whole_gene_list:"+whole_gene_list+" is not found in QF_single_end_process.py, exit."
            log_error.write("whole_gene_list not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
	
        if not os.path.exists(genome_fa):
	    print "Warning: genome_fa:"+genome_fa+" is not found in QF_single_end_process.py, exit."
            log_error.write("genome_fa not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(tophat_genome_ref+".1.bt2"):
	    print "Warning: tophat_genome_ref:"+tophat_genome_ref+" is not found in QF_single_end_process.py, exit."
            log_error.write("tophat_genome_ref not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
	if not os.path.exists(file_prefix):
	    print "Warning: intermediate file_prefix:"+file_prefix+" is not found in QF_single_end_process.py, exit."
            log_error.write("intermediate file_prefix not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(bam_fd):
	    print "Warning: folder of bams:"+bam_fd+" is not found in QF_single_end_process.py, exit."
            log_error.write("folder of bams not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
        if not os.path.exists(outresult_fd):
	    print "Warning: output result folder:"+outresult_fd+" is not found in QF_single_end_process.py, exit."
            log_error.write("output result folder not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
                 

        #if not os.path.exists(out_folder):
	    #mkdir the folder
        #   os.makedirs(out_folder)
                
        if not os.path.exists(QF_path):
	    print "Warning: QueryFuse_path:"+QF_path+" is not found in QF_single_end_process.py, exit."
            log_error.write("QueryFuse_path not found in QF_single_end_process.py, exit.\n"); sys.exit(1)
        
	resume_stat_loc=resume_stat
	
        #hardcoded input file names
        SINGLETON_SORT=file_prefix+"singleton_sorted"
        SINGLE_TO_WHOLE_GENE_BED=bam_fd+"singleton_sorted_whole_gene_anno.bed"
        SINGLE_ON_QUERY_BED=SINGLETON_SORT+"_to_query.bed"
        SINGLE_ON_QUERY_ID=SINGLETON_SORT+"_to_query_ID.txt"
        SINGLE_ON_OTHER_BED=SINGLETON_SORT+"_to_other.bed"
	SINGLE_ON_OTHER_BED_SORT=SINGLE_ON_OTHER_BED+"_sort_pl"
        SINGLE_ON_OTHER_ID=SINGLETON_SORT+"_to_other_ID.txt"
        SINGLETON_SORT_BAM=bam_fd+"singleton_sorted.bam"
        SINGLE_ON_OTHER_ONLY_ID=SINGLETON_SORT+"_to_other_only_ID.txt"
        SINGLE_ON_QUERY_ONLY_ID=SINGLETON_SORT+"_to_query_only_ID.txt"
        SINGLE_ON_BOTH_QUERY_OTHER_ID=SINGLETON_SORT+"_to_both_query_and_other_ID.txt"
        SINGLE_ON_QUERY_ONLY_ID_1=SINGLETON_SORT+"_to_query_only_ID_first_mate.txt"
        SINGLE_ON_QUERY_ONLY_ID_2=SINGLETON_SORT+"_to_query_only_ID_second_mate.txt"
        SINGLE_ON_OTHER_ONLY_ID_1=SINGLETON_SORT+"_to_other_only_ID_first_mate.txt"
        SINGLE_ON_OTHER_ONLY_ID_2=SINGLETON_SORT+"_to_other_only_ID_second_mate.txt"
        SINGLE_ON_BOTH_QUERY_OTHER_ID_1=SINGLETON_SORT+"_to_both_query_and_other_ID_first_mate.txt"
        SINGLE_ON_BOTH_QUERY_OTHER_ID_2=SINGLETON_SORT+"_to_both_query_and_other_ID_second_mate.txt"
	UNMAP_PAIR=file_prefix+"unmapped.bam_unmapped_pair"
	QUERY_FA=outresult_fd+"query_gene.fa"
        QUERY_BED=outresult_fd+"query_gene.bed"
	UNMAP_1_BAM=bam_fd+"unmapped.bam_first_mate.bam"
        UNMAP_2_BAM=bam_fd+"unmapped.bam_second_mate.bam"
        UNMAP_1_FA=bam_fd+"unmapped.bam_first_mate.fa"
        UNMAP_2_FA=bam_fd+"unmapped.bam_second_mate.fa"
        UNMAP_1_ON_QUERY_PSL=file_prefix+"unmapped.bam_first_mate_on_query.psl"
        UNMAP_2_ON_QUERY_PSL=file_prefix+"unmapped.bam_second_mate_on_query.psl"
#        UNMAP_1_ON_QUERY_SPLIT_ID=UNMAP_1_ON_QUERY_PSL+"_split_ID.txt"
        UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_BED=UNMAP_1_ON_QUERY_PSL+"_split_ID_subtract.bed"
        UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_ID=UNMAP_1_ON_QUERY_PSL+"_split_ID_subtract_ID.txt"
        UNMAP_1_ON_QUERY_SPAN_MATE_NOT_QUERY_ID=UNMAP_1_ON_QUERY_PSL+"_span_mate_not_query_ID.txt"
#        UNMAP_2_ON_QUERY_SPLIT_ID=UNMAP_2_ON_QUERY_PSL+"_split_ID.txt"
        UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_ID=UNMAP_2_ON_QUERY_PSL+"_split_ID_subtract_ID.txt"
        UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_BED=UNMAP_2_ON_QUERY_PSL+"_split_ID_subtract.bed"
        UNMAP_2_ON_QUERY_SPAN_MATE_NOT_QUERY_ID=UNMAP_2_ON_QUERY_PSL+"_span_mate_not_query_ID.txt"
        
	repMatch_query=int(round(float(1024*11)/size_query))
	repMatch_other=int(round(float(1024*11)/size_other))
	
        #preprocess unmap to query
        step_name="blat query on unmapped first end, in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="blat query on unmapped second end, in QF_single_end_process.py"
	#to make sure this can match the alignment percentage requirement for short alignment as short as 25bp.
	align_mis=int(read_len-round(read_len)*Align_percent/100)
	spe_Align_percent=int(round(25-align_mis)/25*100)
        blat_1_cmd="blat -stepSize="+str(size_query)+" -repMatch="+str(repMatch_query)+" -minIdentity="+str(spe_Align_percent)+" "+QUERY_FA+" "+UNMAP_1_FA+" "+UNMAP_1_ON_QUERY_PSL
        blat_1_cmd_status=QF_all_modules.resume_func(blat_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=blat_1_cmd_status[0]
	log_whole.write("blat -stepSize="+str(size_query)+" -repMatch="+str(repMatch_query)+' is used\n')
	QF_all_modules.optional_step_check(blat_1_cmd_status, log_whole, log_error, "blatting query on unmapped first end, in QF_single_end_process.py")         
	

        step_name="blat query on unmapped second end, in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Generate more unmapped related files of first end, in QF_single_end_process.py"
        blat_2_cmd="blat -stepSize="+str(size_query)+" -repMatch="+str(repMatch_query)+" -minIdentity="+str(spe_Align_percent)+" "+QUERY_FA+" "+UNMAP_2_FA+" "+UNMAP_2_ON_QUERY_PSL
        blat_2_cmd_status=QF_all_modules.resume_func(blat_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=blat_2_cmd_status[0]
	log_whole.write("blat -stepSize="+str(size_query)+" -repMatch="+str(repMatch_query)+' is used\n')
	QF_all_modules.optional_step_check(blat_2_cmd_status, log_whole, log_error, "blatting query on unmapped second end, in QF_single_end_process.py")         
	

        #Generate more unmapped related files.
        step_name="Generate more unmapped related files of first end, in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Generate more unmapped related files of second end, in QF_single_end_process.py"
        psl_to_unmap_1_cmd=QF_path+"QF_psl_to_unmap_process.sh "+UNMAP_1_ON_QUERY_PSL+" "+str(read_len)+" "+QF_path
        psl_to_unmap_1_cmd_status=QF_all_modules.resume_func(psl_to_unmap_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=psl_to_unmap_1_cmd_status[0]
	QF_all_modules.optional_step_check(psl_to_unmap_1_cmd_status, log_whole, log_error, "Generating more unmapped related files of first end, in QF_single_end_process.py")         
	
        
        step_name="Generate more unmapped related files of second end, in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Get singletons aligned to query in QF_single_end_process.py"
        psl_to_unma_2_cmd=QF_path+"QF_psl_to_unmap_process.sh "+UNMAP_2_ON_QUERY_PSL+" "+str(read_len)+" "+QF_path
        psl_to_unma_2_cmd_status=QF_all_modules.resume_func(psl_to_unma_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=psl_to_unma_2_cmd_status[0]
	QF_all_modules.optional_step_check(psl_to_unma_2_cmd_status, log_whole, log_error, "Generating more unmapped related files of second end, in QF_single_end_process.py")         
	

        step_name="Get singletons aligned to query in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Get singletons aligned to others in QF_single_end_process.py"
        #to seperate the SINGLE_TO_WHOLE_GENE_BED into either query or other.
	single_to_query_other_cmd="python "+QF_path+"/grep_from_a_col.py "+SINGLE_TO_WHOLE_GENE_BED+" 17 "+query_ID+" "+SINGLE_ON_QUERY_BED+" "+SINGLE_ON_OTHER_BED
        single_to_query_other_cmd_status=QF_all_modules.resume_func(single_to_query_other_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_to_query_other_cmd_status[0]
	QF_all_modules.key_step_check(single_to_query_other_cmd_status, step_name, log_whole, log_error)
	log_whole.write("finished getting singletons aligned to query"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))+'\n')
        log_whole.write("===================="+'\n')
    
        step_name="Get singletons aligned to others in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Get singletons aligned to both or uniquely to query and others in QF_single_end_process.py"
        single_cut_ID_cmd="cut -f4 "+SINGLE_ON_QUERY_BED+" > "+SINGLE_ON_QUERY_ID+" ; cut -f4 "+SINGLE_ON_OTHER_BED+" > "+SINGLE_ON_OTHER_ID
        single_cut_ID_cmd_status=QF_all_modules.resume_func(single_cut_ID_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_cut_ID_cmd_status[0]
	QF_all_modules.optional_step_check(single_cut_ID_cmd_status, log_whole, log_error, "Getting singletons aligned to others in QF_single_end_process.py")         
	

        step_name="Get singletons aligned to both or uniquely to query and others in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Separate the singleton_only query ID to first and second mate in QF_single_end_process.py"
        single_to_only_other_cmd="python "+QF_path+"/inter_setdiff_geneID_list.py "+SINGLE_ON_OTHER_ID+" "+SINGLE_ON_QUERY_ID+" "+SINGLE_ON_BOTH_QUERY_OTHER_ID+" "+SINGLE_ON_OTHER_ONLY_ID+" "+SINGLE_ON_QUERY_ONLY_ID
        single_to_only_other_cmd_status=QF_all_modules.resume_func(single_to_only_other_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_to_only_other_cmd_status[0]
	QF_all_modules.key_step_check(single_to_only_other_cmd_status, step_name, log_whole, log_error)
	log_whole.write("finished getting singletons aligned to only others"+'\n')
        log_whole.write(time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time())))
        log_whole.write("===================="+'\n')

        #singletons(4)
        step_name="Separate the singleton_only query ID to first and second mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Separate the singleton_only other ID to first and second mate in QF_single_end_process.py"
        sep_only_query_cmd="grep '/1' "+SINGLE_ON_QUERY_ONLY_ID+" | cut -f1 -d '/' > "+SINGLE_ON_QUERY_ONLY_ID_1+" ; grep '/2' "+SINGLE_ON_QUERY_ONLY_ID+" | cut -f1 -d '/' > "+SINGLE_ON_QUERY_ONLY_ID_2
        sep_only_query_cmd_status=QF_all_modules.resume_func(sep_only_query_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=sep_only_query_cmd_status[0]
	QF_all_modules.optional_step_check(sep_only_query_cmd_status, log_whole, log_error, "Separating the singleton_only query ID to first and second mate in QF_single_end_process.py")         
	

        #singletons(5)
        step_name="Separate the singleton_only other ID to first and second mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="Separate the singleton_both ID to first and second mate in QF_single_end_process.py"
        sep_only_other_cmd="grep '/1' "+SINGLE_ON_OTHER_ONLY_ID+" | cut -f1 -d '/' > "+SINGLE_ON_OTHER_ONLY_ID_1+" ; grep '/2' "+SINGLE_ON_OTHER_ONLY_ID+" | cut -f1 -d '/' > "+SINGLE_ON_OTHER_ONLY_ID_2
        sep_only_other_cmd_status=QF_all_modules.resume_func(sep_only_other_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=sep_only_other_cmd_status[0]
	QF_all_modules.optional_step_check(sep_only_other_cmd_status, log_whole, log_error, "Separating the singleton_only other ID to first and second mate in QF_single_end_process.py")         
	
	
        #singletons(6)
        step_name="Separate the singleton_both ID to first and second mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="sort in perl for singleton and unmap_psl files in QF_single_end_process.py"
        sep_both_cmd="grep '/1' "+SINGLE_ON_BOTH_QUERY_OTHER_ID+" | cut -f1 -d '/' > "+SINGLE_ON_BOTH_QUERY_OTHER_ID_1+" ; grep '/2' "+SINGLE_ON_BOTH_QUERY_OTHER_ID+" | cut -f1 -d '/' > "+SINGLE_ON_BOTH_QUERY_OTHER_ID_2
        sep_both_cmd_status=QF_all_modules.resume_func(sep_both_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=sep_both_cmd_status[0]
	QF_all_modules.optional_step_check(sep_both_cmd_status, log_whole, log_error, "Separating the singleton_both ID to first and second mate in QF_single_end_process.py")         
	
	
        step_name="sort in perl for singleton and unmap_psl files in QF_single_end_process.py"
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_split_mate_to_query on first mate in QF_single_end_process.py"
        SINGLE_ON_OTHER_BED_SORT=SINGLE_ON_OTHER_BED+"_sort_pl"
        perl_sort_cmd=QF_path+"/QF_single_end_perl_sort.sh "+UNMAP_1_ON_QUERY_PSL+" "+UNMAP_2_ON_QUERY_PSL+" "+SINGLE_ON_QUERY_BED+" "+SINGLE_ON_OTHER_BED+" "+QF_path
        perl_sort_cmd_status=QF_all_modules.resume_func(perl_sort_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=perl_sort_cmd_status[0]
	QF_all_modules.optional_step_check(perl_sort_cmd_status, log_whole, log_error, "sorting in perl for singleton and unmap_psl files in QF_single_end_process.py")         
	
        
        step_name="process single_end_split_mate_to_query on first mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_split_mate_to_query on second mate in QF_single_end_process.py"
        single_end_split_mate_to_query_1_cmd="python "+QF_path+"/QF_single_end_split_mate_query.py -i "+file_prefix+" -E first -q "+UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_ID+" -Q "+SINGLE_ON_QUERY_ONLY_ID_2+" -u "+UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_BED+" -m "+UNMAP_1_FA+" -U "+UNMAP_1_ON_QUERY_PSL+" -l "+str(read_len)+" -s "+SINGLE_ON_QUERY_BED+" -w "+whole_gene_list+" -t "+tophat_genome_ref+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -a "+str(Align_percent)
        single_end_split_mate_to_query_1_cmd_status=QF_all_modules.resume_func(single_end_split_mate_to_query_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_end_split_mate_to_query_1_cmd_status[0]
	QF_all_modules.optional_step_check(single_end_split_mate_to_query_1_cmd_status, log_whole, log_error, "processing single_end_split_mate_to_query on first mate in QF_single_end_process.py")         
	

        step_name="process single_end_split_mate_to_query on second mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_split_mate_to_other on first mate in QF_single_end_process.py"
        single_end_split_mate_to_query_2_cmd="python "+QF_path+"/QF_single_end_split_mate_query.py -i "+file_prefix+" -E second -q "+UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_ID+" -Q "+SINGLE_ON_QUERY_ONLY_ID_1+" -u "+UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_BED+" -m "+UNMAP_2_FA+" -U "+UNMAP_2_ON_QUERY_PSL+" -l "+str(read_len)+" -s "+SINGLE_ON_QUERY_BED+" -w "+whole_gene_list+" -t "+tophat_genome_ref+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -a "+str(Align_percent)
        single_end_split_mate_to_query_2_cmd_status=QF_all_modules.resume_func(single_end_split_mate_to_query_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_end_split_mate_to_query_2_cmd_status[0]
	QF_all_modules.optional_step_check(single_end_split_mate_to_query_2_cmd_status, log_whole, log_error, "processing single_end_split_mate_to_query on second mate in QF_single_end_process.py")         
	
        
        step_name="process single_end_split_mate_to_other on first mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_split_mate_to_other on second mate in QF_single_end_process.py"
        single_end_split_mate_to_other_1_cmd="python "+QF_path+"/QF_single_end_split_mate_other.py -i "+file_prefix+" -E first -q "+UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_ID+" -Q "+SINGLE_ON_OTHER_ONLY_ID_2+" -u "+UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_BED+" -m "+UNMAP_1_FA+" -U "+UNMAP_1_ON_QUERY_PSL+" -l "+str(read_len)+" -s "+SINGLE_ON_OTHER_BED_SORT+" -w "+whole_gene_list+" -t "+genome_fa+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -a "+str(Align_percent)+" -O "+str(size_other)
        single_end_split_mate_to_other_1_cmd_status=QF_all_modules.resume_func(single_end_split_mate_to_other_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_end_split_mate_to_other_1_cmd_status[0]
	QF_all_modules.optional_step_check(single_end_split_mate_to_other_1_cmd_status, log_whole, log_error, "processing single_end_split_mate_to_other on first mate in QF_single_end_process.py")         
	

        step_name="process single_end_split_mate_to_other on second mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_split_mate_to_both on first mate in QF_single_end_process.py"
        single_end_split_mate_to_other_2_cmd="python "+QF_path+"/QF_single_end_split_mate_other.py -i "+file_prefix+" -E second -q "+UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_ID+" -Q "+SINGLE_ON_OTHER_ONLY_ID_1+" -u "+UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_BED+" -m "+UNMAP_2_FA+" -U "+UNMAP_2_ON_QUERY_PSL+" -l "+str(read_len)+" -s "+SINGLE_ON_OTHER_BED_SORT+" -w "+whole_gene_list+" -t "+genome_fa+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -a "+str(Align_percent)+" -O "+str(size_other)
        single_end_split_mate_to_other_2_cmd_status=QF_all_modules.resume_func(single_end_split_mate_to_other_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_end_split_mate_to_other_2_cmd_status[0]
	QF_all_modules.optional_step_check(single_end_split_mate_to_other_2_cmd_status, log_whole, log_error, "processing single_end_split_mate_to_other on second mate in QF_single_end_process.py")         
	
        
        step_name="process single_end_split_mate_to_both on first mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_split_mate_to_both on second mate in QF_single_end_process.py"
        single_end_split_mate_to_both_1_cmd="python "+QF_path+"/QF_single_end_split_mate_both.py -i "+file_prefix+" -E first -q "+UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_ID+" -Q "+SINGLE_ON_BOTH_QUERY_OTHER_ID_2+" -u "+UNMAP_1_ON_QUERY_SPLIT_ID_SUBTRACT_BED+" -m "+UNMAP_1_FA+" -U "+UNMAP_1_ON_QUERY_PSL+" -l "+str(read_len)+" -s "+SINGLE_ON_QUERY_BED+" -w "+whole_gene_list+" -t "+tophat_genome_ref+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -a "+str(Align_percent)
        single_end_split_mate_to_both_1_cmd_status=QF_all_modules.resume_func(single_end_split_mate_to_both_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_end_split_mate_to_both_1_cmd_status[0]
	QF_all_modules.optional_step_check(single_end_split_mate_to_both_1_cmd_status, log_whole, log_error, "processing single_end_split_mate_to_both on first mate in QF_single_end_process.py")         
	
        step_name="process single_end_split_mate_to_both on second mate in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="process single_end_span_mate_to_other on first mate in QF_single_end_process.py"
        single_end_split_mate_to_both_2_cmd="python "+QF_path+"/QF_single_end_split_mate_both.py -i "+file_prefix+" -E second -q "+UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_ID+" -Q "+SINGLE_ON_BOTH_QUERY_OTHER_ID_1+" -u "+UNMAP_2_ON_QUERY_SPLIT_ID_SUBTRACT_BED+" -m "+UNMAP_2_FA+" -U "+UNMAP_2_ON_QUERY_PSL+" -l "+str(read_len)+" -s "+SINGLE_ON_QUERY_BED+" -w "+whole_gene_list+" -t "+tophat_genome_ref+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path+" -a "+str(Align_percent)
        single_end_split_mate_to_both_2_cmd_status=QF_all_modules.resume_func(single_end_split_mate_to_both_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        resume_stat_loc=single_end_split_mate_to_both_2_cmd_status[0]
	QF_all_modules.optional_step_check(single_end_split_mate_to_both_2_cmd_status, log_whole, log_error, "processing single_end_split_mate_to_both on second mate in QF_single_end_process.py")         
	
         
        #step_name="process single_end_span_mate_to_other on first mate in QF_single_end_process.py"
        #log_whole.write("===================="+'\n')
        #log_whole.write(step_name+'\n')
        #next_step_name="process single_end_span_mate_to_other on second mate in QF_single_end_process.py"
        #single_end_span_mate_to_other_1_cmd="python "+QF_path+"/QF_single_end_span_mate_other.py -i "+file_prefix+" -E first -q "+UNMAP_1_ON_QUERY_SPAN_MATE_NOT_QUERY_ID+" -S "+SINGLE_ON_OTHER_ONLY_ID_2+" -s "+SINGLE_ON_OTHER_BED+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path
        #_status=QF_all_modules.resume_func(single_end_span_mate_to_other_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        #print "finished processing single_end_span_mate_to_other on first mate"
        #print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        #log_whole.write("===================="+'\n')
        
        #step_name="process single_end_span_mate_to_other on second mate in QF_single_end_process.py"
        #log_whole.write("===================="+'\n')
        #log_whole.write(step_name+'\n')
        #next_step_name="process single_end_span_mate_to_both on first mate in QF_single_end_process.py"
        #single_end_span_mate_to_other_2_cmd="python "+QF_path+"/QF_single_end_span_mate_other.py -i "+file_prefix+" -E second -q "+UNMAP_2_ON_QUERY_SPAN_MATE_NOT_QUERY_ID+" -S "+SINGLE_ON_OTHER_ONLY_ID_1+" -s "+SINGLE_ON_OTHER_BED+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path
        #_status=QF_all_modules.resume_func(single_end_span_mate_to_other_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        #print "finished processing single_end_span_mate_to_other on second mate"
        #print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        #log_whole.write("===================="+'\n')
        
        #step_name="process single_end_span_mate_to_both on first mate in QF_single_end_process.py"
        #log_whole.write("===================="+'\n')
        #log_whole.write(step_name+'\n')
        #next_step_name="process single_end_span_mate_to_both on second mate in QF_single_end_process.py"
        #single_end_span_mate_to_both_1_cmd="python "+QF_path+"/QF_single_end_span_mate_other_for_both.py -i "+file_prefix+" -E first -q "+UNMAP_1_ON_QUERY_SPAN_MATE_NOT_QUERY_ID+" -S "+SINGLE_ON_BOTH_QUERY_OTHER_ID_2+" -s "+SINGLE_ON_OTHER_BED+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path
        #_status=QF_all_modules.resume_func(single_end_span_mate_to_both_1_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        #print "finished processing single_end_span_mate_to_both on first mate for single_to_both"
        #print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        #log_whole.write("===================="+'\n')
        #
        #step_name="process single_end_span_mate_to_both on second mate in QF_single_end_process.py"
        #log_whole.write("===================="+'\n')
        #log_whole.write(step_name+'\n')
        #next_step_name="process unmapped pairs in QF_single_end_process.py"
        #single_end_span_mate_to_both_2_cmd="python "+QF_path+"/QF_single_end_span_mate_other_for_both.py -i "+file_prefix+" -E second -q "+UNMAP_2_ON_QUERY_SPAN_MATE_NOT_QUERY_ID+" -S "+SINGLE_ON_BOTH_QUERY_OTHER_ID_1+" -s "+SINGLE_ON_OTHER_BED+" -g "+LOG_F+" -r "+str(resume_stat)+" -F "+QF_path
        #_status=QF_all_modules.resume_func(single_end_span_mate_to_both_2_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        #print "finished processing single_end_span_mate_to_both on second mate for single_to_both"
        #print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
        #log_whole.write("===================="+'\n')

	step_name="process unmapped pairs in QF_single_end_process.py"
        log_whole.write("===================="+'\n')
        log_whole.write(step_name+'\n')
        next_step_name="This is the last step in this script"
	#next_step_name="split read summary in QF_single_end_process.py"
        #unmap_pair_cmd=QF_path+"QF_unmapped_pair.sh"+" "+UNMAP_1_FA+" "+UNMAP_2_FA+" "+file_prefix+" "+QUERY_FA+" "+str(read_len)+" "+whole_gene_list+" "+tophat_genome_ref+" "+`SCC`+" "+LOG_ERR+" "+R_script+" "+QF_path+" "+Perl_script+" "+blat_script_other
        #unmap_pair_cmd_status=QF_all_modules.resume_func(unmap_pair_cmd, resume_stat_loc, step_name, next_step_name, LOG_OUT)
        #resume_stat_loc=unmap_pair_cmd_status[0]
	#QF_all_modules.optional_step_check(unmap_pair_cmd_status, log_whole, log_error, "processing unmapped pairs in QF_single_end_process.py")         
	
	

