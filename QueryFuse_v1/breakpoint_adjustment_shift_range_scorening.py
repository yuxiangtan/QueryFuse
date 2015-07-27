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
This script is a QueryFuse subfunction in summary process after supporting reads merge
It is for breakpoint correction for merged events also with shift range check.
Also including dinucleotide entropy calculation.
=============================
Usage: python breakpoint_adjustment_shift_range_scorening.py
-h help

-o outfile full name												*[No default value]

-O graph folder, a needed input											*[No default value]

-w whole_gene_list.bed				                                                                *[No default value]

-T tophat_genome_reference_fa - the path of the genome fa file (such as hg19.fa)                                *[No default value]

-g LOG_error file used in the summary function									*[No default value]

-F QueryFuse_path                                                                                               *[No default value]

-l Read length                                                                   				[default value is 99]

-a Align_percent: min percentage of alignment                                                                   [default value is 98]

-Q query_gene.bed file												*[No default value]

-q query_gene.fa file												*[No default value]

-i size step para for blat on query										[default value 5]

-I size step para for blat on other										[default value 11]

-m minscore parameter for blat (which will define the min alignment length allowed)				[default value is 11]

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

def rc(dna):
    import string
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count
	
def entropy_cal(ref_temp_seq,bp_left,range_e,len_DI_COMB,DI_COMB):
    ref_len=len(ref_temp_seq)
    #build seq for dinucletide test by checking range_e bp around the breakpoint on each side.
    if bp_left-range_e-1<0:
	test_seq_l=ref_temp_seq[:(bp_left-1)]
    else:
	test_seq_l=ref_temp_seq[(bp_left-range_e-1):(bp_left-1)]
    if bp_left+range_e-1>ref_len:
	test_seq_r=ref_temp_seq[(bp_left-1):]
    else:
	test_seq_r=ref_temp_seq[(bp_left-1):(bp_left+range_e-1)]
    test_seq=test_seq_l+test_seq_r
    len_test_seq=len(test_seq)
    Entropy=0
    for i_DI in range(len_DI_COMB):
	Pi=round(occurrences(test_seq,DI_COMB[i_DI]))/(len_test_seq-1)
	if Pi>0:
	    Entropy-=Pi*math.log(Pi)
    return Entropy

##By Yuxiang Tan
##Contact: yuxiang.tan@gmail.com
##Compatible Python Version:2.7 or above
###Code Framework
	
if __name__ == "__main__":
    ###Python General Module Import	
#	import sys, csv, getopt, re, subprocess, time, string, random
#	import os
#	from itertools import ifilter

    from operator import itemgetter, attrgetter
    import random, sys, os, subprocess, getopt, copy, shutil, math, copy
    
    if len(sys.argv) < 10:
	print __doc__
	sys.exit(3)


    ###set default value
    #data_seq=None
    #old_fusion_key=None
    #outsupp_fld=None
    out_name=None
    outgraph_fld=None
    whole_gene_list=None
    genome_fa=None
    QUERY_FA=None
    ref_template=None
    query_bed=None
    QF_path=None
    file_prefix=None
    LOG_ERR=None
    R_script="Rscript"
    
    size_query=5
    size_other=11
    Align_percent=98
    read_len=99
    MIN_SCORE="11"
    
    ###get arguments(parameters)
    #all the names of parameters must in the optlist.
    optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:B:o:O:w:g:t:T:F:f:l:r:R:P:Y:a:b:Q:S:d:k:i:I:q:m:z')
    for opt in optlist:
	if opt[0] == '-h':
	    print __doc__; sys.exit(2)
	elif opt[0] == '-f': file_prefix = opt[1]
            ##elif opt[0] == '-B': bam_fd = opt[1]
	elif opt[0] == '-o': out_name = opt[1]
	elif opt[0] == '-O': outgraph_fld = opt[1]
	elif opt[0] == '-w': whole_gene_list = opt[1]
	elif opt[0] == '-g': LOG_ERR = opt[1]
	elif opt[0] == '-t': ref_template =opt[1]
	elif opt[0] == '-T': genome_fa = opt[1]
	elif opt[0] == '-F': QF_path = opt[1]
	elif opt[0] == '-l': read_len = int(opt[1])
	#elif opt[0] == '-r': resume_stat = int(opt[1])
	#elif opt[0] == '-R': R_script =opt[1]
	##elif opt[0] == '-P': Perl_path =opt[1]
	#elif opt[0] == '-Y': Python_path =opt[1]
	elif opt[0] == '-a': Align_percent =int(opt[1])
	#elif opt[0] == '-b': blat_path =opt[1]
	elif opt[0] == '-Q': query_bed =opt[1]
	#elif opt[0] == '-d': data_seq =opt[1]
	#elif opt[0] == '-k': old_fusion_key = opt[1]
	elif opt[0] == '-i': size_query = int(opt[1])
	elif opt[0] == '-I': size_other =int(opt[1])
	elif opt[0] == '-q': QUERY_FA = opt[1]
	elif opt[0] == '-m': MIN_SCORE = opt[1]

    #for parameter input needed.
    if whole_gene_list==None:
	print "whole_gene_list is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
    
    if genome_fa==None:
	print "genome_fa path is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
	
    if out_name==None:
	print "out file name is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
	
    if ref_template==None:
	print "all_ref_template.fa is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
    
    if query_bed==None:
	print "query_gene.bed path is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
    
    if QUERY_FA==None:
	print "query_gene.fa path is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
	
    if QF_path==None:
	print "Warning: QueryFuse_path info is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)  
    
    if LOG_ERR==None:
	print "LOG ERR path is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
    
    if outgraph_fld==None:
	print "fusion_graph folder path is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)
	
    if file_prefix==None:
	print "intermedia folder path is not provided in breakpoint_adjustment_shift_range_scorening.py, exit"
	sys.exit(1)    
                
    if not os.path.exists(whole_gene_list):
	print "Warning: "+whole_gene_list+" whole gene list is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)
    
    if not os.path.exists(genome_fa):
	print "Warning: "+genome_fa+" genome.fa is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)
    
    if not os.path.exists(query_bed):
	print "Warning: "+query_bed+" query_gene.bed is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)
  
    if not os.path.exists(QUERY_FA):
	print "Warning: "+QUERY_FA+" query_gene.fa is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)
	    
    if not os.path.exists(QF_path):
	print "Warning: "+QF_path+" QF_path is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)

    if not os.path.exists(outgraph_fld):
	print "Warning: "+outgraph_fld+" out graph folder path is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)
    
    if not os.path.exists(file_prefix):
	print "Warning: "+file_prefix+" intermediate folder path is not found in breakpoint_adjustment_shift_range_scorening.py, exit."
	sys.exit(1)





    #to correct for multialignment, need to use all fa, if blat to the partner only, it is no way to get multialignment on other genes. As a result, I am using the template to do blat to the whole geneome and to the query also.
    
    #this should be in another file and step.
    repMatch_query=int(round(float(1024*11)/size_query))
    repMatch_other=int(round(float(1024*11)/size_other))
    
    UNMAP_1_ON_QUERY_SPLIT_MATE_OTHER_ID=file_prefix+"/unmapped.bam_first_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_other_ID.txt"
    UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID=file_prefix+"/unmapped.bam_second_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_other_ID.txt"
    UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED=file_prefix+"unmapped.bam_first_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_query_ID_subtract_gene_anno.bed"
    UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED=file_prefix+"unmapped.bam_second_mate_on_query.psl_split_ID_subtract_ID.txt_split_mate_to_query_ID_subtract_gene_anno.bed"
    UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID=file_prefix+"/unmapped.bam_all_split_mate_to_other_ID.txt"
    SINGLE_ON_OTHER_BED=file_prefix+"/singleton_sorted_to_other.bed"
    SINGLE_ON_OTHER_BED_ALL_SPLIT=file_prefix+"/singleton_sorted_to_other_all_split.bed"
    SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE=file_prefix+"/singleton_sorted_to_other_all_split_gene_ID.txt"
    SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT=file_prefix+"/singleton_sorted_to_other_all_split_gene_ID.txt_sort"
    whole_gene_list_split_temp=file_prefix+"/singleton_sorted_to_other_all_split_gene.bed_temp"
    whole_gene_list_split=file_prefix+"/singleton_sorted_to_other_all_split_gene.bed"
    whole_gene_list_split_query_filter=file_prefix+"/singleton_sorted_to_other_all_split_gene_filter_query.bed"
    genome_split_fa=file_prefix+"/singleton_sorted_to_other_all_split_gene.fa"
    
    if os.path.exists(UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID):
	os.remove(UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID)

    #the following need to use subprocess to run
    #step 1 generate gene list bed for only the related other genes.
    #in fact, some of the gene are only come from to query_anno.bed, as a result, they are also considered.
    if os.path.exists(UNMAP_1_ON_QUERY_SPLIT_MATE_OTHER_ID) or os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID) or os.path.exists(UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED) or os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED):
	cmd1=""
	if os.path.exists(UNMAP_1_ON_QUERY_SPLIT_MATE_OTHER_ID):
	    cmd1=cmd1+"cat "+UNMAP_1_ON_QUERY_SPLIT_MATE_OTHER_ID+" >> "+UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID+";"
	    if os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID):
	        cmd1=cmd1+"cat "+UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID+" >> "+UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID+";"
	    cmd1=cmd1+QF_path+"/grep_reads_from_Bed_by_ID.sh "+UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID+" "+SINGLE_ON_OTHER_BED+" "+SINGLE_ON_OTHER_BED_ALL_SPLIT+" "+R_script+" "+QF_path+" "+LOG_ERR+";"
	    cmd1=cmd1+"cut -f17 "+SINGLE_ON_OTHER_BED_ALL_SPLIT+" | sort -u > "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
	    if os.path.exists(UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED):
	        cmd1=cmd1+"cut -f17 "+UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u >> "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
	    if os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED):
	        cmd1=cmd1+"cut -f17 "+UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u >> "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
	    cmd1=cmd1+"sort -u "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+">"+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+";"
	    cmd1=cmd1+QF_path+"/grep_reads_from_whole_gene_list_Bed_by_ID.sh "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+" "+whole_gene_list+" "+whole_gene_list_split_temp+" "+R_script+" "+QF_path+" "+LOG_ERR+";"
	    #to advoid duplicate in gene name, which is not possible in ENSG
	    cmd1=cmd1+"cut -f1,2,3,5 "+whole_gene_list_split_temp+" > "+whole_gene_list_split
	else:
	    if os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID):
		cmd1=cmd1+"cat "+UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID+" >> "+UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID+";"
		cmd1=cmd1+QF_path+"/grep_reads_from_Bed_by_ID.sh "+UNMAP_ALL_ON_QUERY_SPLIT_MATE_OTHER_ID+" "+SINGLE_ON_OTHER_BED+" "+SINGLE_ON_OTHER_BED_ALL_SPLIT+" "+R_script+" "+QF_path+" "+LOG_ERR+";"
		cmd1=cmd1+"cut -f17 "+SINGLE_ON_OTHER_BED_ALL_SPLIT+" | sort -u > "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
		if os.path.exists(UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED):
		    cmd1=cmd1+"cut -f17 "+UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u >> "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
		if os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED):
		    cmd1=cmd1+"cut -f17 "+UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u >> "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
		cmd1=cmd1+"sort -u "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+">"+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+";"
		cmd1=cmd1+QF_path+"/grep_reads_from_whole_gene_list_Bed_by_ID.sh "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+" "+whole_gene_list+" "+whole_gene_list_split_temp+" "+R_script+" "+QF_path+" "+LOG_ERR+";"
		#to advoid duplicate in gene name, which is not possible in ENSG
		cmd1=cmd1+"cut -f1,2,3,5 "+whole_gene_list_split_temp+" > "+whole_gene_list_split
	    else:
		if os.path.exists(UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED):
		    cmd1=cmd1+"cut -f17 "+UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u > "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
		    if os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED):
			cmd1=cmd1+"cut -f17 "+UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u >> "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
		    cmd1=cmd1+"sort -u "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+">"+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+";"
		    cmd1=cmd1+QF_path+"/grep_reads_from_whole_gene_list_Bed_by_ID.sh "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+" "+whole_gene_list+" "+whole_gene_list_split_temp+" "+R_script+" "+QF_path+" "+LOG_ERR+";"
		    #to advoid duplicate in gene name, which is not possible in ENSG
		    cmd1=cmd1+"cut -f1,2,3,5 "+whole_gene_list_split_temp+" > "+whole_gene_list_split
		else:
		    if os.path.exists(UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED):
			cmd1=cmd1+"cut -f17 "+UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED+" | sort -u > "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+";"
			cmd1=cmd1+"sort -u "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE+">"+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+";"
			cmd1=cmd1+QF_path+"/grep_reads_from_whole_gene_list_Bed_by_ID.sh "+SINGLE_ON_OTHER_BED_ALL_SPLIT_GENE_SORT+" "+whole_gene_list+" "+whole_gene_list_split_temp+" "+R_script+" "+QF_path+" "+LOG_ERR+";"
			#to advoid duplicate in gene name, which is not possible in ENSG
			cmd1=cmd1+"cut -f1,2,3,5 "+whole_gene_list_split_temp+" > "+whole_gene_list_split
		    else:
			#all empty
			print "Warning: both "+UNMAP_1_ON_QUERY_SPLIT_MATE_OTHER_ID+" and "+UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID+" and "+UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED+" and "+UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED+" are empty in breakpoint_adjustment_shift_range_scorening.py, use the whole_gene_list.bed."
			#to advoid duplicate in gene name, which is not possible in ENSG
			cmd1="cut -f1,2,3,5 "+whole_gene_list+" > "+whole_gene_list_split
    else:
	print "Warning: both "+UNMAP_1_ON_QUERY_SPLIT_MATE_OTHER_ID+" and "+UNMAP_2_ON_QUERY_SPLIT_MATE_OTHER_ID+" and "+UNMAP_1_ON_QUERY_SPLIT_MATE_QUERY_BED+" and "+UNMAP_2_ON_QUERY_SPLIT_MATE_QUERY_BED+" are empty in breakpoint_adjustment_shift_range_scorening.py, use the whole_gene_list.bed."
	#to advoid duplicate in gene name, which is not possible in ENSG
	cmd1="cut -f1,2,3,5 "+whole_gene_list+" > "+whole_gene_list_split
    
    #remove potential query gene in the list.
    #this is for the scenario that each query.bed has only one query
    handle_query=open(query_bed,"r")
    query_ENSID=handle_query.readline().strip().split("\t")[4]
    handle_query.close()
    cmd1=cmd1+";"+"grep -v "+query_ENSID+" "+whole_gene_list_split+" > "+whole_gene_list_split_query_filter
    
    #blat to other and query
    cmd2="fastaFromBed -fi "+genome_fa+" -bed "+whole_gene_list_split_query_filter+" -fo "+genome_split_fa+" -name"
    
    #Adjust the align_percent to match the min match requirement for min length alignment.
    min_len=int(MIN_SCORE)
    mismatch_num=read_len*(100-Align_percent)/100
    new_align_percent=(100*min_len-100*mismatch_num)/min_len
    #cmd3=blat_script+" -noHead -stepSize="+str(size_query)+" -repMatch="+str(repMatch_query)+" -minIdentity="+str(new_align_percent)+" "+QUERY_FA+" "+ref_template+" "+ref_template+"query.psl"
    cmd3="python "+QF_path+"/local_aligner_wrapper_blat.py -i "+ref_template+" -r "+QUERY_FA+" -s "+str(size_query)+" -a "+str(new_align_percent)+" -o "+ref_template+"query.psl"+" -l "+str(read_len)+" -m "+MIN_SCORE
    #cmd4=blat_script+" -noHead -stepSize="+str(size_other)+" -repMatch="+str(repMatch_other)+" -minIdentity="+str(new_align_percent)+" "+genome_split_fa+" "+ref_template+" "+ref_template+"other.psl"
    cmd4="python "+QF_path+"/local_aligner_wrapper_blat.py -i "+ref_template+" -r "+genome_split_fa+" -s "+str(size_query)+" -a "+str(new_align_percent)+" -o "+ref_template+"other.psl"+" -l "+str(read_len)+" -m "+MIN_SCORE
    
    subprocess.call(cmd1, shell = True)
    subprocess.call(cmd2, shell = True)
    subprocess.call(cmd3, shell = True)
    subprocess.call(cmd4, shell = True)
        
    
    #based on the two blat output, need to check the alignment.
    #first, check whether the alignment to any of the side is longer than len-23
    #second, if there are overlap in between, mark it down (it is the shifting region.)
    #third, take the bp in the middle of the shift region.
    
    query_end=ref_template+"query.psl"
    other_end=ref_template+"other.psl"
    h_query_end=open(query_end,"r")
    h_other_end=open(other_end,"r")
    h_out=open(out_name,"w")
    h_out_blat=open(out_name[:-4]+"_blat_detail.bed","w")
    h_out_filter=open(out_name[:-4]+"_filtered_event.bed","w")
    h_out_file_conn=open(out_name[:-4]+"_file_name_conn.txt","w")
    h_out_file_multi=open(out_name[:-4]+"_multi_align_keys.txt","w")
    #the first mutation is 5 penalty, then 3 after that. 3 on insert into query because it is more like mutation. 1 on insert into target and N count.
    #the start point must be within the alignment percentage range of 0 or 99.
    #the following is a good lesson to learn, if the insert is but enough, it means a split already.
    #99      0       0       0       0       0       1       160170  -       chr13_98654220_R_IPO5_ENSG00000065150_chr3_176853296_F_TBL1XR1_ENSG00000177565_1       99       0       99      chr13   115169878       98494015        98654284       235,64,  0,35,   98494015,98654220,
    #if on minus strand, the breakpoint is at the first loc. if on plus strand, it means forward, bp at second loc.
    #read the whole gene list into dictionary that I can easily access.
    h_whole_gene_split=open(whole_gene_list_split_temp,"r")
    GENE_ANNO={}
    for line_gene in h_whole_gene_split:
	G_line=line_gene.strip().split("\t")
	GENE_ANNO[G_line[4]]=G_line
    
    h_whole_gene_split.close()
    
    h_query_anno=open(query_bed,"r")
    QUERY_ANNO={}
    for line_query in h_query_anno:
	Q_line=line_query.strip().split("\t")
	QUERY_ANNO[Q_line[3]]=Q_line
    
    h_query_anno.close()
    
    #construct the dinucleotide combination.
    DI_COMB=["AA","AC","AT","AG","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
    len_DI_COMB=len(DI_COMB)
    
    query_matrix_c0=[line.strip().split("\t") for line in h_query_end]
    other_matrix_c0=[line.strip().split("\t") for line in h_other_end]
    query_matrix_c9=sorted(query_matrix_c0,key=itemgetter(0),reverse=True)
    other_matrix_c9=sorted(other_matrix_c0,key=itemgetter(0),reverse=True)
    query_matrix=sorted(query_matrix_c9,key=itemgetter(9))
    other_matrix=sorted(other_matrix_c9,key=itemgetter(9))
    
    #just use dic for query_best and other_best with fusion_ID as key. In order to allow multi-alignment, use count as the second key. Use for loop to make combination when
    query_best={}
    pre_query_ID=""
    pre_query_score=""
    count_q=0
    for iq in range(len(query_matrix)):
	#len_q=int(query_matrix[iq][10])
	#start_q=int(query_matrix[iq][11])
	#end_q=int(query_matrix[iq][12])
	#if start_q<=round(len_q)*(1-round(Align_percent)/100) or end_q>=round(len_q)*(round(Align_percent)/100):
	matches=int(query_matrix[iq][0])
	mismatches=int(query_matrix[iq][1])
	mismatchescore=0
	if mismatches>1:
	    mismatchescore=2+mismatches*3
	else:
	    mismatchescore=mismatches*5
	
	nCount=int(query_matrix[iq][3])
	query_insert=int(query_matrix[iq][5])
	n_target_insert=int(query_matrix[iq][6])
	#the first mutation is 5 penalty, then 3 after that. 3 on insert into query because it is more like mutation. 1 on insert into target and N count.
	if n_target_insert>0:
	    targe_insert_len=int(query_matrix[iq][7])
	    if targe_insert_len>=50000:
		n_target_insert=0
	SCORE=matches-mismatchescore-nCount*3-n_target_insert
	if query_matrix[iq][9]!=pre_query_ID:
	    count_q=0
	    query_best[query_matrix[iq][9]]={}
	    query_best[query_matrix[iq][9]][count_q]=query_matrix[iq]
	    pre_query_ID=query_matrix[iq][9]
	    pre_query_score=SCORE
	else:
	    if SCORE>pre_query_score:
		pre_query_score=SCORE
		count_q=0
		query_best[query_matrix[iq][9]]={}
		query_best[query_matrix[iq][9]][count_q]=query_matrix[iq]
	    else:
		if SCORE==pre_query_score:
		    count_q+=1
		    query_best[query_matrix[iq][9]][count_q]=query_matrix[iq]
    
       
    other_best={}
    other_score_dic={}
    pre_other_ID=""
    pre_other_score=""
    count_o=0
    for iq in range(len(other_matrix)):
	#len_o=int(other_matrix[iq][10])
	#start_o=int(other_matrix[iq][11])
	#end_o=int(other_matrix[iq][12])
	#if start_o<=round(len_o)*(1-round(Align_percent)/100) or end_o>=round(len_o)*(round(Align_percent)/100):
	matches=int(other_matrix[iq][0])
	mismatches=int(other_matrix[iq][1])
	mismatchescore=0
	if mismatches>1:
	    mismatchescore=2+mismatches*3
	else:
	    mismatchescore=mismatches*5
	
	nCount=int(other_matrix[iq][3])
	other_insert=int(other_matrix[iq][5])
	n_target_insert=int(other_matrix[iq][6])
	#the first mutation is 5 penalty, then 3 after that. 3 on insert into query because it is more like mutation. 1 on insert into target and N count.
	if n_target_insert>0:
	    targe_insert_len=int(other_matrix[iq][7])
	    if targe_insert_len>=50000:
		n_target_insert=0
	SCORE=matches-mismatchescore-nCount*3-n_target_insert
	if other_matrix[iq][9]!=pre_other_ID:
	    count_o=0
	    other_best[other_matrix[iq][9]]={}
	    other_best[other_matrix[iq][9]][count_o]=other_matrix[iq]
	    other_score_dic[other_matrix[iq][9]]=SCORE
	    pre_other_ID=other_matrix[iq][9]
	    pre_other_score=SCORE
	else:
	    if SCORE>pre_other_score:
		pre_other_score=SCORE
		count_o=0
		other_best[other_matrix[iq][9]]={}
		other_best[other_matrix[iq][9]][count_o]=other_matrix[iq]
		other_score_dic[other_matrix[iq][9]]=SCORE
	    else:
		if SCORE==pre_other_score:
		    count_o+=1
		    other_best[other_matrix[iq][9]][count_o]=other_matrix[iq]
    
    #in order to make sure start from the event with most splitting reads subported.
    other_split_best={}
    for fusion_key in other_best.keys():
	split_num=int(fusion_key.split("_")[-2])
	if not split_num in other_split_best.keys():
	    other_split_best[split_num]={}
	    other_split_best[split_num][fusion_key]=""
	else:
	    other_split_best[split_num][fusion_key]=""
    
    other_split_best_key_sort=other_split_best.keys()
    other_split_best_key_sort.sort(reverse=True)
    
    #use other dic as search
    #first check score by
    #check start(11) and end(12) in other_best and query_best, within alignment range.
    #compare alignment length in other_best with query_best, check shifting
    #assign breakpoint.
    #copy graph using the new fusion breakpoint
    #build file to record the convert between old fusion name and the new one.
    de=-1
    while outgraph_fld[de]=="/":
	de-=1
    
    if (de+1)==0:
	outgraph_fld_final=outgraph_fld+"_final/"
    else:
	outgraph_fld_final=outgraph_fld[:de]+"_final/"
    
    if not os.path.exists(outgraph_fld_final):
	os.mkdir(outgraph_fld_final)
    
    out_dup={}
    group_num=0
    
    #in order to make sure start from the event with most splitting reads subported.
    for split_num in other_split_best_key_sort:
	for fusion_key in other_split_best[split_num].keys():
	    if fusion_key in query_best.keys():
		read_len_o=int(other_best[fusion_key][0][10])
		read_len_q=int(query_best[fusion_key][0][10])
		if read_len_o!=read_len_q:
		    print("Warning: the annotation of two files are not matched, exit")
		    sys.exit(3)
		fusion_score=other_score_dic[fusion_key]
		if fusion_score>(read_len_o-int(MIN_SCORE)):
		    h_out_filter.write(fusion_key+" can align well on the fusion partner gene. It is more likely to be a false positive, skip.\n")
		    continue
		
		#build up dictionary for output usage
		h_out_dic={}
		h_out_count=0
		for sub_key_o in other_best[fusion_key].keys():
		    fusion_o_start=int(other_best[fusion_key][sub_key_o][11])
		    fusion_o_end=int(other_best[fusion_key][sub_key_o][12])
		    fusion_o_len=fusion_o_end-fusion_o_start
		    if other_best[fusion_key][sub_key_o][8]=="-":
			fusion_o_bp=int(other_best[fusion_key][sub_key_o][15])
		    else:
			fusion_o_bp=int(other_best[fusion_key][sub_key_o][16])
		    fusion_o_chr_g=GENE_ANNO[other_best[fusion_key][sub_key_o][13]][0]
		    fusion_o_bp_g=int(fusion_o_bp)+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
		    if fusion_o_len>(read_len_o-int(MIN_SCORE)):
			h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" can align long enough on the fusion partner gene. It is more likely to be a false positive, skip.\n")
			continue
		    if fusion_o_len-fusion_score>15:
			h_out_filter.write(fusion_key+" has too low quality to be reliable, skip.\n")
			continue
			
		    for sub_key_q in query_best[fusion_key].keys():
			fusion_q_start=int(query_best[fusion_key][sub_key_q][11])
			fusion_q_end=int(query_best[fusion_key][sub_key_q][12])
			fusion_q_len=fusion_q_end-fusion_q_start
			#check the start and end first
			if fusion_o_start<=fusion_q_start:
			    if fusion_o_start>round(read_len_o)*(1-round(Align_percent)/100):
				h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" can not be fully aligned. It is more likely to be a false positive, skip.\n")
				continue
			    else:
				if fusion_q_end<round(read_len_o)*(round(Align_percent)/100):
				    h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" can not be fully aligned. It is more likely to be a false positive, skip.\n")
				    continue
			    #the start loction of this blat version should add 1 up. for BLAT v. 34x13
			    overlap_len=fusion_o_end-fusion_q_start
			else:
			    if fusion_q_start>round(read_len_o)*(1-round(Align_percent)/100):
				h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" can not be fully aligned. It is more likely to be a false positive, skip.\n")
				continue
			    else:
				if fusion_o_end<round(read_len_o)*(round(Align_percent)/100):
				    h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" can not be fully aligned. It is more likely to be a false positive, skip.\n")
				    continue
			    overlap_len=fusion_q_end-fusion_o_start
			#if the overlap is too big
			if overlap_len>23:
			    #this criteria is too lose, just restrict it
			    #if (fusion_q_len-overlap_len)<23 or (fusion_o_len-overlap_len)<23:
			    h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" has overlap with longer than 23bp and the remain portion is too short. It is more likely to be a false positive because of highly repetitive region, skip.\n")
			    continue
			
			if query_best[fusion_key][sub_key_q][8]=="-":
			    fusion_q_bp=int(query_best[fusion_key][sub_key_q][15])
			else:
			    fusion_q_bp=int(query_best[fusion_key][sub_key_q][16])
			fusion_q_chr_g=QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][0]
			fusion_q_bp_g=int(fusion_q_bp)+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
			
			if overlap_len>0:
			    bp_relocate=int(overlap_len/2)
			    if fusion_o_start<=fusion_q_start:
				bp_left=fusion_q_start+1+overlap_len-bp_relocate
				if query_best[fusion_key][sub_key_q][8]=="-":
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][16])-overlap_len+bp_relocate+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="F"
				    #the start loction of this blat version should add 1 up. for BLAT v. 34x13
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][15])+1+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				else:
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][15])+1+overlap_len-bp_relocate+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="R"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][16])+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				if other_best[fusion_key][sub_key_o][8]=="-":
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][15])+1+bp_relocate+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="R"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][16])+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				else:
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][16])-bp_relocate+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="F"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][15])+1+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				bp_anno=GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]+"_"+str(bp_o_assign_end)+"_"+str(bp_o_assign)+"_"+str(bp_q_assign)+"_"+str(bp_q_assign_end)+"_"+QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]
				bp_anno_rev=QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]+"_"+str(bp_q_assign_end)+"_"+str(bp_q_assign)+"_"+str(bp_o_assign)+"_"+str(bp_o_assign_end)+"_"+GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]
				
			    else:
				bp_left=fusion_o_start+1+overlap_len-bp_relocate
				if query_best[fusion_key][sub_key_q][8]=="-":
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][15])+1+bp_relocate+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="R"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][16])+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				else:
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][16])-bp_relocate+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="F"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][15])+1+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				if other_best[fusion_key][sub_key_o][8]=="-":
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][16])-overlap_len+bp_relocate+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="F"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][15])+1+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				else:
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][15])+1+overlap_len-bp_relocate+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="R"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][16])+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				bp_anno=QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]+"_"+str(bp_q_assign_end)+"_"+str(bp_q_assign)+"_"+str(bp_o_assign)+"_"+str(bp_o_assign_end)+"_"+GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]
				bp_anno_rev=GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]+"_"+str(bp_o_assign_end)+"_"+str(bp_o_assign)+"_"+str(bp_q_assign)+"_"+str(bp_q_assign_end)+"_"+QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]
			    #write the result out
			    fusion_key_list=fusion_key.split("_")
			    out_key=[fusion_o_chr_g,str(bp_o_assign),dir_o_assign,GENE_ANNO[other_best[fusion_key][sub_key_o][13]][3],GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4],fusion_q_chr_g,str(bp_q_assign),dir_q_assign,QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][3],QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4],fusion_key_list[-2],fusion_key_list[-1],str(overlap_len)] 
			    out_dup_key="_".join([fusion_o_chr_g,str(bp_o_assign),dir_o_assign,fusion_q_chr_g,str(bp_q_assign),dir_q_assign])
			    #print out_dup_key
			    if out_dup_key not in out_dup.keys():
				
				h_graph=open(outgraph_fld+"/"+fusion_key+"_alignment_graph.txt","r") 
				graph_line=h_graph.readline().strip().split("^")
				h_graph.close()
				ref_temp_seq=graph_line[1]
				out_seq=ref_temp_seq[:(bp_left-1)]+"|"+ref_temp_seq[(bp_left-1):]
				#calculate dinucleotide entropy. Basically calculate how many times for each 2bp combination is shown in the template sequence. around the breakpoint!!!!
				Entropy=entropy_cal(ref_temp_seq,bp_left,30,len_DI_COMB,DI_COMB)
				h_out_dic[h_out_count]=copy.deepcopy(out_key)
				h_out_dic[h_out_count].append(str(Entropy))
				
				out_dup[out_dup_key]=out_key
				
				#print("_".join(h_out_dic[h_out_count]))
				#print("\t".join(other_best[fusion_key][sub_key_o]))
				#print("\t".join(query_best[fusion_key][sub_key_q]))
				h_out_blat.write("_".join(h_out_dic[h_out_count])+"\t".join(other_best[fusion_key][sub_key_o])+"\t".join(query_best[fusion_key][sub_key_q])+"\n")
				#build the connect file
				h_out_file_conn.write(outgraph_fld+"/"+fusion_key+" >> "+outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"\n")
				#cp the old graph to new name
				shutil.copyfile(outgraph_fld+"/"+fusion_key+"_alignment_graph.txt", outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_graph.txt")
				shutil.copyfile(outgraph_fld+"/"+fusion_key+"_alignment_reverse_graph.txt", outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_reverse_graph.txt")
				h_graph_out=open(outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_graph.txt","a")
				h_graph_out.write(bp_anno+"^"+out_seq+"\n")
				h_graph_out.close()
				h_graph_rev=open(outgraph_fld+"/"+fusion_key+"_alignment_reverse_graph.txt","r") 
				graph_line_rev=h_graph_rev.readline().strip().split("^")
				ref_temp_seq_rev=graph_line_rev[1]
				out_seq_rev=ref_temp_seq_rev[:(read_len_o-bp_left+1)]+"|"+ref_temp_seq_rev[(read_len_o-bp_left+1):]
				h_graph_out_rev=open(outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_reverse_graph.txt","a")
				h_graph_out_rev.write(bp_anno_rev+"^"+out_seq_rev+"\n")
				h_graph_out_rev.close()
				h_out_count+=1
			else:
			    if overlap_len<round(read_len_o)*(round(Align_percent)/100-1):
				h_out_filter.write(fusion_key+" with bp end in "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_o_chr_g+" "+str(fusion_o_bp_g)+" can not connect with "+other_best[fusion_key][sub_key_o][13]+" at "+fusion_q_chr_g+" "+str(fusion_q_bp_g)+", skip.\n")
				continue
			    if fusion_o_start<=fusion_q_start:
				bp_left=fusion_o_end-overlap_len
				if query_best[fusion_key][sub_key_q][8]=="-":
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][16])+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="F"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][15])+1+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				else:
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][15])+1+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="R"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][16])+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				if other_best[fusion_key][sub_key_o][8]=="-":
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][15])+1+overlap_len+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="R"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][16])+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				else:
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][16])-overlap_len+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="F"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][15])+1+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				bp_anno=GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]+"_"+str(bp_o_assign_end)+"_"+str(bp_o_assign)+"_"+str(bp_q_assign)+"_"+str(bp_q_assign_end)+"_"+QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]
				bp_anno_rev=QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]+"_"+str(bp_q_assign_end)+"_"+str(bp_q_assign)+"_"+str(bp_o_assign)+"_"+str(bp_o_assign_end)+"_"+GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]
				
			    else:
				bp_left=fusion_q_end
				if query_best[fusion_key][sub_key_q][8]=="-":
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][15])+1+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="R"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][16])+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				else:
				    bp_q_assign=int(query_best[fusion_key][sub_key_q][16])+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				    dir_q_assign="F"
				    bp_q_assign_end=int(query_best[fusion_key][sub_key_q][15])+1+int(QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][1])
				if other_best[fusion_key][sub_key_o][8]=="-":
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][16])-overlap_len+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="F"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][15])+1+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				else:
				    bp_o_assign=int(other_best[fusion_key][sub_key_o][15])+1+overlap_len+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				    dir_o_assign="R"
				    bp_o_assign_end=int(other_best[fusion_key][sub_key_o][16])+int(GENE_ANNO[other_best[fusion_key][sub_key_o][13]][1])
				bp_anno=QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]+"_"+str(bp_q_assign_end)+"_"+str(bp_q_assign)+"_"+str(bp_o_assign)+"_"+str(bp_o_assign_end)+"_"+GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]
				bp_anno_rev=GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4]+"_"+str(bp_o_assign_end)+"_"+str(bp_o_assign)+"_"+str(bp_q_assign)+"_"+str(bp_q_assign_end)+"_"+QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4]
				
			    #write the result out
			    fusion_key_list=fusion_key.split("_")
			    out_key=[fusion_o_chr_g,str(bp_o_assign),dir_o_assign,GENE_ANNO[other_best[fusion_key][sub_key_o][13]][3],GENE_ANNO[other_best[fusion_key][sub_key_o][13]][4],fusion_q_chr_g,str(bp_q_assign),dir_q_assign,QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][3],QUERY_ANNO[query_best[fusion_key][sub_key_q][13]][4],fusion_key_list[-2],fusion_key_list[-1],str(overlap_len)] 
			    out_dup_key="_".join([fusion_o_chr_g,str(bp_o_assign),dir_o_assign,fusion_q_chr_g,str(bp_q_assign),dir_q_assign])
			    if out_dup_key not in out_dup.keys():
				
				h_graph=open(outgraph_fld+"/"+fusion_key+"_alignment_graph.txt","r") 
				graph_line=h_graph.readline().strip().split("^")
				h_graph.close()
				#print (outgraph_fld+"/"+fusion_key+"_alignment_graph.txt")
				#print graph_line
				ref_temp_seq=graph_line[1]
				out_seq=ref_temp_seq[:(bp_left-1)]+"|"+ref_temp_seq[(bp_left-1):]
				#calculate dinucleotide entropy. Basically calculate how many times for each 2bp combination is shown in the template sequence. around the breakpoint!!!!
				Entropy=entropy_cal(ref_temp_seq,bp_left,30,len_DI_COMB,DI_COMB)
				h_out_dic[h_out_count]=copy.deepcopy(out_key)
				h_out_dic[h_out_count].append(str(Entropy))
				
				out_dup[out_dup_key]=out_key
				
				#h_out.write("\t".join(out_key)+"\t"+str(overlap_len)+"\n")
				h_out_blat.write("_".join(h_out_dic[h_out_count])+"\t".join(other_best[fusion_key][sub_key_o])+"\t".join(query_best[fusion_key][sub_key_q])+"\n")
				#build the connect file
				h_out_file_conn.write(outgraph_fld+"/"+fusion_key+" >> "+outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"\n")
				#cp the old graph to new name
				shutil.copyfile(outgraph_fld+"/"+fusion_key+"_alignment_graph.txt", outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_graph.txt")
				shutil.copyfile(outgraph_fld+"/"+fusion_key+"_alignment_reverse_graph.txt", outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_reverse_graph.txt")
				h_graph_out=open(outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_graph.txt","a")
				h_graph_out.write(bp_anno+"\t"+out_seq+"\n")
				h_graph_out.close()
				h_graph_rev=open(outgraph_fld+"/"+fusion_key+"_alignment_reverse_graph.txt","r") 
				graph_line_rev=h_graph_rev.readline().strip().split("^")
				ref_temp_seq_rev=graph_line_rev[1]
				out_seq_rev=ref_temp_seq_rev[:(read_len_o-bp_left)]+"|"+ref_temp_seq_rev[(read_len_o-bp_left):]
				h_graph_out_rev=open(outgraph_fld_final+"/"+"_".join(h_out_dic[h_out_count])+"_alignment_reverse_graph.txt","a")
				h_graph_out_rev.write(bp_anno_rev+"\t"+out_seq_rev+"\n")
				h_graph_out_rev.close()
				h_out_count+=1
		h_out_dic_len=len(h_out_dic.keys())
		if h_out_dic_len==1:
		    h_out.write("\t".join(h_out_dic[0])+"\t"+str(h_out_dic_len)+"\n")
		else:
		    if h_out_dic_len>1:
			for dic_i in h_out_dic.keys():
			    #header now is: split_num, split-pval, overlap_len/shift_range,di-entropy,multi-alignment
			    h_out_dic[dic_i][-4]=str(round(int(h_out_dic[dic_i][-4])/round(h_out_dic_len),4))
			    h_out.write("\t".join(h_out_dic[dic_i])+"\t"+str(h_out_dic_len)+"\n")
			    h_out_file_multi.write("\t".join(h_out_dic[dic_i])+"\t"+str(h_out_dic_len)+"\t"+str(group_num)+"\n")
			group_num+=1
		    
	    else:
		h_out_filter.write(fusion_key+" is not in blat to query result, skip.\n")
	
    h_out.close()
    h_out_blat.close()
    h_out_filter.close()
    h_out_file_conn.close()
    h_out_file_multi.close()