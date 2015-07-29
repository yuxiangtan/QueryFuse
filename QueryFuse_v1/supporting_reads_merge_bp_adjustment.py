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
This script is a QueryFuse subfunction in summary process after split/span fusion event filtering.
It is for reference template generation for each reliable event and do breakpoint adjustment.

=============================
Usage: python supporting_reads_merge_bp_adjustment.py
-h help

-o result_folder for the support files										*[No default value]

-O result_folder for the graph files										*[No default value]

-k fusion name for this search 											*[No default value]

-d input data file to search from										*[No default value]

-q FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT 									*[No default value]

-Q sub-group reads file												*[No default value]

-F QF_path													*[No default value]

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

##By Yuxiang Tan
##Contact: yuxiang.tan@gmail.com
##Compatible Python Version:2.7 or above
###Code Framework
	
if __name__ == "__main__":
	import random, sys, os, subprocess, getopt
	from operator import itemgetter, attrgetter
	#it is too slow to load this module
	#from scipy.stats import ks_2samp
	
        if len(sys.argv) < 5:
		print __doc__
		sys.exit(3)


	###set default value
	data_seq=None
	sub_read=None
	old_fusion_key=None
	outsupp_fld=None
	outgraph_fld=None
	FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT=None
	MIN_SCORE=11

	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:B:o:O:w:g:t:T:F:f:l:r:R:P:Y:a:b:Q:S:d:k:i:I:q:m:z')
        for opt in optlist:
		if opt[0] == '-h':
			print __doc__; sys.exit(2)
		##elif opt[0] == '-B': bam_fd = opt[1]
		elif opt[0] == '-o': outsupp_fld = opt[1]
		elif opt[0] == '-O': outgraph_fld = opt[1]
		elif opt[0] == '-F': QF_path = opt[1]
		elif opt[0] == '-Q': sub_read =opt[1]
		elif opt[0] == '-d': data_seq =opt[1]
		elif opt[0] == '-k': old_fusion_key = opt[1]
		elif opt[0] == '-q': FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT = opt[1]
		elif opt[0] == '-m': MIN_SCORE =int(opt[1])

        #for parameter input needed.
	if old_fusion_key==None:
            print "fusion key is not provided in supporting_reads_merge_bp_adjustment.py, exit"
            sys.exit(1)
	    
	if data_seq==None:
            print "grouped_seqs.txt is not provided in supporting_reads_merge_bp_adjustment.py, exit"
            sys.exit(1)
        
	if outsupp_fld==None:
            print "fusion_support folder path is not provided in supporting_reads_merge_bp_adjustment.py, exit"
            sys.exit(1)
        
	if outgraph_fld==None:
            print "fusion_graph folder path is not provided in supporting_reads_merge_bp_adjustment.py, exit"
            sys.exit(1)
        
	if sub_read==None:
            print "sub-group reads file is not provided in supporting_reads_merge_bp_adjustment.py, exit"
            sys.exit(1)
	
	#this is a output name
	if FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT==None:
            print "FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT is not provided in supporting_reads_merge_bp_adjustment.py, exit"
            sys.exit(1)
	
	if not os.path.exists(data_seq):
		print "Warning: "+data_seq+" group_seqs.txt is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(sub_read):
		print "Warning: "+sub_read+" group_reads.txt is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
		
	if not os.path.exists(outsupp_fld):
		print "Warning: "+outsupp_fld+" fusion_support folder is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(outgraph_fld):
		os.mkdir(outgraph_fld)

	#use MIN_SCORE*2 as seed_len to search for seed. move base by base to search at the beginning.
	#record all the read where the seed is, this is also how many spaces it should be in front of the read to align to the ref_temp. In this matrix, what reads supporting this template is also recorded.
	#if all the reads have the same seed (by find function), then find the one have the seed at the most left hand side and the one on the most right hand side, merged 3 part together
	#meanwhile mark down the maximum number of reads that have the share and the reads did not have the seed.
	#If the number of supporting is bigger than 60%, move on.
	#if <80% filter this event out and mark template not available!
	
	#because the shifting of reads are recorded here, I use it to generate supporting graph.
	#because grep bed and fastafrom bed is fast, just generate other.fa each time.



h_seq=open(data_seq,"r")
h_reads=open(sub_read,"r")
seq_matrix_c0=[line.strip().split("\t") for line in h_seq]
if len(seq_matrix_c0)==0:
    sys.exit(3)
read_matrix_c12=[line.strip().split("\t") for line in h_reads]
seq_matrix=sorted(seq_matrix_c0,key=itemgetter(0))
read_matrix_c51=sorted(read_matrix_c12,key=itemgetter(51),reverse=True)
read_matrix=sorted(read_matrix_c51,key=itemgetter(12))


ran_time=0
dic_count={}
dic_start={}
dic_whole={}
dic_read={}
dic_len_shift={}
max_runtime=5
seq_dir={}
seq_end={}
seq_ID={}
read_end={}
read_ID={}
if len(seq_matrix)<max_runtime:
	max_runtime=len(seq_matrix)

for k in range(len(seq_matrix)):
	seq_end[k]=seq_matrix[k][0][-1]
	seq_ID[k]=seq_matrix[k][0][:-2]
	read_end[k]=read_matrix[k][51][-1]
	read_ID[k]=read_matrix[k][12]
	if read_ID[k]!=seq_ID[k]:
		print "The format of two input files are not matched"
		sys.exit(3)
	
	if seq_end[k]==1:
		if read_matrix[k][2]=="+":
			seq_dir[k]="F"
		else:
			seq_dir[k]="R"
	else:
		if read_matrix[k][2]=="+":
			seq_dir[k]="R"
		else:
			seq_dir[k]="F"

ran_seq=random.sample(range(len(seq_matrix)),max_runtime)
	
#run a few random read to find the seed	
while ran_time < max_runtime:
	#try to make sure not picking the same read 
	ran_i=ran_seq[ran_time]
	read_len=len(seq_matrix[ran_i][1])
	dic_count[ran_i]={}
	dic_start[ran_i]={}
	dic_whole[ran_i]={}
	dic_read[ran_i]={}
	dic_len_shift[ran_i]={}
	if read_len<(MIN_SCORE*2):
		print "Their is a read too short to be compare, need to check why!"
		sys.exit(1)
	
	if read_ID[ran_time]!=seq_ID[ran_time]:
		print "The format of two input files are not matched"
		sys.exit(3)
	
	if seq_end[ran_time]==read_end[ran_time]:
		print "The end info does not match for this read"
		ran_time+=1
		#just go to the next one, skip this one
		continue
	
	shift_len=0
	#shift the MIN_SCORE*2bp seed range
	while shift_len<(read_len-MIN_SCORE*2):
		count=0
		if seq_dir[ran_i]=="F":
			sub_str=seq_matrix[ran_i][1][shift_len:(shift_len+MIN_SCORE*2)]
		else:
			sub_str=rc(seq_matrix[ran_i][1][shift_len:(shift_len+MIN_SCORE*2)])
		dic_start[ran_i][shift_len]={}
		dic_read[ran_i][shift_len]={}
		dic_len_shift[ran_i][shift_len]={}
		#scan the seed on each read
		for i in range(len(seq_matrix)):
			if seq_dir[i]=="F":
				start = seq_matrix[i][1].find(sub_str) + 1
			else:
				start = rc(seq_matrix[i][1]).find(sub_str) + 1
			
			if start>0:
				count+=1
				#record the id of the read with that start length
				dic_start[ran_i][shift_len][start]=i
				#record all the read where the seed is, this is also how many spaces it should be in front of the read to align to the ref_temp. In this matrix, what reads supporting this template is also recorded.
				dic_read[ran_i][shift_len][i]=start-1
				if start-1 in dic_len_shift[ran_i][shift_len].keys():
					dic_len_shift[ran_i][shift_len][start-1].append(i)
				else:
					dic_len_shift[ran_i][shift_len][start-1]=[i]
		dic_count[ran_i][count]=shift_len
		if max(dic_count[ran_i].keys())==len(seq_matrix):
			dic_whole[ran_i][count]=dic_start[ran_i][shift_len]
			ran_time=max_runtime
			break
		
		dic_whole[ran_i][count]=dic_start[ran_i][shift_len]
		shift_len+=1
	ran_time+=1


#get the seed with max count
max_check={}
for j in dic_count.keys():
	max_check[max(dic_count[j].keys())]=j


count_max=max(max_check.keys())
count_max_ID=max_check[count_max]

#check whether >60% of reads can shared the max sed
if count_max>(0.6*len(seq_matrix)):
	max_left_key=max(dic_whole[count_max_ID][count_max].keys())
	max_right_key=min(dic_whole[count_max_ID][count_max].keys())
	
	max_left_ID=dic_whole[count_max_ID][count_max][max_left_key]
	max_right_ID=dic_whole[count_max_ID][count_max][max_right_key]
	
	if seq_dir[max_left_ID]=="F":
		seed_seq=seq_matrix[max_left_ID][1][(max_left_key-1):(max_left_key+MIN_SCORE*2-1)]
		left_seq=seq_matrix[max_left_ID][1][:(max_left_key-1)]
	else:
		seed_seq=rc(seq_matrix[max_left_ID][1])[(max_left_key-1):(max_left_key+MIN_SCORE*2-1)]
		left_seq=rc(seq_matrix[max_left_ID][1])[:(max_left_key-1)]
	
	if seq_dir[max_right_ID]=="F":
		right_seq=seq_matrix[max_right_ID][1][(max_right_key+MIN_SCORE*2-1):]
	else:
		right_seq=rc(seq_matrix[max_right_ID][1])[(max_right_key+MIN_SCORE*2-1):]
	
	ref_temp=left_seq+seed_seq+right_seq
	#dic_seq_final=dic_read[ran_i][dic_count[count_max_ID][count_max]]
else:
	print "There is no fusion reference template generated. It may a mix of multiple less reliable fusions for the fusion"+old_fusion_key
	sys.exit(3)


#count the split position p-val
len_ref=len(ref_temp)
if len_ref>read_len:
	distri_start=[]
	for z in range(max_right_key-1,max_left_key):
		if z in dic_len_shift[count_max_ID][dic_count[count_max_ID][count_max]].keys():
			distri_start.insert(z,str(len(dic_len_shift[count_max_ID][dic_count[count_max_ID][count_max]][z])))
		else:
			distri_start.insert(z,str(0))
	#convert the matrix in to a string
	distri_start_string="_".join(distri_start)
	split_pval_cmd="Rscript "+QF_path+"/KS_test_for_split_pval.R in_string="+distri_start_string
	proc = subprocess.Popen(split_pval_cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	split_pval=stdout.strip().split(" ")[1]
	
	##old script for doing ks-test in py but it is too slow
	#mean_uni=float(sum(distri_start))/(max_left_key-max_right_key+1)
	#uniform_temp=[mean_uni]*(max_left_key-max_right_key+1)
	#split_pval=ks_2samp(distri_start, uniform_temp)[0]
		
	old_file_key="_".join(old_fusion_key.split("^")[:-1])+"_"+str(count_max)+"_"+str(split_pval)
	out_key="\t".join(old_fusion_key.split("^")[:-1])+"\t"+str(count_max)+"_"+str(split_pval)
else:
	split_pval=1
	old_file_key="_".join(old_fusion_key.split("^")[:-1])+"_1"+"_"+str(split_pval)
	out_key="\t".join(old_fusion_key.split("^")[:-1])+"\t1"+"_"+str(split_pval)


ID_list=old_fusion_key.split("^")
bp_anno=ID_list[4]+"_"+ID_list[1]+"_"+ID_list[1]+"_"+ID_list[6]+"_"+ID_list[6]+"_"+ID_list[9]
len_head=len(bp_anno)
put_len=len("putative_reference_seq")
space_head=""
for x in range(len_head-put_len):
	space_head=space_head+" "

#output main file
h_ref_filtered=open(FUSION_BREAK_POINT_SUM_COUNT_MERGE_REF_FILT,"a")
h_ref_filtered.write(out_key+"\n")
h_ref_filtered.close()


temp_aln_name=outgraph_fld+"/"+old_file_key+"_alignment_graph.txt"
temp_aln_rev_name=outgraph_fld+"/"+old_file_key+"_alignment_reverse_graph.txt"
temp_fa_all=outsupp_fld+"/all_ref_template.fa"
h_temp_aln=open(temp_aln_name,"w")
h_temp_aln_rev=open(temp_aln_rev_name,"w")
h_temp_fa_all=open(temp_fa_all,"a")
h_temp_fa_all.write(">"+old_file_key+"\n"+ref_temp+"\n")
h_temp_fa_all.close()
#now build the supporting graph using the old key
h_temp_aln.write("putative_reference_seq"+space_head+"^"+ref_temp+"\n")
space_ID={}
#generate alignment graph.
dic_len_key_sort=dic_len_shift[count_max_ID][dic_count[count_max_ID][count_max]].keys()
dic_len_key_sort.sort()
for dic_len_key in dic_len_key_sort:
	for key_ID in dic_len_shift[count_max_ID][dic_count[count_max_ID][count_max]][dic_len_key]:
		space_seq=""
		space_ID[key_ID]=""
		#build the output seq
		for y in range(len_head-len(seq_matrix[key_ID][0])):
			space_ID[key_ID]=space_ID[key_ID]+" "
		for t in range(max_left_key-1-dic_len_key):
			space_seq=space_seq+" "
		if seq_dir[key_ID]=="F":
			out_seq=space_seq+seq_matrix[key_ID][1]
		else:
			out_seq=space_seq+rc(seq_matrix[key_ID][1])
		h_temp_aln.write(seq_matrix[key_ID][0]+space_ID[key_ID]+"^"+out_seq+"\n")


h_temp_aln.close()

#how to get reverse complement
rc_ref_temp=rc(ref_temp)
h_temp_aln_rev.write("putative_reference_seq"+space_head+"^"+rc_ref_temp+"\n")

for dic_len_key in dic_len_key_sort:
	for key_ID in dic_len_shift[count_max_ID][dic_count[count_max_ID][count_max]][dic_len_key]:
		rc_space_seq=""
		#build the output seq
		for t in range(dic_len_key):
			rc_space_seq=rc_space_seq+" "
		if seq_dir[key_ID]=="F":
			rc_out_seq=rc_space_seq+rc(seq_matrix[key_ID][1])
		else:
			rc_out_seq=rc_space_seq+seq_matrix[key_ID][1]
		h_temp_aln_rev.write(seq_matrix[key_ID][0]+space_ID[key_ID]+"\t"+rc_out_seq+"\n")

h_temp_aln_rev.close()

