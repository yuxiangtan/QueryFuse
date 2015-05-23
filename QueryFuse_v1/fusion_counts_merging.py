# -*- coding: cp936 -*-
"""
The main purpose of this script is to merge reported events, which are 5bp around, together into one.


=============================
Usage: python fusion_counts_merging.py -c -r -o -O -q -g 
-h help

-c count files to processed                     [default value: fusion_break_point_summary_count.txt]

-r read files to processed                      [default value: fusion_break_point_summary_reads.txt]

-o name of output files                       	[default value: fusion_break_point_summary_count_merged.txt]

-O output folder for subgroups			[default value: fusion_supporting_sub_group]

-q Query gene bed file				[default value: query_gene.bed]

-g grouping range                               [default value: 5]

===================
input description:
input files:
fusion_break_point_summary_count.txt include all the candidate fusions and supporting reads.
fusion_break_point_summary_reads.txt include all the supporting reads for these fusions.
query_gene.bed

======================
output files: fusion_break_point_summary_count_merged.txt only keep the fusions have pass the filter.


============================

Python & Module requirement:
Versions: 2.4 or above
Module: No additional Python Module is required.

============================
Library file requirement:
Not Standalone version, few library file is required.
============================
command line example:

"""

##Copyright
##By Liye Zhang
##Contact: bioliyezhang@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework
'''
command line example: python whole_summary_similarity_check.py -i -o 

'''

def get_row_info(row_num, data_count_sort):
	chrA=data_count_sort[row_num][0]
	bpA=data_count_sort[row_num][1]
	directA=data_count_sort[row_num][2]
	chrB=data_count_sort[row_num][5]
	bpB=data_count_sort[row_num][6]
	directB=data_count_sort[row_num][7]
	return chrA,bpA,directA,chrB,bpB,directB

def renew_pri(row_num, data_count_sort):
	chrA=data_count_sort[row_num][0]
	bpA=data_count_sort[row_num][1]
	directA=data_count_sort[row_num][2]
	chrB=data_count_sort[row_num][5]
	bpB=data_count_sort[row_num][6]
	directB=data_count_sort[row_num][7]
	chrA_pri=chrA
	bpA_pri=bpA
	directA_pri=directA
	chrB_pri=chrB
	bpB_pri=bpB
	directB_pri=directB
	return chrA_pri,bpA_pri,directA_pri,chrB_pri,bpB_pri,directB_pri

def get_subport_reads(row_num, data_count_sort,data_read_sort,out_list):
        chrA=data_count_sort[row_num][0]
        bpA=data_count_sort[row_num][1]
        #directA=data_count_sort[row_num][2]
	ENSGA=data_count_sort[row_num][4]
        chrB=data_count_sort[row_num][5]
        bpB=data_count_sort[row_num][6]
        #directB=data_count_sort[row_num][7]
        #if chrA in read table< chrA, continue, if chrA in read table > chrA, break.
        for row_read_num in range(len(data_read_sort)):
                chrA_read=data_read_sort[row_read_num][0]
                bpA_read=data_read_sort[row_read_num][1]
                bpB_read=data_read_sort[row_read_num][5]
		ENSGA_read=data_read_sort[row_read_num][25]
                if chrA==chrA_read:
                        #check bp, check list empty or not, attend the list
                        if bpA==bpA_read:
                                if int(bpB)==int(bpB_read)+int(headQ):
					if ENSGA==ENSGA_read:
						if out_list=="":
							out_list="\t".join(data_read_sort[row_read_num])+"\n"
						else:
	                                                out_list=out_list+"\t".join(data_read_sort[row_read_num])+"\n"
                else:
                        if chrA<chrA_read:
                                break
        #print(out_list)
        return out_list

if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re, subprocess, time, string, random
	import os
	import math, copy
	from itertools import ifilter
	from operator import itemgetter, attrgetter
	
	###Liye own common function,class loading
	#from Constant_Library import *
	#from General_Library import *
	#from File_Class import *
	#from Sequencing_Library import *
	#
	OUTPUT_SEP_CHAR='\t'
	
	
        ###Start of the program       
	#exit if not enough arguments
	#if len(sys.argv) < 3:
	#	print __doc__
	#	sys.exit(0)
	

	###set default value
	#can use the keys in Constant_Libary as default.
	suffix="txt"
	infile_count="fusion_break_point_summary_count.txt"
	infile_read="fusion_break_point_summary_reads.txt"
        query_bed="query_gene.bed"
	outfile="fusion_break_point_summary_count_merged.txt"
        outfolder="fusion_supporting_sub_group/"
	group_range=5*2
        #log_path=None
	sep_char='\t'
	#sep_gene=','
	#split_n=1
	#span_n=1
	#sum_n=2
	input_path=os.getcwd()
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	#all the names can be only one letter!
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hc:r:q:o:O:g:u:z')
	for opt in optlist:
		if opt[0] == '-h':
			print __doc__; sys.exit(1)
		elif opt[0] == '-c': infile_count = opt[1]
		elif opt[0] == '-r': infile_read = opt[1]
		elif opt[0] == '-q': query_bed =opt[1]
		elif opt[0] == '-o': outfile = opt[1]
		elif opt[0] == '-O': outfolder = opt[1]
		#elif opt[0] == '-s': suffix = opt[1]
		#elif opt[0] == '-d': sep_char =opt[1]
		#elif opt[0] == '-j': infile_skip= int(opt[1])
		#elif opt[0] == '-t': split_n = int(opt[1])
		elif opt[0] == '-g': group_range = int(opt[1])
		#elif opt[0] == '-u': sum_n =int(opt[1])
		#elif opt[0] == '-l': unique_id_length = int(opt[1])
		#elif opt[0] == '-log_p': log_path =opt[1]

        
        #for parameter input needed.
        h_count=open(infile_count,"r")
        h_read=open(infile_read,"r")
	h_query=open(query_bed,"r")
	h_out=open(outfile,"w")
		
        if not os.path.exists(outfolder):
		os.mkdir(outfolder)
		
	#Tasks:
	#Find the events that are close to each other.
	#    read the file in as a matrix
	#    sort by bp and chr of the other gene (col 2 and then col 1)
	#Group by using dictionary function
	#    Scan the sorted matrix, if the chr match, use the first row bp as the key of dictionary. use the bp on query as the sub key and the whole line as content.
	#    If the bp of other gene is closed (10bp, since can consider it as two side)to the first key, check the subkey (10bp also)
	#    If the content follow the same direction (RR RF FR FF), merge them together.
	#Reporting the result after merging
	#    For each key and subkey pair, get the all the content.
	#    sort the supproting read file by chr (if the chr not match, move on)
	#    Scan though the supporting read file to extract all the supporting reads.
	#    Count it and compare it to the sum of the contents (warning if they are not matched)
	#    Use the average of bp in contents as the name of file to put these supporting reads in.
	#    generate fusion_count_merged.txt for all these new fusions.
		    
	#The way I difine direction is in fusion_report_summary_merged.py and breakpoint_asign_summary_merged.py

	#read in the file
	data_count=[line.strip().split("\t") for line in h_count]
	data_read=[line.strip().split("\t") for line in h_read]
	data_query=[line.strip().split("\t") for line in h_query]
	if len(data_query)!=1:
		print "Warning: the query bed file in fusion_counts_merging.py is not farmated as one row as expected , exit."
		sys.exit(1)
	
	#sort the file as wanted (it is already order by other + query)
	data_count_reordered=data_count
	for row_num in range(len(data_count)):
		data_count_reordered[row_num][1]=int(data_count[row_num][1])
	data_read_reordered=data_read
	for row_num in range(len(data_read)):
		data_read_reordered[row_num][1]=int(data_read[row_num][1])
	
	data_count_c2=sorted(data_count_reordered,key=itemgetter(1))
	#order by ensem_ID
	data_count_c1=sorted(data_count_c2,key=itemgetter(0))
	data_count_sort_all=data_count_c1
	data_count_sort=data_count_c1
	data_read_c2=sorted(data_read_reordered,key=itemgetter(1))
	#order by ensem_ID
	data_read_c1=sorted(data_read_c2,key=itemgetter(0))
	data_read_sort_all=data_read_c1
	data_read_sort=data_read_c1
	
	for row_num in range(len(data_count)):
		data_count_sort[row_num][1]=str(data_count_sort_all[row_num][1])
	
	for row_num in range(len(data_read)):
		data_read_sort[row_num][1]=str(data_read_sort_all[row_num][1])
		
	#read in query
	chrQ=data_query[0][0]
	headQ=data_query[0][1]
	endQ=data_query[0][2]
	
	
	dic_out={}
	dic_out_read={}
	for row_num in range(len(data_count_sort)):
		#read the first row and build the basic dic
		if row_num==0:
			chrA_pri,bpA_pri,directA_pri,chrB_pri,bpB_pri,directB_pri=renew_pri(row_num, data_count_sort)
			dic_out[chrA_pri]={}
			dic_out[chrA_pri][bpA_pri]={}
			dic_out[chrA_pri][bpA_pri][bpB_pri]={}
			dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri]={}
			dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]={}
			dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri][row_num]=data_count_sort[row_num]
			#get the supporting reads
			dic_out_read[chrA_pri]={}
			dic_out_read[chrA_pri][bpA_pri]={}
			dic_out_read[chrA_pri][bpA_pri][bpB_pri]={}
			dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri]={}
			dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]=""
			dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri])
			if chrB_pri!=chrQ:
				print "Warning: the query bed file in fusion_counts_merging.py does not match the query gene in fusion_break_point_summary_count , exit."
				sys.exit(1)
		#for the rest rows, need to compare to existing dictionaries
		else:
			chrA,bpA,directA,chrB,bpB,directB=get_row_info(row_num, data_count_sort)
			if chrB!=chrQ:
				print "Warning: the query bed file in fusion_counts_merging.py does not match the query gene in fusion_break_point_summary_count , exit."
				sys.exit(1)
				
			#Scan the sorted matrix, if the chr match, use the first row bp as the key of dictionary. use the bp on query as the sub key, then directionA, then directB and row_number as the 3rd subkey, the whole line as content.
			if chrA==chrA_pri:
				if (int(bpA)-int(bpA_pri))<11:
					#bpB is not sorted, so need to use abs. also, this can be multi locations, as a result, it is a for loop to look for keys
					bpB_key_exist=0
					for bpB_key in dic_out[chrA_pri][bpA_pri].keys():
						if abs(int(bpB)-int(bpB_key))<11:
							if directA==directA_pri:
								if directB==directB_pri:
									#merge
									dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
									dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB])
									bpB_key_exist=bpB_key_exist+1
								#because it is not match to this directB_pri perfectly,
								#but it can mathc to potentially another directB
								else:
									#check whether it fits another directB key:
									if directB in dic_out_read[chrA_pri][bpA_pri][bpB_key][directA].keys():
										dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
										dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB])
										bpB_key_exist=bpB_key_exist+1
									#old test script
									#if directB in dic_out_read[chrA_pri][bpA_pri][bpB_key][directA].keys():
									#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
									#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB])
									#else:
									#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB]={}
									#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
									#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=""
									#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA_pri][directB_pri])
							else:
								if directA in dic_out_read[chrA_pri][bpA_pri][bpB_key].keys():
									if directB in dic_out_read[chrA_pri][bpA_pri][bpB_key][directA].keys():
										dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
										dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB])
										bpB_key_exist=bpB_key_exist+1
									#else:
									#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB]={}
									#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
									#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=""
									#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA_pri][directB_pri])
								#else:
								#	dic_out[chrA_pri][bpA_pri][bpB_key][directA]={}
								#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB]={}
								#	dic_out[chrA_pri][bpA_pri][bpB_key][directA][directB][row_num]=data_count_sort[row_num]
								#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA]={}
								#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=""
								#	dic_out_read[chrA_pri][bpA_pri][bpB_key][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_key][directA_pri][directB_pri])
					#if no key in bpB_key can fit. then just build a new one
					if bpB_key_exist==0:
						if bpB not in dic_out[chrA_pri][bpA_pri].keys():
							dic_out[chrA_pri][bpA_pri][bpB]={}
						if directA not in dic_out[chrA_pri][bpA_pri][bpB].keys():
							dic_out[chrA_pri][bpA_pri][bpB][directA]={}
						if directB not in dic_out[chrA_pri][bpA_pri][bpB][directA].keys():
							dic_out[chrA_pri][bpA_pri][bpB][directA][directB]={}
						dic_out[chrA_pri][bpA_pri][bpB][directA][directB][row_num]=data_count_sort[row_num]
						if bpB not in dic_out_read[chrA_pri][bpA_pri].keys():
							dic_out_read[chrA_pri][bpA_pri][bpB]={}
						if directA not in dic_out_read[chrA_pri][bpA_pri][bpB].keys():
							dic_out_read[chrA_pri][bpA_pri][bpB][directA]={}
						dic_out_read[chrA_pri][bpA_pri][bpB][directA][directB]=""
						dic_out_read[chrA_pri][bpA_pri][bpB][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB][directA][directB])	
				else:
					chrA_pri,bpA_pri,directA_pri,chrB_pri,bpB_pri,directB_pri=renew_pri(row_num, data_count_sort)
					#if bpA_pri not in dic_out[chrA_pri].keys():
					dic_out[chrA_pri][bpA_pri]={}
					dic_out[chrA_pri][bpA_pri][bpB_pri]={}
					dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri]={}
					dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]={}
					dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri][row_num]=data_count_sort[row_num]
					dic_out_read[chrA_pri][bpA_pri]={}
					dic_out_read[chrA_pri][bpA_pri][bpB_pri]={}
					dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA]={}
					dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA][directB]=""
					dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA][directB]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB][directA][directB])
					if chrB_pri!=chrQ:
						print "Warning: the query bed file in fusion_counts_merging.py does not match the query gene in fusion_break_point_summary_count , exit."
						sys.exit(1)
				
			else:
				chrA_pri,bpA_pri,directA_pri,chrB_pri,bpB_pri,directB_pri=renew_pri(row_num, data_count_sort)
				dic_out[chrA_pri]={}
				dic_out[chrA_pri][bpA_pri]={}
				dic_out[chrA_pri][bpA_pri][bpB_pri]={}
				dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri]={}
				dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]={}
				dic_out[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri][row_num]=data_count_sort[row_num]
				#get the supporting reads
				dic_out_read[chrA_pri]={}
				dic_out_read[chrA_pri][bpA_pri]={}
				dic_out_read[chrA_pri][bpA_pri][bpB_pri]={}
				dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri]={}
				dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]=""
				dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri]=get_subport_reads(row_num, data_count_sort,data_read_sort,dic_out_read[chrA_pri][bpA_pri][bpB_pri][directA_pri][directB_pri])
				if chrB_pri!=chrQ:
					print "Warning: the query bed file in fusion_counts_merging.py does not match the query gene in fusion_break_point_summary_count , exit."
					sys.exit(1)
			
	#after merged in dictionaries, the problem now is how to output it. Especially where is the bp should be assigned.
	#output the count_file:
	for chrA_key in dic_out.keys():
		for bpA_key in dic_out[chrA_key].keys():
			for bpB_key in dic_out[chrA_key][bpA_key].keys():
				for directA_key in dic_out[chrA_key][bpA_key][bpB_key].keys():
					for directB_key in dic_out[chrA_key][bpA_key][bpB_key][directA_key].keys():
						#if only one row, just output it in groups
						if len(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key])==1:
							row_key=dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key].keys()[0]
							h_out.write("\t".join(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key])+"\n")
							out_name="_".join(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key])
							h_out_group_ID=open(outfolder+out_name+"grouped_ID.txt","w")
							h_out_group_ID.write("\t".join(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key])+"\n")
							h_out_group_reads=open(outfolder+out_name+"grouped_reads.txt","w")
							h_out_group_reads.write(dic_out_read[chrA_key][bpA_key][bpB_key][directA_key][directB_key])
							h_out_group_ID.close()
							h_out_group_reads.close()
						#if more the one, need to merge
						else:
							counter=0
							bpA_counter=0
							bpB_counter=0
							split_merge=0
							#I clound not find better way then just take average because it is no way to check whether it is correct breakpoint just by this
							for row_key in dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key].keys():
								bpA=int(copy.deepcopy(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key][1]))
								bpB=int(copy.deepcopy(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key][6]))
								split_num=int(copy.deepcopy(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key][10]))
								bpA_counter=bpA_counter+bpA
								bpB_counter=bpB_counter+bpB
								split_merge=split_merge+split_num
								counter=counter+1
							bpA_merge=str(int(bpA_counter/counter))
							bpB_merge=str(int(bpB_counter/counter))
							temp_name=copy.deepcopy(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key])
							temp_name[1]=bpA_merge
							temp_name[6]=bpB_merge
							temp_name[10]=str(split_merge)
							out_name="_".join(temp_name)
							h_out_group_ID=open(outfolder+out_name+"grouped_ID.txt","w")
							for row_key in dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key].keys():
								h_out_group_ID.write("\t".join(dic_out[chrA_key][bpA_key][bpB_key][directA_key][directB_key][row_key])+"\n")
							h_out_group_reads=open(outfolder+out_name+"grouped_reads.txt","w")
							h_out_group_reads.write(dic_out_read[chrA_key][bpA_key][bpB_key][directA_key][directB_key])
							h_out_group_ID.close()
							h_out_group_reads.close()
							h_out.write("\t".join(temp_name)+"\n")
							
	


	
	
	h_out.close()
	
