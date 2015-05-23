# -*- coding: cp936 -*-
"""
This script is used to replace fasta_from_bed and also divede the fa file into two groups: related to what I need and the rest.

This script will finally generate subtarcted fa file for target reads and the non touched reads into a generally two lines fa file.
Example: python /usr3/graduate/ytan7/CBMrepository/utilities/tags/QueryFuse_v2.3/fa_split_subtract_both_unmap.py -i unmapped.bam_first_mate_on_query.psl_split_ID_subtract.bed -I unmapped.bam_second_mate_on_query.psl_split_ID_subtract.bed -f ../../bams/unmapped.bam_first_mate.fa -F ../../bams/unmapped.bam_second_mate.fa -g log_error -o both_unmapped_subtract.fa
=============================
Usage: python QF_single_end_process.py
-h help

-i path of ummapped first end split subtract.bed			                                                 *[No default value]

-I path of ummapped second end split subtract.bed			                                                 *[No default value]

-f path of ummapped first end fa					                                                 *[No default value]

-F path of ummapped second end fa					                                                 *[No default value]

-g path of error log file                                                                                                *[No default value]

-o out merged fa                                               								 *[No default value, suggested: both_unmapped_subtract.fa]

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
	import sys, os, csv, getopt
	import QF_all_modules
	
	###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 6:
		print __doc__
		sys.exit(3)
        
	###set default value
	#can use the keys in Constant_Libary as default.
	bed_file_1=None
	bed_file_2=None
	fa1=None
	fa2=None
	log_error=None
	out_fa=None
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:I:f:F:g:o:z')
        for opt in optlist:
		if opt[0] == '-h':
		    print __doc__; sys.exit(2)
		elif opt[0] == '-i': bed_file_1 = opt[1]
		elif opt[0] == '-I': bed_file_2 = opt[1]
		elif opt[0] == '-f': fa1 = opt[1]
		elif opt[0] == '-F': fa2 = opt[1]
		elif opt[0] == '-g': log_error = opt[1]
		elif opt[0] == '-o': out_fa = opt[1]
		
        if bed_file_1==None:
		print "bed_file_1 is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if bed_file_2==None:
		print "bed_file_2 is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if fa1==None:
		print "fa1 is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if fa2==None:
		print "fa2 is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if log_error==None:
		print "log_error is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if out_fa==None:
		print "out_fa is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if not os.path.exists(bed_file_1):
		print "Warning: "+bed_file_1+" bed_file_1 is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(bed_file_2):
		print "Warning: "+bed_file_2+" bed_file_2 is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(fa1):
		print "Warning: "+fa1+" fa1 is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(fa2):
		print "Warning: "+fa2+" fa2 is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(log_error):
		print "Warning: "+log_error+" error_log is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	
	#read ID into dictionary
	bed_1_dic=QF_all_modules.subtract_bed_to_dic(bed_file_1)
	bed_2_dic=QF_all_modules.subtract_bed_to_dic(bed_file_2)
	
	#read fa into dictionary
	row_int_fa1=QF_all_modules.fa_to_row_distance(fa1,log_error)
	row_int_fa2=QF_all_modules.fa_to_row_distance(fa2,log_error)
	
	fa_1_dic=QF_all_modules.fa_to_dic(fa1, row_int_fa1)
	fa_2_dic=QF_all_modules.fa_to_dic(fa2, row_int_fa2)
	
	#get the overlap of two sets
	bed_1_set=set(bed_1_dic.keys())
	bed_2_set=set(bed_2_dic.keys())
	
	bed_1_2_set=bed_1_set & bed_2_set
	
	#output subtract fa and the rest fa for each end.
	#because this is for the both_unmapped, I will add the end annotation on the output fa.
	#also no output of the rest fa needed because this is expected to done in the split_to_query step
	h_out=open(out_fa,"w")
	h_out_ID=open(out_fa[:-3]+"_ID.txt","w")
	#these two are not needed, it is just because not the input bed is not correct.
	h_out_ID1=open(out_fa[:-3]+"_ID1.bed","w")
	h_out_ID2=open(out_fa[:-3]+"_ID2.bed","w")
	
	for interset_ID in bed_1_2_set:
		h_out_ID.write(interset_ID+"\n")
		#for end1
		for sub_key in bed_1_dic[interset_ID].keys():
			row_temp=sub_key.split("^")
			#-1 because python use 0 as the start point
			END1=int(row_temp[1])-1
			#did not -1 because this is how to get to the end.
			END2=int(row_temp[2])
			out_ID=">"+interset_ID+"/1"
			sub_fa=fa_1_dic[interset_ID][END1:END2]
			h_out.write(out_ID+"\n"+sub_fa+"\n")
			h_out_ID1.write(interset_ID+"\t"+str(END1+2)+"\t"+str(END2)+"\n")
		#for end2
		for sub_key in bed_2_dic[interset_ID].keys():
			row_temp=sub_key.split("^")
			#-1 because python use 0 as the start point
			END1=int(row_temp[1])-1
			#did not -1 because this is how to get to the end.
			END2=int(row_temp[2])
			out_ID=">"+interset_ID+"/2"
			sub_fa=fa_2_dic[interset_ID][END1:END2]
			h_out.write(out_ID+"\n"+sub_fa+"\n")
			h_out_ID2.write(interset_ID+"\t"+str(END1+2)+"\t"+str(END2)+"\n")
			
	h_out.close()
	h_out_ID.close()