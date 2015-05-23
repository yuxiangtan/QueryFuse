# -*- coding: cp936 -*-
"""
The main purpose of this script is to check the similarity of two ends around the breakpoint region.

=============================
Usage: python whole_summary_similarity_check.py -i -o 
-h help

-i files to processed                           *[No default value]

-I input file path				[default value: current folder]

-o name of output files                       	*[No default value]

-s suffix for the files to be processed         [default value "txt"]

-d sep_char within among columns                [default value '\t']

-j skip header lines                            [default value 0]
	if skip header lines = 0 indicates that there is no header line
	if skip header lines = n indicates that first n lines are header lines

-l unique_id length for the infile		[default value 2]

#-t number of split_reads			[default value 1]
#-a number of span_reads			[default value 1]
#-u number of sum of supporting reads		[default value 2]


===================
input description:
input files: whole_summary_all.txt include all the candidate fusions and supporting reads.


======================
output files: whole_summary_filtered.txt only keep the fusions have pass the filter.


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

#To do this check, first, need to extract n bases around the breakpoint (3n+3n).(because after real-case testing, )
#As a result, need to generate bed file to extract from the genemoe.fa
#Second, use bed tools to get this two n*6 long fa.
#use these two strings from fa and try to compare similarity (which means the use??? by shifting left and right from the breakpoint for n bp, for the next 2n bp (same direction on both end) can get the same m nucletide at that position.(m do not need to be consecutive, for example, in the TBL1_PDL2 defuse 97 case, on TBL1 is TTACCAG, so m=2 for CxG)
#m is the reported similarity. Pick the biggest one as the one to use.
#Note: for both direction, only one direction will have this property.
#example from TL03: TBL1_PDL2 the most reliable one (this model generally work, but not for RP11 case, from that, I found I need 2n extended. This example is still n extended.)
#PDL2 end: n=3: TTTGAT|CTGAAA
#TBL1 end: n=3: GCTTAC|CTGGAA
#As a result, now, the similarity m=3 on right side.
#If shift 1bp to left on PDL2, it will be TTTGA|TCTGAA then n=0, same as left 2 and 3 bp or right 1,2,3bp.
#If shift 1bp to left on TBL1, it will be GCTTA|CCTGGA then n=0, same as others.
#More complicated situation is move on both end. In which, both move 1bp to right will have n=2.
#In all, it is n*2 pairs of fragments comparing to n*2 pairs of fragments (each move have left and right side.) In these n*2^2 pairs, find the biggest m!

#two things need to be careful, first is string, second is direction.
#default of n is 3 because generally this is enough and 3 is the codon length. n should always be a multiple of 3.(If I want to check n bp around, I need 2n more bp as extended in the input.)

#because I will do the fragment extract, I will do the binucletide check here too.

def span_split_check(infile, outfile, sep_char,frag_size,read_len,fragsize_std):
	#check the parameters
	if infile==None:
		print "infile is not provided"; sys.exit(0)       
	
	if outfile==None:
		print "outfile is not provided"; sys.exit(0)       
	
	#check whether the files provide is there.
	if not os.path.exists(infile):
		print "infile not found"; sys.exit(0)
	
	h_infile=open(infile,"r")
	h_outfile=open(outfile,"w")
	
	#read in line by line and filter.
	for line in h_infile:
		items=line.strip().split(sep_char)
		span_n_in=int(items[11])
		split_n_in=int(items[10])
		sum_n_in=int(items[12])
		if span_n_in>=span_n and split_n_in>=split_n and sum_n_in>=sum_n:
			h_outfile.write(line)
	
	h_infile.close()
	h_outfile.close()


if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re
	import os
	import math
	from itertools import ifilter
	
	###Liye own common function,class loading
	#from Constant_Library import *
	#from General_Library import *
	#from File_Class import *
	#from Sequencing_Library import *
	#
	OUTPUT_SEP_CHAR='\t'
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 3:
		print __doc__
		sys.exit(0)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	suffix="txt"
	infile=None
	outfile=None
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
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:I:o:s:d:j:l:t:a:u:z')
	for opt in optlist:
		if opt[0] == '-h':
			print __doc__; sys.exit(0)
		elif opt[0] == '-i': infile = opt[1]
		#elif opt[0] == '-I': input_path =opt[1]
		elif opt[0] == '-o': outfile = opt[1]
		#elif opt[0] == '-s': suffix = opt[1]
		#elif opt[0] == '-d': sep_char =opt[1]
		#elif opt[0] == '-j': infile_skip= int(opt[1])
		#elif opt[0] == '-t': split_n = int(opt[1])
		#elif opt[0] == '-a': span_n = int(opt[1])
		#elif opt[0] == '-u': sum_n =int(opt[1])
		#elif opt[0] == '-l': unique_id_length = int(opt[1])
		#elif opt[0] == '-log_p': log_path =opt[1]


	
	#if not os.path.exists(log_path):
	#	print "log_path not found"; sys.exit(0) 
	
	#prepare files.

	#errorfile=log_path+"error.log"
	#runlogfile=log_path+"run.log"
	#h_errorfile=open(errorfile,"a")
	#h_runfile=open(runlogfile,"a")
	
	
	span_split_check(infile, outfile, sep_char)
	
	#copy the things from shell directly to here.
	#use this to call existing functions.
	#sam_extract_cmd="samtools view "+bam_infile + " " + coor_region +"> "+  outfile_sam
	#subprocess.call(sam_extract_cmd,shell=True)
                
	

        
        
        
                


                  
                              
                        
               
                               

		

        
	
	
	
