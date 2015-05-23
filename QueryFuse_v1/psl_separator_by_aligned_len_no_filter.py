# -*- coding: cp936 -*-
"""
The main purpose of this script is to separate the psl file into three part basing on the first column (aligned number) and 2 cutoff number.

=============================
Usage: python psl_separator_by_aligned_len_no_filter.py -i -o 
-h help

-i files (with whole path) to processed         *[No default value]

-a name of output a                       	*[No default value]

-b name of output b				*[No default value]

-c name of output c				*[No default value]

-n cutoff number1				*[No default value]

-N cutoff number2				*[No default value]

-d sep_char within among columns                [default value '\t']

-j skip header lines                            [default value 0]
	if skip header lines = 0 indicates that there is no header line
	if skip header lines = n indicates that first n lines are header lines



===================
input description:
input files: psl file with general format like the psl file from blat.


======================
output files: (if num_1 and num_2 is two number gave as cutoff, num_1 must be <= num_2)
file a: psl with aligned_len <=num_1
file b: psl with aligned_len >num_1 <=num_2
file c: psl with aligned_len >num_2

============================

Python & Module requirement:
Versions: 2.4 or above
Module: No additional Python Module is required.

============================
Library file requirement:
Not Standalone version, few library file is required.
============================
command line example: python /usr3/graduate/ytan7/CBMrepository/utilities/tags/QueryFuse_v2.0/psl_separator_by_aligned_len_no_filter.py -i unmapped.bam_first_mate_on_query.psl -a test_low -b test_mid -c test_up -n 30 -N 54 -j 5

"""

##Copyright
##By Yuxiang Tan
##Contact: yuxiang.tan@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework
'''
command line example: python /usr3/graduate/ytan7/CBMrepository/utilities/tags/QueryFuse_v2.0/psl_separator_by_aligned_len_no_filter.py -i unmapped.bam_first_mate_on_query.psl -a test_low -b test_mid -c test_up -n 30 -N 54 -j 5

'''

if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re
	import os
	import math
	from itertools import ifilter
	
	OUTPUT_SEP_CHAR='\t'
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 7:
		print __doc__
		sys.exit(3)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	infile=None
	outfileA=None
	outfileB=None
	outfileC=None
	cutoff1=None
	cutoff2=None
	sep_char='\t'
	infile_skip=0
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	#all the names can be only one letter!
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:a:b:c:d:j:n:N:z')
	for opt in optlist:
		if opt[0] == '-h':
			print __doc__; sys.exit(2)
		elif opt[0] == '-i': infile = opt[1]
		elif opt[0] == '-a': outfileA = opt[1]
		elif opt[0] == '-b': outfileB = opt[1]
		elif opt[0] == '-c': outfileC = opt[1]
		elif opt[0] == '-n': cutoff1 = int(opt[1])
		elif opt[0] == '-N': cutoff2 = int(opt[1])
		elif opt[0] == '-d': sep_char =opt[1]
		elif opt[0] == '-j': infile_skip= int(opt[1])
		


	#check the parameters
	if infile==None:
		print "infile is not provided in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2)       
	
	if outfileA==None:
		print "outfileA is not provided in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2)       
	
	if outfileB==None:
		print "outfileB is not provided in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2)       
	
	if outfileC==None:
		print "outfileC is not provided in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2)       
	
	if cutoff1==None:
		print "cutoff1 is not provided in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2)
		
	if cutoff2==None:
		print "cutoff2 is not provided in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2) 
	
	if cutoff1>cutoff2:
		print "cutoff1 is greater than cutoff2 which indicate error!!! in psl_separator_by_aligned_len_no_filter.py, exit."; sys.exit(2)
	
	#check whether the files provide is there.
	if not os.path.exists(infile):
		print infile+" not found in whole_summary_filter.py, exit."; sys.exit(2) 
	
	#prepare files.
	h_infile=open(infile,"r")
	h_outfileA=open(outfileA,"w")
	h_outfileB=open(outfileB,"w")
	h_outfileC=open(outfileC,"w")
	
	#to skip the top n lines
	for i in range(infile_skip):
		h_infile.readline()

	#read in line by line and separate.
	for line in h_infile:
		items=line.strip().split("\t")
		#in case the first column of this file is not number as psl
		if items[0].isdigit():
			align_len=int(items[0])
		else:
			sys.exit("The column of the inputfile is not a number in psl_separator_by_aligned_len_no_filter.py, exit.")

		if align_len>cutoff2:
			h_outfileC.write(line)
		else:
			if align_len>cutoff1:
				h_outfileB.write(line)
			else:
				h_outfileA.write(line)
	
	h_infile.close()
	h_outfileA.close()
	h_outfileB.close()
	h_outfileC.close()
	
        
        
        
                


                  
                              
                        
               
                               

		

        
	
	
	
