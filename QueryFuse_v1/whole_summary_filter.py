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
The main purpose of this script is to filter the input summary file by cutoffs on split, span and total read number.

=============================
Usage: python whole_summary_filter.py -i -o 
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

-t number of split_reads			[default value 1]
-a number of span_reads			[default value 1]
-u number of sum of supporting reads		[default value 2]


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
##By Yuxiang Tan
##Contact: yuxiang.tan@gmail.com
##Compatible Python Version:2.4 or above

###Code Framework
'''
command line example: python whole_summary_filter.py -i -o 

'''

if __name__ == "__main__":
	###Python General Module Import	
	import sys, csv, getopt, re
	import os
	import math
	from itertools import ifilter
	
	##Liye own common function,class loading
	from Constant_Library import *
	from General_Library import *
	from File_Class import *
	from Sequencing_Library import *
	
	OUTPUT_SEP_CHAR='\t'
	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 3:
		print __doc__
		sys.exit(3)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	suffix="txt"
	infile=None
	outfile=None
	#log_path=None
	sep_char='\t'
	#sep_gene=','
	split_n=1
	span_n=1
	sum_n=2
	input_path=os.getcwd()
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	#all the names can be only one letter!
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:I:o:s:d:j:l:t:a:u:z')
	for opt in optlist:
		if opt[0] == '-h':
			print __doc__; sys.exit(2)
		elif opt[0] == '-i': infile = opt[1]
		#elif opt[0] == '-I': input_path =opt[1]
		elif opt[0] == '-o': outfile = opt[1]
		#elif opt[0] == '-s': suffix = opt[1]
		#elif opt[0] == '-d': sep_char =opt[1]
		#elif opt[0] == '-j': infile_skip= int(opt[1])
		elif opt[0] == '-t': split_n = int(opt[1])
		elif opt[0] == '-a': span_n = int(opt[1])
		elif opt[0] == '-u': sum_n =int(opt[1])
		#elif opt[0] == '-l': unique_id_length = int(opt[1])
		#elif opt[0] == '-log_p': log_path =opt[1]

	#check the parameters
	if infile==None:
		print infile+" is not provided in whole_summary_filter.py, exit."; sys.exit(2)       
	
	if outfile==None:
		print outfile+" is not provided in whole_summary_filter.py, exit."; sys.exit(2)       
	
	#if log_path==None:
	#	print "log_path is not provided"; sys.exit(0)       
	
	#check whether the files provide is there.
	if not os.path.exists(infile):
		print infile+" not found in whole_summary_filter.py, exit."; sys.exit(2) 
	
	#if not os.path.exists(log_path):
	#	print "log_path not found"; sys.exit(0) 
	
	#prepare files.
	h_infile=open(infile,"r")
	h_outfile=open(outfile,"w")
	#errorfile=log_path+"error.log"
	#runlogfile=log_path+"run.log"
	#h_errorfile=open(errorfile,"a")
	#h_runfile=open(runlogfile,"a")
	
	
	#read in line by line and filter.
	for line in h_infile:
		items=line.strip().split("\t")
		span_n_in=int(items[11])
		split_n_in=int(items[10])
		sum_n_in=int(items[12])
		if span_n_in>=span_n and split_n_in>=split_n and sum_n_in>=sum_n:
			h_outfile.write(line)
	
	h_infile.close()
	h_outfile.close()
	
	#copy the things from shell directly to here.
	#use this to call existing functions.
	#sam_extract_cmd="samtools view "+bam_infile + " " + coor_region +"> "+  outfile_sam
	#subprocess.call(sam_extract_cmd,shell=True)
                
	

        
        
        
                


                  
                              
                        
               
                               

		

        
	
	
	
