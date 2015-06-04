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
This script is wrap out bowtie as aligner and generate related bam file. 
In order to be replacable by other aligner, I use this type of bam as the standard output now

This script will finally generate a bam file
Example: python /usr3/graduate/ytan7/CBMrepository/utilities/tags/QueryFuse_v2.3/bowtie_wrapper.py -i unmapped.bam_first_mate_on_query.psl_split_ID_subtract.bed -I unmapped.bam_second_mate_on_query.psl_split_ID_subtract.bed -f ../../bams/unmapped.bam_first_mate.fa -F ../../bams/unmapped.bam_second_mate.fa -g log_error -o both_unmapped_subtract.fa
=============================
Usage: python bowtie_wrapper.py
-h help

-i in fa file			   					                                                 *[No default value]

-r reference with chr as annotation					                                                 *[No default value]

-o out SAM file                                               								 *[No default value]

-S SCC module                                               								 [default value: 0]

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
	import sys, os, csv, getopt, subprocess
	
	###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 6:
		print __doc__
		sys.exit(3)
        
	###set default value
	#can use the keys in Constant_Libary as default.
	fa_file=None
	genome_ref=None
	out_sam=None
	SCC="0"
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:r:S:o:z')
        for opt in optlist:
		if opt[0] == '-h':
		    print __doc__; sys.exit(2)
		elif opt[0] == '-i': fa_file = opt[1]
		elif opt[0] == '-r': genome_ref = opt[1]
		elif opt[0] == '-S': SCC = opt[1]
		elif opt[0] == '-o': out_sam = opt[1]
		
        if fa_file==None:
		print "fa_file is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if genome_ref==None:
		print "genome_ref is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if out_sam==None:
		print "out_sam is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if not os.path.exists(fa_file):
		print "Warning: "+fa_file+" input fa is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(genome_ref+".1.bt2"):
		print "Warning: "+genome_ref+" genome reference is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if SCC=="0":
		bowtie2_cmd="bowtie2 --very-sensitive-local -f -x "+genome_ref+" -U "+fa_file+" > "+out_sam
        else:
		bowtie2_cmd=SCC+"; bowtie2 --very-sensitive-local -f -x "+genome_ref+" -U "+fa_file+" > "+out_sam
        
	subprocess.call(bowtie2_cmd,shell=True)
