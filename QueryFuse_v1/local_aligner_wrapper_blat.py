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
This script is wrap out blat as aligner and generate related psl file. 
In order to be replacable by other aligner, I use psl as the standard output.
However, the BLAT v. 34x13 has problem on its output location, which is good for fastaFromBed, however, I think the number is not human readable, as a result, I will change it to a uniform, understandable one.
This script will finally generate a psl file
Example: python /usr3/graduate/ytan7/CBMrepository/utilities/tags/QueryFuse_v2.3/local_aligner_wrapper_blat.py 
=============================
Usage: python local_aligner_wrapper_blat.py
-h help

-i in fa file			   					                                                 *[No default value]

-r target fa file							                                                 *[No default value]

-s stepSize value 							                                                 [default value: 11]

-a alignment percentage							                                                 [default value: 98]
		
-o out blat file                                               								 *[No default value]

-l read length                                             								 *[No default value]

-H With header or not													[default: FALSE/TRUE]

-m minscore parameter for blat (which will define the min alignment length allowed)					[default value: 11]
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
	if len(sys.argv) < 5:
		print __doc__
		sys.exit(3)
        
	###set default value
	#can use the keys in Constant_Libary as default.
	fa_file=None
	target_fa=None
	out_psl=None
	read_len=None
	sizeStep="11"
	repmatch="1024"
	Align_percent="98"
	HEADER="FALSE"
	MIN_SCORE=11
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:r:s:a:c:H:l:o:m:z')
        for opt in optlist:
		if opt[0] == '-h':
		    print __doc__; sys.exit(2)
		elif opt[0] == '-i': fa_file = opt[1]
		elif opt[0] == '-r': target_fa = opt[1]
		elif opt[0] == '-s': sizeStep = opt[1]
		elif opt[0] == '-a': Align_percent = opt[1]
		elif opt[0] == '-H': HEADER = opt[1]
		elif opt[0] == '-l': read_len = int(opt[1])
		elif opt[0] == '-m': MIN_SCORE = int(opt[1])
		elif opt[0] == '-o': out_psl = opt[1]
		
        if fa_file==None:
		print "fa_file is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if target_fa==None:
		print "target_fa is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	if out_psl==None:
		print "out_psl is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
		
	if read_len==None:
		print "read length is not provided in supporting_reads_merge_bp_adjustment.py, exit"
		sys.exit(1)
	
	
	if not os.path.exists(fa_file):
		print "Warning: "+fa_file+" input fa is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	if not os.path.exists(target_fa):
		print "Warning: "+target_fa+" target fa is not found in supporting_reads_merge_bp_adjustment.py, exit."
		sys.exit(1)
	
	repmatch=str(int(1024*11/round(int(sizeStep))))
	
	#to make sure this can match the alignment percentage requirement for short alignment as short as MIN_SCORE bp.
	align_mis=int(read_len-round(read_len)*int(Align_percent)/100)
	spe_Align_percent=str(int(round(MIN_SCORE-align_mis)/MIN_SCORE*100))
        
	
		
	if HEADER=="TRUE":
		blat_1_cmd="blat -minScore="+str(MIN_SCORE)+" -stepSize="+sizeStep+" -repMatch="+repmatch+" -minIdentity="+spe_Align_percent+" "+target_fa+" "+fa_file+" "+out_psl
	else:
		blat_1_cmd="blat -minScore="+str(MIN_SCORE)+" -stepSize="+sizeStep+" -repMatch="+repmatch+" -minIdentity="+spe_Align_percent+" -noHead "+target_fa+" "+fa_file+" "+out_psl
	subprocess.check_call(blat_1_cmd,shell=True)


	#the following section should be used when the output format is not the same as BLAT v. 34x13 psl output (it can be annotation different or even because the output from the aligner is not psl.)
	#using the following way, convert the output back to consistant psl.
		#import uuid
		#unique_filename = str(uuid.uuid4())
		#if HEADER=="TRUE":
		#	blat_1_cmd="blat -minScore=11 -stepSize="+sizeStep+" -repMatch="+repmatch+" -minIdentity="+spe_Align_percent+" "+target_fa+" "+fa_file+" "+out_psl+unique_filename
		#else:
		#	blat_1_cmd="blat -minScore=11 -stepSize="+sizeStep+" -repMatch="+repmatch+" -minIdentity="+spe_Align_percent+" -noHead "+target_fa+" "+fa_file+" "+out_psl+unique_filename
		#
		#subprocess.check_call(blat_1_cmd,shell=True)
		##modify the temp output psl into the one same as BLAT v. 34x13
		##need to add the script when there is a version like this, but there is no such version yet.
		#
		##remove the intermedia temp file.
		#os.remove(out_psl+unique_filename)
	