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

import sys, os
from operator import itemgetter, attrgetter
import string, subprocess
import time
import random
#split the accept.bam into 3 subgroups.
#move the unmapped.bam to the bamfolder.
#sort these bam files by ID

#check num of parameter
if len(sys.argv) < 6:
	print "Warning: Usage: python QF_preprocess.py accept.bam unmapped.bam bam_folder QF_path whole_gene_list"
	sys.exit(1)
    
accept_bam = sys.argv[1]
unmapped_bam = sys.argv[2]
outfd = sys.argv[3]
QF_path = sys.argv[4]
whole_gene_list= sys.argv[5]
#L = int(sys.argv[3])

if not os.path.exists(accept_bam):
	print "Warning: "+accept_bam+" is not found in QF_preprocess.py, exist."
	sys.exit(1)

if not os.path.exists(unmapped_bam):
	print "Warning: "+unmapped_bam+" is not found in QF_preprocess.py, exist."
	sys.exit(1)
	
if not os.path.exists(outfd):
	os.mkdir(outfd)

delimiter = '/'
name_list=accept_bam.split("/")

#split aligned bam
split_bam_cmd="python "+QF_path+"/split_bam_no_align.py "+accept_bam+" -o "+outfd
subprocess.check_call(split_bam_cmd,shell=True)

#sorting bam files
SINGLETON_BAM=outfd+"singleton.bam"
SINGLETON_SORT=outfd+"singleton_sorted"
SINGLETON_SORT_BAM=SINGLETON_SORT+".bam"
SINGLE_TO_WHOLE_GENE_BED=SINGLETON_SORT+"_whole_gene_anno.bed"
PAIR_BAM=outfd+"paired.bam"
PAIR_SORT=outfd+"paired_sorted"
UNMAP_1_BAM=outfd+"unmapped.bam_first_mate.bam"
UNMAP_1_FA=outfd+"unmapped.bam_first_mate.fa"
UNMAP_2_BAM=outfd+"unmapped.bam_second_mate.bam"
UNMAP_2_FA=outfd+"unmapped.bam_second_mate.fa"

print "sort pair_aligned_end"
print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
sort_pair_cmd="samtools sort -n "+PAIR_BAM+" "+PAIR_SORT
subprocess.check_call(sort_pair_cmd,shell=True)

print "sort single_aligned_end"
print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
sort_single_cmd="samtools sort -n "+SINGLETON_BAM+" "+SINGLETON_SORT
subprocess.check_call(sort_single_cmd,shell=True)

#print "sort single_aligned_end"
#print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
annotate_reads_cmd="intersectBed -abam "+SINGLETON_SORT_BAM+" -b "+whole_gene_list+" -bed -wo > "+SINGLE_TO_WHOLE_GENE_BED
subprocess.check_call(annotate_reads_cmd,shell=True)

#this seperate can be replaced by py script. 
print "separate unmap bam into two ends"
print time.strftime('%Y-%m-%d %A %X %Z',time.localtime(time.time()))
get_unmap1_cmd = "bamtools filter -in "+unmapped_bam+" -out "+UNMAP_1_BAM+" -isMateMapped true -isFirstMate true"
get_unmap2_cmd = "bamtools filter -in "+unmapped_bam+" -out "+UNMAP_2_BAM+" -isMateMapped true -isSecondMate true"
conver_fa1_cmd = "bamtools convert -format fasta -in "+UNMAP_1_BAM+" -out "+UNMAP_1_FA
conver_fa2_cmd = "bamtools convert -format fasta -in "+UNMAP_2_BAM+" -out "+UNMAP_2_FA
separate_unmap_cmd=get_unmap1_cmd+"; "+get_unmap2_cmd+"; "+conver_fa1_cmd+"; "+conver_fa2_cmd
subprocess.check_call(separate_unmap_cmd,shell=True)
