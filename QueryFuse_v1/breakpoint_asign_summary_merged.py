#dictionary & file I/O are in lab08

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

import sys

#look at the psl, see what strand, + means forward, if 0 on this end, the breakpoint will be the right loc, if 99 on this end, the bp will be the left loci. If on - strand then opposite.
#then based on the side of psl, you can know in bed is the other side and based on strand, you can know which loc (left or right) is the breakpoint.
#new anno is at the front of each row as chr1, bp1, strand, the other end1, direction1 [5 5pi(0) of read is  on this end,3 3pi (99 or read_len) of read is on this end], bp2, strand, the other end2, direct2
#chr2 is not showed here because it is always the query, also loc2 is relative location on query rather than on the genome.

#check num of parameter
if len(sys.argv) < 5:
	print "Usage: python breakpoint_asign_summary_merged.py clean_bed_file read_len outputfile Align_percent"
	sys.exit(1)

datafile = sys.argv[1]
read_len = int(sys.argv[2])
outfile_name = sys.argv[3]
Align_percent = int(sys.argv[4])

#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
handle=open(datafile,"r")
outfile=open(outfile_name,"w")

for line in handle:
	items=line.strip().split("\t")
	strand_other=items[5]
	chr1=items[0]
	loc1_left=items[1]
	loc1_right=items[2]
	match_other=int(items[17])
	strand_query=items[26]
	start_query=int(items[29])
	end_query=int(items[30])
	loc2_left=items[33]
	loc2_right=items[34]
	#in case some situation did not think about
	direct1="NA"
	direct2="NA"
	bp1="NA"
	other_end1="NA"
	bp2="NA"
	other_end2="NA"
	#check start close to 0 or not.
	if start_query<=(read_len*2/100):
		direct2="5"
		direct1="3"
		#check strand then decide breakpoint
		if strand_query=="+":
			bp2=loc2_right
			other_end2=loc2_left
		if strand_query=="-":
			bp2=loc2_left
			other_end2=loc2_right
		if strand_other=="+":
			bp1=loc1_left
			other_end1=loc1_right
		if strand_other=="-":
			bp1=loc1_right
			other_end1=loc1_left
		line_out=chr1+"\t"+bp1+"\t"+strand_other+"\t"+other_end1+"\t"+direct1+"\t"+bp2+"\t"+strand_query+"\t"+other_end2+"\t"+direct2+"\t"+line
		outfile.write(line_out)
	else:
		if end_query>=(read_len*Align_percent/100):
			direct2="3"
			direct1="5"
			if strand_query=="+":
				bp2=loc2_left
				other_end2=loc2_right
			if strand_query=="-":
				bp2=loc2_right
				other_end2=loc2_left
			if strand_other=="+":
				bp1=loc1_right
				other_end1=loc1_left
			if strand_other=="-":
				bp1=loc1_left
				other_end1=loc1_right
			line_out=chr1+"\t"+bp1+"\t"+strand_other+"\t"+other_end1+"\t"+direct1+"\t"+bp2+"\t"+strand_query+"\t"+other_end2+"\t"+direct2+"\t"+line
			outfile.write(line_out)


handle.close()
outfile.close()

