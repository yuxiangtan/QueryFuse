#dictionary & file I/O are in lab08

import sys

#look at the match of other, then look at the start and end of query, if start and end not close to read_len*98/100, then filter.
#if match of other + end-start of query bigger than read len, then filter.

#check num of parameter
if len(sys.argv) < 5:
	print "Warning! Usage: python eliminate_wrong_split_pair_match_for_summary_merged.py dup_removed_bed_file read_len outputfile"
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
	match_other=int(items[17])
	strand_query=items[26]
	start_query=int(items[29])
	end_query=int(items[30])
	#check whether range over size
	if (end_query-start_query+match_other)<=read_len:
		#check start close to 0 or not.
		if start_query<=(read_len*(100-Align_percent)/100):
			outfile.write(line)
		else:
			#check end close to read_len or not
			if end_query>=(read_len*Align_percent/100):
				outfile.write(line)
			
handle.close()
outfile.close()	

