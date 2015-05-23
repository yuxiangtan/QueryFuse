import sys
import string
#read lines, excludes lines with first col smaller or equal to read_len*98/100

#check num of parameter
if len(sys.argv) < 5:
	print "Warning! Usage: python eliminate_mismatch_read_psl.py psl_file read_len outputfile Align_percent"
	sys.exit(1)
    
datafile = sys.argv[1]
read_len = string.atoi(sys.argv[2])
outfile_name = sys.argv[3]
Align_percent = int(sys.argv[4])

#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
handle=open(datafile,"r")
outfile=open(outfile_name,"w")


for line in handle:
	items=line.strip().split("\t")
	match=string.atoi(items[0])
	#check whether match is bigger than read_len*98/100
	if match > read_len*Align_percent/100:
		outfile.write(line)

handle.close()
outfile.close()	

