import sys, os
from operator import itemgetter, attrgetter
import string, subprocess
import time
import random
#This script is to check in a specific column, whether a row contains a string, out put and split the result into contains file and non-contains file.
#input need: column id (start from 1,2 ....) , target file (must be tab deliminated at this stage), string of interest
#output: contains_file, not contains file

#check num of parameter
if len(sys.argv) < 6:
	print "Warning: Usage: python grep_from_a_col.py target_file column_id string_of_interest outfile_contains outfile_non_contains"
	sys.exit(1)
    
target_file = sys.argv[1]
column_id = int(sys.argv[2])
str_target= sys.argv[3]
outfile_contains = sys.argv[4]
outfile_non_contains = sys.argv[5]
#L = int(sys.argv[3])

#print no_module_load
#print module_cmd
#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
#out_split=open(outfile_count,"w")
#out_span=open(outfile_read,"w")

if not os.path.exists(target_file):
	print target_file+" is not found in grep_from_a_col.py, exit"
	sys.exit(1)

h_file=open(target_file,"r")

h_out_con=open(outfile_contains,"w")
h_out_non=open(outfile_non_contains,"w")

for line in h_file:
	items=line.strip().split("\t")
	if line.startswith('#')or not line.split():
		continue
	
	#check whether this string is in this column.
	if str_target in items[column_id-1]:
		h_out_con.write(line)
	else:
		h_out_non.write(line)
	

h_file.close()
h_out_con.close()
h_out_non.close()

