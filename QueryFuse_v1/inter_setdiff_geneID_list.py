import sys, os
#from sets import Set

#This script is to get intersect, A-only and B-only from listA and B
#input need: only one col of IDs in file A and file B
#output: intersect, A-only, B-only

#check num of parameter
if len(sys.argv) < 6:
	print "Warning! Usage: python inter_setdiff_geneID_list.py listA listB outlist_intersect outlist_A_only outlist_B_only"
	print "Example: python inter_setdiff_geneID_list.py SINGLE_ON_QUERY_ID SINGLE_ON_OTHER_ID SINGLE_ON_BOTH_QUERY_OTHER_ID SINGLE_ON_QUERY_ONLY_ID SINGLE_ON_OTHER_ONLY_ID"
	sys.exit(1)
    
listA = sys.argv[1]
listB = sys.argv[2]
outlist_intersect= sys.argv[3]
outlist_A_only = sys.argv[4]
outlist_B_only = sys.argv[5]
#L = int(sys.argv[3])

#print no_module_load
#print module_cmd
#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
#out_split=open(outfile_count,"w")
#out_span=open(outfile_read,"w")

if not os.path.exists(listA):
	print "file of A list is not found"
	sys.exit(1)
if not os.path.exists(listB):
	print "file of B list is not found"
	sys.exit(1)
	
#read file into a list:
with open(listA) as f_A:
	lines_A = f_A.read().splitlines()

with open(listB) as f_B:
	lines_B = f_B.read().splitlines()
	
#set up sets
A_set=set(lines_A)
B_set=set(lines_B)

A_B_inter=A_set.intersection(B_set)
A_only_set=A_set.difference(B_set)
B_only_set=B_set.difference(A_set)

#change sets back to list
A_B_inter_list= list(A_B_inter)
A_only_list= list(A_only_set)
B_only_list= list(B_only_set)

#open out file to output
h_out_inter=open(outlist_intersect,"w")
h_out_A=open(outlist_A_only,"w")
h_out_B=open(outlist_B_only,"w")

for line in A_B_inter_list:
	h_out_inter.write(line+"\n")
	
for line in A_only_list:
	h_out_A.write(line+"\n")

for line in B_only_list:
	h_out_B.write(line+"\n")

h_out_inter.close()
h_out_A.close()
h_out_B.close()

