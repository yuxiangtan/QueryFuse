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

#use dictionary to record and count
#chr1_bp1_side1_geneID1_ensemID1 as lv1 (side means F or R, F is 5&+ or 3&-, which means breakpoint on the right, opposite to R)
#chr2_bp2_side2_geneID2_ensemID2 as lv2
#check whether ID1 is in fact a query gene, if yes, filtered it
#check whether two side is too closed to each other (hard coded <50000bp) [The average gene density in the human genome is about one per 40-45 kb of DNA. Assuming a mean size of, say, 10-15 kb, human genes should be separated by about 30 kb of nongenic DNA on average ]
#dic_count records count for each dic
#dic_read records all the reads for each dic
#print out dic1 with dicname + count into one file
#print out dic2 with only reads into read file



#check num of parameter
if len(sys.argv) < 5:
	print "Usage: python fusion_report_summary_merged.py bp_anno_bed_file query.bed outputfile_count outputfile_reads"
	sys.exit(1)
    
datafile = sys.argv[1]
query_bed= sys.argv[2]
outfile_count = sys.argv[3]
outfile_read = sys.argv[4]

#L = int(sys.argv[3])

#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
handle=open(datafile,"r")
han_query=open(query_bed,"r")
outcount=open(outfile_count,"w")
outread=open(outfile_read,"w")

query_line=han_query.readline()
query_items=query_line.strip().split("\t")
chr_q=query_items[0]
start_q=query_items[1]
end_q=query_items[2]
ID_q=query_items[3]
Ens_q=query_items[4]

dic_count={}
dic_read={}
for line in handle:
	items=line.strip().split("\t")
	#read in informations
	chr1=items[0]
	bp1=items[1]
	strand1=items[2]
	end1=items[3]
	direct1=items[4]
	bp2=items[5]
	strand2=items[6]
	end2=items[7]
	direct2=items[8]
	ID1=items[24]
	Ens1=items[25]
	side1="NA"
	side2="NA"
	
	if strand1=="+":
		if direct1=="5":
			side1="F"
		if direct1=="3":
			side1="R"
	if strand1=="-":
		if direct1=="5":
			side1="R"
		if direct1=="3":
			side1="F"
	if strand2=="+":
		if direct2=="5":
			side2="F"
			bp2_g=int(bp2)+int(start_q)
			end2_g=int(end2)+int(start_q)
		if direct2=="3":
			side2="R"
			bp2_g=int(bp2)+int(start_q)
			end2_g=int(end2)+int(start_q)
	if strand2=="-":
		if direct2=="5":
			side2="R"
			bp2_g=int(bp2)+int(start_q)
			end2_g=int(end2)+int(start_q)
		if direct2=="3":
			side2="F"
			bp2_g=int(bp2)+int(start_q)
			end2_g=int(end2)+int(start_q)
	lv1=chr1+"\t"+bp1+"\t"+side1+"\t"+ID1+"\t"+Ens1
	lv2=chr_q+"\t"+str(bp2_g)+"\t"+side2+"\t"+ID_q+"\t"+Ens_q
	
	#check whether ID1 is in fact a query gene, if yes, filtered it
	if ID1 == ID_q:
		continue
	
	#check whether two side is too closed to each other (hard coded <50000bp) 	
	if chr1 == chr_q:
		if abs(int(bp1)-bp2_g)<50000:
			continue
		
	#check whether lv1 (breakpoint on other) exist
	if lv1 not in dic_count.keys():
		#each step create a new dictionary, if not, can not fill content in for the sub-levels
		dic_count[lv1]={}
		dic_count[lv1][lv2]=1
		dic_read[lv1]={}
		dic_read[lv1][lv2]=line
	else:
		#check whether lv2 exist in dic of lv1
		if lv2 not in dic_count[lv1].keys():
			dic_count[lv1][lv2]=1
			dic_read[lv1][lv2]=line
		else:
			dic_count[lv1][lv2]+=1
			dic_read[lv1][lv2]=dic_read[lv1][lv2]+line
		

for lv1_name in dic_count:
	for lv2_name in dic_count[lv1_name]:
		print_line=lv1_name+"\t"+lv2_name+"\t"+str(dic_count[lv1_name][lv2_name])+"\n"
		outcount.write(print_line)

for lv1_name in dic_read:
	for lv2_name in dic_read[lv1_name]:
		#further modify can consider put header to each groupline, such as >>
		outread.write(dic_read[lv1_name][lv2_name])			
			
handle.close()
han_query.close()
outcount.close()
outread.close()

