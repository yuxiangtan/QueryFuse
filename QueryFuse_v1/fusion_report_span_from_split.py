import sys, os
from operator import itemgetter, attrgetter
import string
import time
#first, sort the fusion_break_point_summary_count.txt
#second, sort the paired_sorted_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed
#when sorting sort the query loc1 first, then sort the loc1 of other and mainly on chr of other.
#look at the chr of other in split and then scan through the file with same chr in span.
#get the direction from loc and F/R, and compare with span range.
#if fit, put the rowID into exist dic, and add count and write in the split_span file. (here fit means the boundary can across the breakpoint for 5bp at most. Generally shifting is ~3bp, so use 5bp for safety.)
#if not, check the exist dic, if there, check noexist dic, if also there, delete that dic.
#if not, and also not in the exist dic, then check noexist dic, if not there, then put it in.


#check num of parameter
if len(sys.argv) < 8:
	print "Warning! Usage: python fusion_report_span_from_split.py fusion_break_point_summary_count.txt paired_sorted_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed read_length read_std whole_fusion_sum.txt fusion_split_span_support.txt fusion_span_only.bed "
	sys.exit("Warning! Usage: python fusion_report_span_from_split.py fusion_break_point_summary_count.txt paired_sorted_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed read_length read_std whole_fusion_sum.txt fusion_split_span_support.txt fusion_span_only.bed ")
    
data_split_n = sys.argv[1]
data_span_n = sys.argv[2]
read_len = int(sys.argv[3])
read_std = int(sys.argv[4])
outfile_count = sys.argv[5]
outfile_split_span = sys.argv[6]
outfile_span = sys.argv[7]

#minus std is not consider
if read_std < 0 :
	read_std=0

#L = int(sys.argv[3])

#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
#out_split=open(outfile_count,"w")
#out_span=open(outfile_read,"w")
if not os.path.exists(data_split_n):
	print "Warning: "+data_split_n+" is not found, exit."
        sys.exit("Warning: "+data_split_n+" is not found, exit.")
	    

if not os.path.exists(data_span_n):
	print "Warning: "+data_span_n+" is not found"
	#generate the out_split_span file with no spanning read added.
	h_split=open(data_split_n,"r")
	data_split=[line.strip().split("\t") for line in h_split]
	for row_split_i in data_split:
		if len(row_split_i)==1:
			data_split.remove(row_split_i)
	out_count=open(outfile_count,"w")
	for row_split in range(len(data_split)):
		line_split='\t'.join(data_split[row_split])
		split_num=data_split[row_split][10]
		out_count.write(line_split+"\t0\t"+str(split_num)+"\n")
        sys.exit("Warning: "+data_span_n+" is not found")


#open files
h_split=open(data_split_n,"r")
h_span=open(data_span_n,"r")
out_count=open(outfile_count,"w")
out_split_span=open(outfile_split_span,"w")
out_span=open(outfile_span,"w")

#read the files in
data_split=[line.strip().split("\t") for line in h_split]
data_span=[line.strip().split("\t") for line in h_span]

#to get rid to empty row if it happens.
for row_split_i in data_split:
	if len(row_split_i)==1:
		data_split.remove(row_split_i)


for row_span_i in data_span:
	if len(row_span_i)==1:
		data_span.remove(row_span_i)


#out_split_tem=open("temp_split.txt","w")
#out_span_tem=open("temp_span.txt","w")
	
print("sort split file by col2 and the col1")
print(time.strftime('%Y-%m-%d,%H:%M:%S',time.localtime(time.time())))
data_split_c7=sorted(data_split,key=itemgetter(6))
data_split_c2=sorted(data_split_c7,key=itemgetter(1))
data_split_sort=sorted(data_split_c2,key=itemgetter(0))
#for row_split_j in data_split_sort:
#	out_split_tem.write("%s\n" % row_split_j)	
print("sort span file by col2 and the col1")
print(time.strftime('%Y-%m-%d,%H:%M:%S',time.localtime(time.time())))
data_span_c20=sorted(data_span,key=itemgetter(19))
data_span_c2=sorted(data_span_c20,key=itemgetter(1))
data_span_sort=sorted(data_span_c2,key=itemgetter(0))
#for row_span_j in data_span_sort:
#	out_span_tem.write("%s\n" % row_span_j)

#output the sortedlist into file (no need)
#for i in data_split_sort:
#    k='\t'.join(i)
#    out_split.write(k+"\n")

#for i in data_span_sort:
#    k='\t'.join(i)
#    out_span.write(k+"\n")

#look at the chr of other in split and then scan through the file with same chr in span.
#get the direction from loc and F/R, and compare with span range.
#if fit, put the rowID into exist dic, and add count and write in the split_span file.
#if not, check the exist dic, if there, check noexist dic, if also there, delete that dic.
#if not, and also not in the exist dic, then check noexist dic, if not there, then put it in.

print("match the split and span")
print(time.strftime('%Y-%m-%d,%H:%M:%S',time.localtime(time.time())))
dic_exist={}
dic_count={}
dic_noexist={}
row_ini=0
COUNT=0
row_span=0
for row_split in range(len(data_split_sort)):
	line_split='\t'.join(data_split_sort[row_split])
	chr_other=data_split_sort[row_split][0]
	loc_other=data_split_sort[row_split][1]
	dir_other=data_split_sort[row_split][2]
	chr_query=data_split_sort[row_split][5]
	loc_query=data_split_sort[row_split][6]
	dir_query=data_split_sort[row_split][7]
	shift_range=int(data_split_sort[row_split][12])
	dic_count[line_split]=0
	#To make sure the span list will go back to the top if there two continuous same chr.
	row_span=row_ini
	#for row_span in range(len(data_span_sort)):
	while (row_span < len(data_span_sort)):
		COUNT+=1
		line_span='\t'.join(data_span_sort[row_span])
		chr_other_span=data_span_sort[row_span][0]
		loc1_other_span=data_span_sort[row_span][1]
		loc2_other_span=data_span_sort[row_span][2]	
		chr_query_span=data_span_sort[row_span][18]
		loc1_query_span=data_span_sort[row_span][19]
		loc2_query_span=data_span_sort[row_span][20]
		row_span+=1
		if chr_other == chr_other_span:
			if 'row_last' in dir():
				del row_last
			if dir_other == "F":
				if int(loc2_other_span) >= (int(loc_other)-2*(read_len+read_std)-5-shift_range) and int(loc2_other_span) <= (int(loc_other)+5+shift_range) :
					OTHER="Y"
				else:
					OTHER="N"
			if dir_other == "R":
				if int(loc1_other_span) <= (int(loc_other)+2*(read_len+read_std)+5+shift_range) and int(loc1_other_span) >= (int(loc_other)-5-shift_range) :
					OTHER="Y"
				else:
					OTHER="N"
			if chr_query == chr_query_span:
				if dir_query == "F":
					if int(loc2_query_span) >= (int(loc_query)-2*(read_len+read_std)-5-shift_range) and int(loc2_query_span) <= (int(loc_query)+5+shift_range) :
						QUERY="Y"
					else:
						QUERY="N"
				if dir_query == "R":
					if int(loc1_query_span) <= (int(loc_query)+2*(read_len+read_std)+5+shift_range) and int(loc1_query_span) >= (int(loc_query)-5-shift_range) :
						QUERY="Y"
					else:
						QUERY="N"
			else:
				QUERY="N"
			#if both yes, means fit:
			if QUERY == "Y" and OTHER == "Y":
				if (row_span-1) not in dic_exist.keys():
					dic_exist[row_span-1]={}
					if (row_span-1) in dic_noexist.keys():
						del dic_noexist[row_span-1]
				#it is not in the if statement because I want the spanning reads be shared.
				out_split_span.write(line_split+'\t'+line_span+"\n")
				dic_count[line_split]+=1
			else:
				if (row_span-1) in dic_exist.keys():
					if (row_span-1) in dic_noexist.keys():
						del dic_noexist[(row_span-1)]
				else:
					if (row_span-1) not in dic_noexist.keys(): 
						dic_noexist[row_span-1]=line_span+"\n"
		else:
			if chr_other < chr_other_span:
				#if chr_other > chr_pre_span:
				#	dic_count[line_split]=0
				if (row_span-1) in dic_exist.keys():
					if (row_span-1) in dic_noexist.keys():
						del dic_noexist[row_span-1]
				else:
					if (row_span-1) not in dic_noexist.keys():
						dic_noexist[row_span-1]=line_span+"\n"
				row_last=row_span
				row_span=row_ini		
				break
			else:
				row_ini=row_span-1
				if (row_span-1) not in dic_exist.keys():
					if (row_span-1) not in dic_noexist.keys():
						dic_noexist[row_span-1]=line_span+"\n"
				if 'row_last' in dir():
					del row_last
				#continue
		#row_span+=1
		#chr_pre_span=chr_other_span
#	if row_split > 0:
#		break
		

if 'row_last' in dir():
	print("last_row_number:"+str(row_last))
	for last_rows in range(row_last,len(data_span_sort)):
		line_span='\t'.join(data_span_sort[last_rows])
		dic_noexist[last_rows]=line_span+"\n"
else:
	print("No last_row_number. There are no extra row to process")

print (COUNT)
print(time.strftime('%Y-%m-%d,%H:%M:%S',time.localtime(time.time())))

#build the multi_align dictionary to merge multi-align ones back.
h_split_multi=open(data_split_n[:-4]+"_multi_align_keys.txt","r")
multi_key_dic={}
group_dic={}

for line_multi in h_split_multi:
	multi_line=line_multi.strip().split("\t")
	multi_key="\t".join(multi_line[:-1])
	multi_key_dic[multi_key]=multi_line[-1]
	#filter out events with no spanning reads.should not have it before filtering.
	#if dic_count[multi_key]!=0:
	if multi_line[-1] not in group_dic.keys():
		group_dic[multi_line[-1]]={}
		group_dic[multi_line[-1]][multi_key]=dic_count[multi_key]
	else:
		group_dic[multi_line[-1]][multi_key]=dic_count[multi_key]

h_split_multi.close()

#write the split_span count into outfile
for KEYS in dic_count:
	if KEYS not in multi_key_dic.keys():
		#print the key and edd the count at the end
		items=KEYS.strip().split("\t")
		split_num=items[10]
		split_pval=items[11]
		shift_length=items[12]
		di_entropy=items[13]
		multi_align=items[14]
		sum_num=int(split_num)+int(dic_count[KEYS])
		out_count.write("\t".join(items[:-4])+"\t"+str(dic_count[KEYS])+"\t"+str(sum_num)+"\t"+split_pval+"\t"+shift_length+"\t"+di_entropy+"\t"+multi_align+"\n")	
	else:
		if multi_key_dic[KEYS] in group_dic.keys():
			if KEYS in group_dic[multi_key_dic[KEYS]].keys():
				items=KEYS.strip().split("\t")
				split_num=float(items[10])
				split_pval=items[11]
				shift_length=items[12]
				di_entropy=items[13]
				multi_align_ori=int(items[14])
				multi_align=len(group_dic[multi_key_dic[KEYS]].keys())
				split_n_adj=int(round(split_num*multi_align_ori/multi_align))
				sum_num=split_n_adj+int(dic_count[KEYS])
				#header now is: split_num, split-pval, overlap_len/shift_length,di-entropy,
				out_count.write("\t".join(items[:-5])+"\t"+str(split_n_adj)+"\t"+str(dic_count[KEYS])+"\t"+str(sum_num)+"\t"+split_pval+"\t"+shift_length+"\t"+di_entropy+"\t"+str(multi_align)+"\n")
		else:
			items=KEYS.strip().split("\t")
			split_num=int(round(float(items[10])))
			split_pval=items[11]
			shift_length=items[12]
			di_entropy=items[13]
			multi_align=items[14]
			sum_num=split_num+int(dic_count[KEYS])
			out_count.write("\t".join(items[:-5])+"\t"+str(split_num)+"\t"+str(dic_count[KEYS])+"\t"+str(sum_num)+"\t"+split_pval+"\t"+shift_length+"\t"+di_entropy+"\t"+multi_align+"\n")


#write the span only reads into span_only output
#also, group them and give the count. This will be done by another script.
for KEYS in dic_noexist:
	#print the key and edd the count at the end
	out_span.write(dic_noexist[KEYS])	



h_split.close()
h_span.close()
out_count.close()
out_split_span.close()
out_span.close()

#out_split_tem.close()
#out_span_tem.close()