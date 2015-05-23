# -*- coding: cp936 -*-
"""
This script is a subfunction in result summary in QueryFuse.
It use the sorted span file as input, it should have two section, first one is aligned to other genes, second one is aligned to query gene.

The result should be sort by gene name.
In each gene group, order by loc, use read length*3 as searching range. Read end with strand as direction. Now, query end is not sorting by any of these and just group by direction.
=============================
Usage: python QF_span_only_grouping.py
-h help

-i input file of span only reads                                                                                    *[No default value, but generally should be outfd/fusion_span_only.bed]

-o output file index                                                                                                *[No default value, but generally should be outfd/span_grouped_summary without suffix]

-g LOG_folder                                                                                                       *[No default value]

-f filter for span supporting reads number									    [default value 2]

-l Read_length - length of reads		                                                                    [default value 99]

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
	import sys, csv, getopt, re, subprocess, time, string, random
	import os
	from itertools import ifilter
        from operator import itemgetter, attrgetter

	
	
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 3:
		print __doc__
		sys.exit(0)
	
	###set default value
	#can use the keys in Constant_Libary as default.
	file_in=None
        file_out=None
	LOG_F=None
        
        read_len="99"
	filter_num="2"
	
	###get arguments(parameters)
	#all the names of parameters must in the optlist.
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:o:g:l:f:z')
        for opt in optlist:
		if opt[0] == '-h':
		    print __doc__; sys.exit(0)
		elif opt[0] == '-i': file_in = opt[1]
		elif opt[0] == '-o': file_out = opt[1]
		elif opt[0] == '-g': LOG_F = opt[1]
		elif opt[0] == '-l': read_len = opt[1]
		elif opt[0] == '-f': filter_num = opt[1]

   
        if LOG_F==None:
		print "Warning: LOG_F is not provided"; sys.exit(0)
        
        if not os.path.exists(LOG_F):
		os.mkdir(LOG_F)
                
        LOG_ERR=LOG_F+"/error.log"
        LOG_OUT=LOG_F+"/out.log"
        
        log_error=open(LOG_ERR,"a")
        
        #for parameter input needed.
        if file_in==None:
		print "Warning: file_in is not provided in QF_span_only_grouping.py, exit."
		log_error.write("Warning: file_in is not provided\n"); sys.exit(1)
        
        if file_out==None:
		print "Warning: file_out is not provided in QF_span_only_grouping.py, exit."
		log_error.write("Warning: file_out is not provided\n"); sys.exit(1)
                
                
	if not os.path.exists(file_in):
		print "Warning: "+file_in+" file_in not found in QF_span_only_grouping.py, exit."
		log_error.write("Warning: "+file_in+" file_in not found in QF_span_only_grouping.py, exit.\n"); sys.exit(1)

	if os.stat(file_in).st_size==0:
		print "Note: "+file_in+" is empty which means no more spanning reads left in QF_span_only_grouping.py, exit."
		log_error.write("Note: "+file_in+" is empty which means no more spanning reads left in QF_span_only_grouping.py, exit.\n"); sys.exit(0)

	
        #if not os.path.exists(out_folder):
	    #mkdir the folder
        #   os.makedirs(out_folder)
                
#        if not os.path.exists(QF_path):
#		print "QueryFuse_path not found"
#		log_error.write("QueryFuse_path not found\n"); sys.exit(1)

        data_in = file_in
        group_range = int(read_len)*4 
        outfile = file_out
       
#test location:
#data_in="/restricted/projectnb/montilab-p/LinGA_unprotected/ytan/TL_PCNSL/align/TL_03_pairend_tophat_out/gene_specific_fusion_KIAA1432/SUMMARY/fusion_span_only.bed"
#outfile = "/restricted/projectnb/montilab-p/LinGA_unprotected/ytan/TL_PCNSL/align/TL_03_pairend_tophat_out/gene_specific_fusion_KIAA1432/SUMMARY/fusion_span_only_grouped_summary"
	h_span=open(data_in,"r")
	h_out=open(outfile+".txt","w")
	h_out_detail=open(outfile+"_detail.txt","w")
	h_out_filter=open(outfile+"_filtered.txt","w")
	h_out_detail_filter=open(outfile+"_detail_filtered.txt","w")
	
	#read in the file
	data_span=[line.strip().split("\t") for line in h_span]
	#sort the file as wanted (it is already order by other + query)
	data_span_reordered=data_span
	#col need to be redefined.
	#print("sort span file by col2 and 17")
	#print(time.localtime(time.time()))
	#order by loc
	data_span_c2=sorted(data_span_reordered,key=itemgetter(1))
	#order by ensem_ID
	data_span_c17=sorted(data_span_c2,key=itemgetter(16))
	data_span_sort=data_span_c17
	
	
	#read line by line to put into groups (not consider query end loc, assume for each group should have only one fusion break)
	#within the searching range (read_len*2*2 because of two side. Group all the reads within this range into groups of diff direction (F means breakpoint on the right and R means breakpoint on the left)
	dic_out={}
	dic_detail={}
	dic_out_final={}
	dic_detail_final={}
	gene_query_pre=data_span_sort[0][34]
	for row_num in range(len(data_span_sort)):
		line_split='\t'.join(data_span_sort[row_num])+"\n"
		#read the first row
		if row_num==0:
			gene_id_pre=data_span_sort[row_num][16]
			#check gene id whether is esemble ID begins with ENS
			if gene_id_pre[0:3]!="ENS":
				print("The ID is not esembleID")
				break
			loc_pri=data_span_sort[row_num][1]
			#check whether the loc is bigger than 0
			if int(loc_pri)<0:
				print("The coordinate col is not correct")
				break
			end_other=(data_span_sort[row_num][3].split("/"))[1]
			str_other=data_span_sort[row_num][5]
			end_query=(data_span_sort[row_num][21].split("/"))[1]
			str_query=data_span_sort[row_num][23]
			if str_other=="+":
				if end_other=="1":
					dir_other="F"
					if end_query=="2":
						if str_query=="+":
							dir_query="R"
						elif str_query=="-":
							dir_query="F"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
				elif end_other=="2":
					dir_other="R"
					if end_query=="1":
						if str_query=="+":
							dir_query="F"
						elif str_query=="-":
							dir_query="R"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
			elif str_other=="-":
				if end_other=="1":
					dir_other="R"
					if end_query=="2":
						if str_query=="+":
							dir_query="R"
						elif str_query=="-":
							dir_query="F"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
				elif end_other=="2":
					dir_other="F"
					if end_query=="1":
						if str_query=="+":
							dir_query="F"
						elif str_query=="-":
							dir_query="R"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
			else:
				print"This row: "+line_split.strip("\n")+" has strand problem."
			
			if (dir_other=="F" or dir_other=="R") and (dir_query=="F" or dir_query=="R"):
				dic_out[gene_id_pre]={}
				dic_detail[gene_id_pre]={}
				dic_name = loc_pri+"_"+dir_other+"_"+dir_query
				dic_out[gene_id_pre][dic_name]=1
				dic_detail[gene_id_pre][dic_name]=line_split
			else:
				print"Can not extract directional information from this row: "+line_split.strip("\n")
		else:
			loc_other=data_span_sort[row_num][1]
			gene_id=data_span_sort[row_num][16]
			end_other=(data_span_sort[row_num][3].split("/"))[1]
			str_other=data_span_sort[row_num][5]
			end_query=(data_span_sort[row_num][21].split("/"))[1]
			str_query=data_span_sort[row_num][23]
			if str_other=="+":
				if end_other=="1":
					dir_other="F"
					if end_query=="2":
						if str_query=="+":
							dir_query="R"
						elif str_query=="-":
							dir_query="F"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
				elif end_other=="2":
					dir_other="R"
					if end_query=="1":
						if str_query=="+":
							dir_query="F"
						elif str_query=="-":
							dir_query="R"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
			elif str_other=="-":
				if end_other=="1":
					dir_other="R"
					if end_query=="2":
						if str_query=="+":
							dir_query="R"
						elif str_query=="-":
							dir_query="F"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
				elif end_other=="2":
					dir_other="F"
					if end_query=="1":
						if str_query=="+":
							dir_query="F"
						elif str_query=="-":
							dir_query="R"
						else:
							print("This row: "+line_split.strip("\n")+" has strand problem.")
							continue
					else:
						print("This row: "+line_split.strip("\n")+" has been merged incorrectly.")
						continue
			else:
				print("This row: "+line_split.strip("\n")+" has strand problem.")
			
			if (dir_other=="F" or dir_other=="R") and (dir_query=="F" or dir_query=="R"):
				if gene_id==gene_id_pre:
					if int(loc_other)-int(loc_pri)>=0:
						if int(loc_other)-int(loc_pri)<group_range:
							dic_name = loc_pri+"_"+dir_other+"_"+dir_query	
							if dic_name in dic_out[gene_id_pre].keys():
								dic_out[gene_id_pre][dic_name]+=1
								dic_detail[gene_id_pre][dic_name]=dic_detail[gene_id_pre][dic_name]+line_split
							else:
								dic_out[gene_id_pre][dic_name]=1
								dic_detail[gene_id_pre][dic_name]=line_split
						else:
							loc_pri=loc_other
							dic_name = loc_pri+"_"+dir_other+"_"+dir_query
							dic_out[gene_id_pre][dic_name]=1
							dic_detail[gene_id_pre][dic_name]=line_split
					else:
						print("this file is not correctly sorted!!!!!")
						break
				else:
                                        if gene_id!=gene_query_pre:
                                                #summarize the previous group with same gene_id
                                                for k in dic_out[gene_id_pre].keys():
                                                        dir_other_pre=k.split("_")[1]
                                                        if dir_other_pre=="F":
                                                                loc_other_col=2
                                                                loc_other_up=0
                                                        else:
                                                                loc_other_col=1
                                                                loc_other_up=10000000000000000000000000000
                                                        dir_query_pre=k.split("_")[2]
                                                        if dir_query_pre=="F":
                                                                loc_query_col=20
                                                                loc_query_up=0
                                                        else:
                                                                loc_query_col=19
                                                                loc_query_up=10000000000000000000000000000
                                                        reads=dic_detail[gene_id_pre][k].split("\n")
                                                        del reads[len(reads)-1]
                                                        for r in reads:
                                                                fields=r.split("\t")
                                                                loc_other_now=int(fields[loc_other_col])
                                                                loc_query_now=int(fields[loc_query_col])
                                                                gene_name_other=fields[15]
                                                                gene_name_query=fields[33]
                                                                chr_other=fields[0]
                                                                chr_query=fields[18]
                                                                gene_query=fields[34]
                                                                if gene_query!=gene_query_pre:
                                                                        print("This file has problem at query gene field, should double check!!!!")
                                                                        break
                                                                if dir_other_pre=="F":
                                                                        if loc_other_now>loc_other_up:
                                                                                loc_other_up=loc_other_now
                                                                else:
                                                                        if loc_other_now<loc_other_up:
                                                                                loc_other_up=loc_other_now
                                                                if dir_query_pre=="F":
                                                                        if loc_query_now>loc_query_up:
                                                                                loc_query_up=loc_query_now
                                                                else:
                                                                        if loc_query_now<loc_query_up:
                                                                                loc_query_up=loc_query_now
                                                        dic_name_update=gene_name_other+"_"+gene_id_pre+"_"+chr_other+"_"+str(loc_other_up)+"_"+dir_other_pre+"_"+gene_name_query+"_"+gene_query+"_"+chr_query+"_"+str(loc_query_up)+"_"+dir_query_pre
                                                        dic_out_final[dic_name_update]=dic_out[gene_id_pre][k]
                                                        dic_detail_final[dic_name_update]=dic_detail[gene_id_pre][k]
                                                
                                                #start the info of a new group
                                                gene_id_pre=gene_id
                                                loc_pri=loc_other
                                                dic_name = loc_pri+"_"+dir_other+"_"+dir_query	
                                                dic_out[gene_id_pre]={}
                                                dic_detail[gene_id_pre]={}
                                                dic_name = loc_pri+"_"+dir_other+"_"+dir_query
                                                dic_out[gene_id_pre][dic_name]=1
                                                dic_detail[gene_id_pre][dic_name]=line_split
			else:
				print("Can not extract directional information from this row: "+line_split.strip("\n"))
	
	
	#summarize the last group
	gene_id_pre=gene_id
        if gene_id!=gene_query_pre:
                for k in dic_out[gene_id_pre].keys():
                        dir_other_pre=k.split("_")[1]
                        if dir_other_pre=="F":
                                loc_other_col=2
                                loc_other_up=0
                        else:
                                loc_other_col=1
                                loc_other_up=10000000000000000000000000000
                        dir_query_pre=k.split("_")[2]
                        if dir_query_pre=="F":
                                loc_query_col=20
                                loc_query_up=0
                        else:
                                loc_query_col=19
                                loc_query_up=10000000000000000000000000000
                        reads=dic_detail[gene_id_pre][k].split("\n")
                        del reads[len(reads)-1]
                        for r in reads:
                                fields=r.split("\t")
                                loc_other_now=int(fields[loc_other_col])
                                loc_query_now=int(fields[loc_query_col])
                                gene_name_other=fields[15]
                                gene_name_query=fields[33]
                                chr_other=fields[0]
                                chr_query=fields[18]
                                gene_query=fields[34]
                                if gene_query!=gene_query_pre:
                                        print("This file has problem at query gene field, should double check!!!!")
                                        break
                                if dir_other_pre=="F":
                                        if loc_other_now>loc_other_up:
                                                loc_other_up=loc_other_now
                                else:
                                        if loc_other_now<loc_other_up:
                                                loc_other_up=loc_other_now
                                if dir_query_pre=="F":
                                        if loc_query_now>loc_query_up:
                                                loc_query_up=loc_query_now
                                else:
                                        if loc_query_now<loc_query_up:
                                                loc_query_up=loc_query_now
                        dic_name_update=gene_name_other+"_"+gene_id_pre+"_"+chr_other+"_"+str(loc_other_up)+"_"+dir_other_pre+"_"+gene_name_query+"_"+gene_query+"_"+chr_query+"_"+str(loc_query_up)+"_"+dir_query_pre
                        dic_out_final[dic_name_update]=dic_out[gene_id_pre][k]
                        dic_detail_final[dic_name_update]=dic_detail[gene_id_pre][k]
		
	for key_1 in dic_out_final.keys():
		fusion_location='\t'.join(key_1.split("_"))+"\t"+str(dic_out_final[key_1])+"\n"	
		h_out.write(fusion_location)
		h_out_detail.write(">>>"+key_1+"\t"+str(dic_out_final[key_1])+"\n"+dic_detail_final[key_1])
		if dic_out_final[key_1]>=int(filter_num):
			h_out_filter.write(fusion_location)
			h_out_detail_filter.write(">>>"+key_1+"\t"+str(dic_out_final[key_1])+"\n"+dic_detail_final[key_1])
							
			
	
	h_span.close()
	h_out.close()
	h_out_detail.close()
	h_out_filter.close()
	h_out_detail_filter.close()