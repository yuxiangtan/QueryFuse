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

#look at the chr and coordinate of ends in a pair read, if all the same to the existing one, then consider it duplicate and eliminate it.
#the input file must be sorted? [harsh talbe should be more robust]
#or I can just use harsh table to save? 
#use dictionary in python, so it should be dic in dic, because of 3 sectors.

#check num of parameter
if len(sys.argv) < 3:
	print "Warning! Usage: python eliminate_dup_for_summary_merged.py bed_file outputfile"
	sys.exit(1)
    
datafile = sys.argv[1]
outfile_name = sys.argv[2]
#L = int(sys.argv[3])

#handle=open(datafile,"r",encoding='utf-8')
#outfile=open(outfile_name,"w",encoding='utf-8')
#encoding is not working for python 2.x
handle=open(datafile,"r")
outfile=open(outfile_name,"w")

End_other={}
#End_query={}
for line in handle:
	items=line.strip().split("\t")
	sec1=items[0]+','+items[1]+','+items[2]
	sec2=items[33]+','+items[34]
	sec3=items[39]+','+items[40]+','+items[41]
	#check whether sec1 exist
	if sec1 not in End_other.keys():
		#each step create a new dictionary, if not, can not fill content in for the sub-levels
		End_other[sec1]={}
		End_other[sec1][sec2]={}
		End_other[sec1][sec2][sec3]=1
		outfile.write(line)
	else:
		#check whether sec2 exist in dic of sec1
		if sec2 not in End_other[sec1].keys():
			End_other[sec1][sec2]={}
			End_other[sec1][sec2][sec3]=1
			outfile.write(line)
		else:
			#check whether all 3 secs are continously together, if yes add 1 count.
			if sec3 not in End_other[sec1][sec2].keys() :
				End_other[sec1][sec2][sec3]=1
				outfile.write(line)
			else:
				End_other[sec1][sec2][sec3]+=1

handle.close()
outfile.close()	

