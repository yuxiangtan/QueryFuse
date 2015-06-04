##Generally use this R version: /share/apps/R/R-2.12.2_gnu412_x86_64_vanilla/bin/Rscript --vanilla 

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


#this script is to find the intersect of two gene ID list
#use linux grep command is extremely slow
#first, read in the file and turn into matrix
#second, get the intersect
#third, output in to table
#intersect_gene_list.R file.list1 file.list2 file.out

#check arguments
for (e in commandArgs()) {
        ta = strsplit(e,"=",fixed=TRUE)
        if(! is.na(ta[[1]][2])) {
                temp = ta[[1]][2]
                if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
                temp = as.integer(temp)
                }
        if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
                temp = as.numeric(temp)
                }
        assign(ta[[1]][1],temp)
        } else {
        assign(ta[[1]][1],TRUE)
        }
}

#check whether file in is exist
if (!exists("file.list1")) {
	stop("\n\nWarning: Usage: gene_ID list1 is not exist, please check the path. \n\n")
}

#check whether file in is exist
if (!exists("file.list2")) {
	stop("\n\nWarning: Usage: gene_ID list 2 is not exist, please check the path. \n\n")
}



#read two list of ID
list1<-read.delim2(file=file.list1,header=FALSE,stringsAsFactors=FALSE)
list2<-read.delim2(file=file.list2,header=FALSE,stringsAsFactors=FALSE)

#check list format
list1_check=substr(list1[1,1],nchar(list1[1,1])-1,nchar(list1[1,1])-1)=="/"
list2_check=substr(list2[1,1],nchar(list2[1,1])-1,nchar(list2[1,1])-1)=="/"

if (list1_check!=list2_check){
    if (list1_check=="TRUE") {
        split_out<- function(x, c1) substr(x,1,nchar(x)-2)
        list1<-apply(list1,1,split_out)
        #list2<-as.matrix(list2)
    } else {
        split_out<- function(x, c1) substr(x,1,nchar(x)-2)
        list2<-apply(list2,1,split_out)
        #list1<-as.matrix(list1)
    }
}

list1<-as.matrix(list1)
list2<-as.matrix(list2)


#find the intersect
intersect_list<-intersect(list1,list2)

#output the file
write.table(intersect_list, file=`file.out`,  row.names = FALSE , col.names = FALSE, quote=FALSE)