#/share/apps/R/R-2.12.2_gnu412_x86_64_vanilla/bin/Rscript whole_summary_score_ranking.R infile="list of files to be merged" file.out= split_n=1 span_n=1 sum_n=2

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

#check whether files exist
if (!exists("infile")) {
	stop("\n\nUsage: infile is not exist, please check the path. \n\n")
}

HEADER_ROW=c("CHR_PARTNER_GENE","BREAKPOINT_PARTNER_GENE","DIRECTION_PARTNER_GENE","PARTNER_GENE_NAME","PARTNER_GENE_ID","CHR_QUERY_GENE","BREAKPOINT_QUERY_GENE","DIRECTION_QUERY_GENE","QUERY_GENE_NAME","QUERY_GENE_ID","SPLIT_NUM","SPAN_NUM","SUPPORT_SUM_NUM","SPLIT_PVAL","SHIFT_RANGE","DINUCLEOTIDE_ENTROPY","MULTI-ALIGN_NUM","RANK_SCORE")
colClasses_TYPE=c("character","character","character","character","character","character","character","character","character","character","numeric","integer","numeric","numeric","integer","numeric","integer")

file_in=read.delim(infile,header=FALSE,colClasses=colClasses_TYPE)
#colnames(file_in)=HEADER_ROW

#FOR ENTROPY, 16 cases, the max is -1/16*log(1/16)*16=2.77. The min is 0. -1/2*log(1/2)=0.35 means half cases are in one type, As a result, I use 0.7 as the cutoff.
KEEP_ID=which(file_in[,11]>=as.numeric(split_n) & file_in[,12]>=as.numeric(span_n) & file_in[,13]>=as.numeric(sum_n) & file_in[,16]>0.7)
RANK_SCORE_MATRIX=file_in[KEEP_ID,]

RANK_split_num=rank(-RANK_SCORE_MATRIX[,11])
RANK_span_num=rank(-RANK_SCORE_MATRIX[,12])
RANK_sum_num=rank(-RANK_SCORE_MATRIX[,13])
RANK_split_pval=rank(-RANK_SCORE_MATRIX[,14])
RANK_sphif_length=rank(-RANK_SCORE_MATRIX[,15])
RANK_entropy=rank(-RANK_SCORE_MATRIX[,16])
RANK_multi_align=rank(RANK_SCORE_MATRIX[,17])
RANK_ALL=RANK_split_num+RANK_span_num*2+RANK_sum_num+RANK_split_pval+RANK_sphif_length+RANK_entropy+RANK_multi_align

RANK_SCORE_MATRIX[,18]=RANK_ALL
RANK_SCORE_MATRIX_OUT=RANK_SCORE_MATRIX[order(RANK_SCORE_MATRIX[,18]),]
colnames(RANK_SCORE_MATRIX_OUT)=HEADER_ROW

write.table(RANK_SCORE_MATRIX_OUT, file=file.out, sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE)


