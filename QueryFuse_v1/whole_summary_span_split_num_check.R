#whole_summary_span_split_num_check.R file.in= file.out= HEAD_LINE=TRUEorFALSE frag_size= read_len= fragsize_std= pval_cut
#Get the row boundary of a sorted psl from lower and upper boundary of match number 

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
if (!exists("file.in")) {
	stop("\n\nUsage: file.in is not exist, please check the path. \n\n")
}

#test location and parameters:
#/restricted/projectnb/montilab-p/LinGA_unprotected/ytan/TL_PCNSL/align/TL_03_pairend_tophat_out/gene_specific_fusion_KIAA1432/SUMMARY
#file.in="./whole_fusion_sum_all.txt"
##not sure how to solve the \\t problem from entry, so just use default sep.
sep_char="\t"
#HEAD_LINE="FALSE"
#frag_size="150"
#read_len="99"
#fragsize_std="50"
#pval_cut=pnorm(2)

#need to make sure sep_char is string and HEAD_LINE is boolen
#sep_char=as.character(sep_char)
HEAD_LINE=as.logical(HEAD_LINE)
frag_size=as.numeric(frag_size)
read_len=as.numeric(read_len)
fragsize_std=as.numeric(fragsize_std)
pval_cut=as.numeric(pval_cut)

#read in the file after checking the file whether it is empty
if ((file.info(file.in)$size==0)|is.na(file.info(file.in)$size)){
    stop("\n\n file.in is not exist or empty. \n\n")
} else {
summary_table<-read.delim2(file=file.in, sep=sep_char, header = HEAD_LINE)
}

row_num=dim(summary_table)[1]
col_num=dim(summary_table)[2]

E_dis=frag_size/(read_len*2)-1
std_dis=fragsize_std/(read_len*2)

#count the pval of number of span/split reads following frag_size and read_len distribution.
summary_table[,col_num+1]=pnorm((as.numeric(summary_table[,col_num-1])/as.numeric(summary_table[,col_num-2])),E_dis,std_dis)
#use the pval_cut off to assign yes/no, then sort by yes/no
summary_table[,col_num+2]=(abs(summary_table[,col_num+1]-0.5)<=abs(pval_cut-0.5))

if (HEAD_LINE){
    colnames(summary_table)[col_num+1] <- "p-val of span/split read number"
    colnames(summary_table)[col_num+2] <- "pass of span/split read number check" 
}


#output the file
write.table(summary_table,file=`file.out`,quote=FALSE,row.name=FALSE,col.name=HEAD_LINE)

