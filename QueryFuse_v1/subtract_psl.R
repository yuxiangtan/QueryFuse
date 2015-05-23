#subtract_psl.R
#Modify the psl into its subtract
#First separate the psl into 5' aligned and 3' aligned group, and throw the rest away.
#modify the each group together.

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
	stop("\n\nWarning: file.in is not exist in subtract_psl.R, please check the path, exit. \n\n")
}


#read in the file
ori_psl<-read.csv(file=file.in,sep = "\t",quote="",header=FALSE)

ori_psl_5pi_1<- as.matrix(ori_psl[which(as.numeric(ori_psl[,12])<=5 ),])
if (dim(ori_psl_5pi_1)[2]>0){
if (dim(ori_psl_5pi_1)[2]<=1){
	ori_psl_5pi_1<- t(ori_psl_5pi_1)
}
}
ori_psl_5pi<- as.matrix(ori_psl_5pi_1[which(as.numeric(ori_psl_5pi_1[,13])<(as.numeric(read_length)-5)),])
if (dim(ori_psl_5pi)[2]>0){
if (dim(ori_psl_5pi)[2]<=1){
	ori_psl_5pi<- t(ori_psl_5pi)
}
}
#ori_psl_5pi[,12]<-as.numeric(ori_psl_5pi[,13])+1
#the reason not to plus 1 is because this file it to generate bed for fastafromBed, and it take 0 as the original start point.
ori_psl_5pi[,12]<-as.numeric(ori_psl_5pi[,13])
ori_psl_5pi[,13]<-read_length

#check whether the end is close to 3'
ori_psl_3pi_1<- as.matrix(ori_psl[which(as.numeric(ori_psl[,13])>=(as.numeric(read_length)-5)),])
if (dim(ori_psl_3pi_1)[2]>0){
if (dim(ori_psl_3pi_1)[2]<=1){
	ori_psl_3pi_1<- t(ori_psl_3pi_1)
}
}
ori_psl_3pi<-as.matrix(ori_psl_3pi_1[which(as.numeric(ori_psl_3pi_1[,12])>5),])
if (dim(ori_psl_3pi)[2]>0){
if (dim(ori_psl_3pi)[2]<=1){
	ori_psl_3pi<- t(ori_psl_3pi)
}
#ori_psl_3pi[,13]<-as.numeric(ori_psl_3pi[,12])-1
#this is for blat v. 34x13
ori_psl_3pi[,13]<-as.numeric(ori_psl_3pi[,12])
#if for online blat version, it should minus 1, because it will give the true correct start point 
ori_psl_3pi[,12]<-0
}
ori_psl_out<-rbind(ori_psl_5pi,ori_psl_3pi)
ori_bed_out<-as.matrix(ori_psl_out[,c(10,12,13,10)])
if (dim(ori_bed_out)[2]==1){
	ori_bed_out<- t(ori_bed_out)
}
#output the file
write.table(ori_bed_out,file=`file.out`,quote=FALSE,row.name=FALSE,col.name=FALSE,sep = "\t")

