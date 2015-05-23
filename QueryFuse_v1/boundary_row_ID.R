#boundary_row_ID.R
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
	stop("\n\nWarning: file.in is not exist for boundary_row_ID.R, please check the path, exit. \n\n")
}


#read in the file
#since it can only have a column of numbers, read.table is OK in this script
match_frame<-read.table(file=file.in)

lower_b<-which(match_frame<lower_bound)
upper_b<-which(match_frame>upper_bound)
#if  (length(upper_b)==0) {
#	upper_b=dim(match_frame)[1]
#}

out_matr<-c(0,0)
out_matr[1]<-length(lower_b)
out_matr[2]<-upper_b[1]

#output the file
write.table(out_matr,file=`file.out`,quote=FALSE,row.name=FALSE,col.name=FALSE)

