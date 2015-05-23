#match_row_ID.R
#Get the row number of rows whose ID is in the ID list
#read in the bed_ID_list and read_ID_list
#which function to get the rows

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
if (!exists("file.bed")) {
	stop("\n\nWarning: file.bed is not exist for match_row_ID.R, please check the path, exit. \n\n")
}


if (!exists("file.list")) {
	stop("\n\nWarning: file.list is not exist for match_row_ID.R, please check the path, exit. \n\n")
}

#read in the file
bed_ID<-as.matrix(read.delim2(file=`file.bed`,header=FALSE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE))
list_ID<-as.matrix(read.delim2(file=`file.list`,header=FALSE,stringsAsFactors=FALSE,blank.lines.skip=FALSE,fill=TRUE))

row_ID<-which(bed_ID %in% list_ID)

#output the file
write.table(row_ID, file=`file.out`,  row.names = FALSE , col.names = FALSE, quote=FALSE)