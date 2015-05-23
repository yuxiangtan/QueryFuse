#/share/apps/R/R-2.12.2_gnu412_x86_64_vanilla/bin/Rscript KS_test_for_split_pval.R in_string="matrix of barplot table"

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

in_dist=as.numeric(strsplit(in_string,"_")[[1]])
len_dist=length(in_dist)
temp_dist=rep(mean(in_dist)/len_dist,len_dist)

out_val=ks.test(in_dist,temp_dist)$p.value

print(out_val)