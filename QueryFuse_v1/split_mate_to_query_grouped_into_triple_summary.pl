use strict;

#use of hash table
#the anno_sort readID should be subset of ID in mate table.
#look at the anno_sort,get an readID
#read in that ID from unmapped_to_query to get the other subset, use the one with most match.
#read in that ID from mate bed, get the one with most match and from query, for each multi alignment, use whole alignment as content for that ensembleID. if have different location of a gene, then attach them together
#check the current read whether 3 subs are there and first two has more than 99 matched, if yes,print this pair alignment
#warning: three files must be sorted by readID

#check if there are 5 parameters
die "Usage: perl split_mate_to_query_grouped_into_triple_summary.pl SORTED_ANNO_BED mate_singlet_to_query_only_sorted unmapped.bam_1/2_mate_on_query.psl_sorted Read_length Align_percent:def_98" if @ARGV < 5;

my $PSL_SUM=$ARGV[0];
my $MATE_OTHER=$ARGV[1];
my $QUERY_PSL=$ARGV[2];
my $READ_LEN=$ARGV[3];
my $ALIGN_PER=$ARGV[4];

my $OUT_BED=$PSL_SUM."_merge_with_mate_summary.bed";

#print $MATE_OTHER."\n";
#Open files
open(FILEHANDLE1,$PSL_SUM)  || die("Can't open the file£º$!");
open(FILEHANDLE2,$MATE_OTHER)  || die("Can't open the file£º$!");
open(FILEHANDLE3,$QUERY_PSL)  || die("Can't open the file£º$!");

open(OUTFILE, ">$OUT_BED") or die "Cannot write $OUT_BED: $!.\n";
my %ID_LIST;
#read first line
my $keyword_line=<FILEHANDLE1>;
$keyword_line=~ s/\n//g;
my @LINESPLIT1=split(/\t/,$keyword_line); 
my $readID1=@LINESPLIT1[3];
my $pre_readID=$readID1;
my $ensembl1=@LINESPLIT1[15];
my $match1=@LINESPLIT1[17];

#in order to get rid of potential extra /1 or /2 end mark in the ID (if they do not have, this work but also not fail.)
$readID1=~s/\/1$//;
$readID1=~s/\/2$//;

#print "read1 information\n";
#print $readID1."\n";
#print $ensembl1."\n";
#print $match1."\n"; 

my $keyword_line2;
my @LINESPLIT2;
my @readID_st;
my $readID2;
my $ensembl2;

my $keyword_line3;
my @LINESPLIT3;
my @readID_st3;
my $readID3;
my $match3;
my %MATCH_LIST;
my %READ_LIST;
my $query;
while (<FILEHANDLE3>) {
	$keyword_line3=$_;
	$keyword_line3=~ s/\n//g;
	@LINESPLIT3=split(/\t/,$keyword_line3);
	#@readID_st=split(/\//,@LINESPLIT2[9]);
	#print @readID_st."\n";
	$readID3=@LINESPLIT3[9];
	#print $readID3."\n";
	#get rid of the potential problem
	$readID3=~s/\/1$//;
	$readID3=~s/\/2$//;
	$match3=@LINESPLIT3[0];
	$query=@LINESPLIT3[13];
	#print $match3."\n";

	if (($readID3)eq($readID1)){

		if ($READ_LEN-$match3>=$match1){
			#print "read3 information, read3=read1\n";
			#print $readID3."\n";
			#print $query."\n";
			#print $match3."\n"; 
			if (exists$MATCH_LIST{$readID3}){
				if ($MATCH_LIST{$readID3}<$match3){
					$MATCH_LIST{$readID3}=$match3;
					$READ_LIST{$readID3}=$keyword_line3;
					#print "read3 information\n";
					#print $readID3."\n";
					#print $query."\n";
					#print $match3."\n"; 
				}
				#print $ID_LIST{$ensembl2}."\n";
			} else {
				$MATCH_LIST{$readID3}=$match3;
				$READ_LIST{$readID3}=$keyword_line3;
				#print "read3 information\n";
				#print $readID3."\n";
				#print $query."\n";
				#print $match3."\n"; 
				#print $READ_LIST{$readID3}."\n";
				#print $ID_LIST{$ensembl2}."\n";
			}
		} 
	}else {
			if (($readID3)gt($readID1)){
				#print $keyword_line3."\n";
				#print "this ID bigger than existing"."\n";
				last;
			}		
	}
}

#print $readID3."\n";
#print "\n";
#print "\n";
#print "\n";
#print "\n";

while (<FILEHANDLE2>) {
	$keyword_line2=$_;
	$keyword_line2=~ s/\n//g;
	@LINESPLIT2=split(/\t/,$keyword_line2);
	@readID_st=split(/\//,@LINESPLIT2[3]);
	#print $keyword_line2."\n";
	#print @LINESPLIT2."\n";
	#print @readID_st."\n";
	$readID2=$readID_st[0];
	#print $readID2."\n";
	$ensembl2=@LINESPLIT2[15];
	if (($readID2)eq($readID1)){
		#print "read2 information\n";
#print $readID2."\n";
#print $ensembl2."\n";
 
		if (($ensembl2)eq($query)){

			if (exists$ID_LIST{$ensembl2}){
				$ID_LIST{$ensembl2} .= "\t".$keyword_line2;
				#print $ID_LIST{$ensembl2}."\n";
			} else {
				$ID_LIST{$ensembl2}=$keyword_line2;
				#print $ID_LIST{$ensembl2}."\n";
			}
					#print "read2 information\n";
#print $readID2."\n";
#print $ensembl2."\n";
#print $ID_LIST{$ensembl2}."\n";		
		}
	} else {
		if (($readID2)gt($readID1)){
			#print $keyword_line."\n";
			#print "this ID bigger than existing"."\n";
			last;
		}
	}
	#print $ID_LIST."\n";
}

#print "condition\n";
#print $ID_LIST{$query}."\n";
#print $READ_LIST{$readID1}."\n";		
if ((exists$ID_LIST{$query})&&(exists$READ_LIST{$readID1})){
	#print "read1 information when both exist\n";
	#print $match1."\n";
	#print $MATCH_LIST{$readID1}."\n";
	#print $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$query}."\n";
	if (($match1+$MATCH_LIST{$readID1})>=($READ_LEN*$ALIGN_PER/100)){
			#print "yes!\n";
			print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$query}."\n";
	}
}



my $count=1;
#read the following lines
while (<FILEHANDLE1>) {
	#read in a line
	$keyword_line=$_;
	$keyword_line=~ s/\n//g;
	@LINESPLIT1=split(/\t/,$keyword_line); 
	$readID1=@LINESPLIT1[3];
	$ensembl1=@LINESPLIT1[15];
	$match1=@LINESPLIT1[17];
	$readID1=~s/\/1$//;
	$readID1=~s/\/2$//;
	#if ($count==1){
	#print "read1 information\n";
	#print $readID1."\n";
	#print $pre_readID."\n";
	#print $ensembl1."\n";
	#print $match1."\n"; 
	#$count=2;
	#}
	#check whether it is same ID as previous one
	if (($readID1)eq($pre_readID)){
		#if ($count==2){
		#print "read1 information\n";
		#print $ID_LIST{$query}."\n";
		#print $READ_LIST{$readID1}."\n";
		#$count=3;
		#}
				#even it is same ID, if it full fill all the criteria, consider it as multilocation.
				if ((exists$ID_LIST{$query})&&(exists$READ_LIST{$readID1})){
					#if ($count==3){
					#print "read1 information when both exist\n";
					#print $match1."\n";
					#print $MATCH_LIST{$readID1}."\n";
					#$count=4;
					#}
					if (($match1+$MATCH_LIST{$readID1})>=($READ_LEN*$ALIGN_PER/100)){
						#print "same as previous ID, but pass criterias!\n";
							print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$query}."\n";
					}
				}
	} else {
		#not the same as previous ID, if ID1 equal to ID3
		if (($readID3)eq($readID1)){
			#if ($count==4){
			#print "read3==readID1, loop1\n";
			#print "read3 information\n";
			#print $readID3."\n";
			#print $match3."\n";
			#print $query."\n";
			#print $match1."\n"; 
			#$count=5;
			#}
#				check whether match3 can be the biggest one. Do I need to sort match before sorting ID?
				if ($READ_LEN-$match3>=$match1){
					if (exists$MATCH_LIST{$readID3}){
						if ($MATCH_LIST{$readID3}<$match3){
							$MATCH_LIST{$readID3}=$match3;
							$READ_LIST{$readID3}=$keyword_line3;
						}
						#print $ID_LIST{$ensembl2}."\n";
					} else {
						$MATCH_LIST{$readID3}=$match3;
						$READ_LIST{$readID3}=$keyword_line3;
						#print $ID_LIST{$ensembl2}."\n";
					}
				}
			} else {
				#ID3 is not equal to ID1
					#if ID3 is bigger than ID1 than just quite the loop
					if (($readID3)gt($readID1)){
						#print $keyword_line."\n";
						#print "this ID bigger than existing"."\n";
						last;
					} else {
			#ID3 is smaller than ID1, the need to read forward to try to get the same ID
		#print "working\n".$readID3;
					while (<FILEHANDLE3>) {
						#print "working!!!!!!\n";
						$keyword_line3=$_;
						$keyword_line3=~ s/\n//g;
						@LINESPLIT3=split(/\t/,$keyword_line3);
						#@readID_st=split(/\//,@LINESPLIT2[9]);
						#print @readID_st."\n";
						$readID3=@LINESPLIT3[9];
						#get rid of the potential problem
						$readID3=~s/\/1$//;
						$readID3=~s/\/2$//;
						#print $readID2."\n";
						$match3=@LINESPLIT3[0];
						$query=@LINESPLIT3[13];
	
						#print "read3 information\n";
						#print $readID3."\n";
						#print $match3."\n";
						#print $query."\n";
						#print $match1."\n"; 
	
						#in the loop, check whether ID3 catch up ID1
						if (($readID3)eq($readID1)){
							#if ($count==4){
							#print "read3==readID1,loop2\n";
							#print "read3 information\n";
							#print $readID3."\n";
							#print $match3."\n";
							#print $query."\n";
							#print $match1."\n"; 
							#$count=5;
							#}
							#get max match3 in this case.
							if ($READ_LEN-$match3>=$match1){
								if (exists$MATCH_LIST{$readID3}){
									if ($MATCH_LIST{$readID3}<$match3){
										$MATCH_LIST{$readID3}=$match3;
										$READ_LIST{$readID3}=$keyword_line3;
									}
									#print $ID_LIST{$ensembl2}."\n";
								} else {
									$MATCH_LIST{$readID3}=$match3;
									$READ_LIST{$readID3}=$keyword_line3;
									#print $ID_LIST{$ensembl2}."\n";
								}
							}
						} else {
							#if ID3 past ID1, then quite the loop
							if (($readID3)gt($readID1)){
								#print $keyword_line."\n";
								#print "this ID bigger than existing"."\n";
								last;
							} 
						}
					}
				}
			}
		if (($readID2)lt($readID1)){
			%ID_LIST = ();
			while (<FILEHANDLE2>) {
				$keyword_line2=$_;
				$keyword_line2=~ s/\n//g;
				@LINESPLIT2=split(/\t/,$keyword_line2);
				@readID_st=split(/\//,@LINESPLIT2[3]);
				#print @readID_st."\n";
				$readID2=$readID_st[0];
				#print $readID2."\n";
				$ensembl2=@LINESPLIT2[15];
				if (($readID2)eq($readID1)){
					if (($ensembl2)eq($query)){
						if (exists$ID_LIST{$ensembl2}){
							$ID_LIST{$ensembl2} .= "\t".$keyword_line2;
							#print $ID_LIST{$ensembl2}."\n";
						} else {
							$ID_LIST{$ensembl2}=$keyword_line2;
							#print $ID_LIST{$ensembl2}."\n";
						}
					}
				} else {
					if (($readID2)gt($readID1)){
						#print $keyword_line."\n";
						#print "this ID bigger than existing"."\n";
						last;
					}
				}
			}

			if ((exists$ID_LIST{$query})&&(exists$READ_LIST{$readID1})){
				if (($match1+$MATCH_LIST{$readID1})>=($READ_LEN*$ALIGN_PER/100)){
					#print "yes3!\n";
						print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$query}."\n";
				}
			}

		} else {
			if (($readID2)eq($readID1)){
				%ID_LIST = ();
				$ID_LIST{$ensembl2}=$keyword_line2;
				while (<FILEHANDLE2>) {
					$keyword_line2=$_;
					$keyword_line2=~ s/\n//g;
					@LINESPLIT2=split(/\t/,$keyword_line2);
					@readID_st=split(/\//,@LINESPLIT2[3]);
					#print @readID_st."\n";
					$readID2=$readID_st[0];
					#print $readID2."\n";
					$ensembl2=@LINESPLIT2[15];
					if (($readID2)eq($readID1)){
						if (($ensembl2)eq($query)){
							if (exists$ID_LIST{$ensembl2}){
								$ID_LIST{$ensembl2} .= "\t".$keyword_line2;
								#print $ID_LIST{$ensembl2}."\n";
							} else {
								$ID_LIST{$ensembl2}=$keyword_line2;
								#print $ID_LIST{$ensembl2}."\n";
							}
						}
					} else {
						if (($readID2)gt($readID1)){
							#print $keyword_line."\n";
							#print "this ID bigger than existing"."\n";
							last;
						}
					}
				}
				if ((exists$ID_LIST{$query})&&(exists$READ_LIST{$readID1})){
					if (($match1+$MATCH_LIST{$readID1})>=($READ_LEN*$ALIGN_PER/100)){
						#print "yes4!\n";
							print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$query}."\n";
					}
				}
			} else {
				%ID_LIST = ();
			}
		}
	}
	$pre_readID=$readID1;
}
