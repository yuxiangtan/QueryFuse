use strict;

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

#use of hash table
#the split_to_mate_summary_table readID should be subset of ID in mate table.
#look at the summary_table,get an readID
#read in that ID from unmapped_to_query to get the other subset, use the with most match.
#read in that ID from mate and build hash table for each multi alignment, use whole alignment as content for that ensembleID. if have different location of a gene, then attach them together
#check the current read whether aligned to one of the mate, if yes,print this pair alignment
#read the summary table for the read multi alignment for this read.
#warning: three files must be sorted by readID

#check if there are 3 parameters
die "Warning: Usage: perl split_mate_to_other_grouped_into_triple_summary.pl split_mate_to_other_blat.psl_summary_sorted mate_singlet_to_other_only_sorted unmapped.bam_1/2_mate_on_query.psl_sorted" if @ARGV < 3;

my $PSL_SUM=$ARGV[0];
my $MATE_OTHER=$ARGV[1];
my $QUERY_PSL=$ARGV[2];

my $OUT_BED=$PSL_SUM."_merge_with_mate_summary.bed";

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
my $readID1=@LINESPLIT1[9];
#my $readID1_temp=@LINESPLIT1[9];
my $pre_readID=$readID1;
my $ensembl1=@LINESPLIT1[13];
my $keyword_line2;
my @LINESPLIT2;
my @readID_st;
my $readID2;
my $ensembl2;

#in order to get rid of potential extra /1 or /2 end mark in the ID (if they do not have, this work but also not fail.)
$readID1=~s/\/1$//;
$readID1=~s/\/2$//;

#my $split_key="on";
##check whether readID1's length changed by get rid of /1 or /2, if yes set split / on read2 off.
#if (length($readID1)!=length($readID1_temp) ) {
#	$split_key="off";
#}

		
#print "file1_line ".$keyword_line."\n";
#print "readID in file1: ".$readID1."\n";
#my $count2=0;

while (<FILEHANDLE2>) {
	$keyword_line2=$_;
	$keyword_line2=~ s/\n//g;
	@LINESPLIT2=split(/\t/,$keyword_line2);
	#check whether ID1 has /1 or /2 and determine whether to split.
	#if ($split_key eq "on") {
	#	@readID_st=split(/\//,@LINESPLIT2[3]);
	#	$readID2=$readID_st[0];
	#} else {
	#	$readID2=@LINESPLIT2[3];
	#}
	@readID_st=split(/\//,@LINESPLIT2[3]);
	$readID2=$readID_st[0];
	#print @readID_st."\n";
	#print $readID2."\n";
	$ensembl2=@LINESPLIT2[16];
	#if ( $count2==0 ){
	#	print "readID in file2: ".$readID2."\n";
	#}
	#$count2 +=1;

	if (($readID2)eq($readID1)){
		if (exists$ID_LIST{$ensembl2}){
			$ID_LIST{$ensembl2} .= "\t".$keyword_line2;
			#print $ID_LIST{$ensembl2}."\n";
		} else {
			$ID_LIST{$ensembl2}=$keyword_line2;
			#print $ID_LIST{$ensembl2}."\n";
		}
	} else {
		if (($readID2)gt($readID1)){
			#print $keyword_line."\n";
			#print "this ID bigger than existing"."\n";
			last;
		}
	}
}

my $keyword_line3;
my @LINESPLIT3;
my @readID_st3;
my $readID3;
my $match3;
my %MATCH_LIST;
my %READ_LIST;

#my $count3=0;

while (<FILEHANDLE3>) {
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
	
	#if ( $count3==0 ){
	#	print "readID in file3: ".$readID3."\n";
	#}
	#$count3 +=1;
	
	if (($readID3)eq($readID1)){
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
	} else {
		if (($readID3)gt($readID1)){
			#print $keyword_line."\n";
			#print "this ID bigger than existing"."\n";
			last;
		}
	}
}
if ((exists$ID_LIST{$ensembl1})&&(exists$READ_LIST{$readID1})){
	print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$ensembl1}."\n";
}



#read the following lines
while (<FILEHANDLE1>) {
	$keyword_line=$_;
	$keyword_line=~ s/\n//g;
	@LINESPLIT1=split(/\t/,$keyword_line); 
	$readID1=@LINESPLIT1[9];
	#get rid of the potential problem
	$readID1=~s/\/1$//;
	$readID1=~s/\/2$//;
	$ensembl1=@LINESPLIT1[13];
	if (($readID1)eq($pre_readID)){
		if ((exists$ID_LIST{$ensembl1})&&(exists$READ_LIST{$readID1})){
			print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$ensembl1}."\n";
		}
	} else {
		while (<FILEHANDLE3>) {
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
			if (($readID3)eq($readID1)){
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
			} else {
				if (($readID3)gt($readID1)){
					$MATCH_LIST{$readID3}=$match3;
					$READ_LIST{$readID3}=$keyword_line3;
					#print $keyword_line."\n";
					#print "this ID bigger than existing"."\n";
					last;
				}
			}
		}
		if (($readID2)lt($readID1)){
			%ID_LIST = ();
			while (<FILEHANDLE2>) {
				$keyword_line2=$_;
				$keyword_line2=~ s/\n//g;
				@LINESPLIT2=split(/\t/,$keyword_line2);
				#if ($split_key eq "on") {
				#	@readID_st=split(/\//,@LINESPLIT2[3]);
				#	$readID2=$readID_st[0];
				#} else {
				#	$readID2=@LINESPLIT2[3];
				#}
				@readID_st=split(/\//,@LINESPLIT2[3]);
				$readID2=$readID_st[0];
				$ensembl2=@LINESPLIT2[16];
				if (($readID2)eq($readID1)){
					if (exists$ID_LIST{$ensembl2}){
						$ID_LIST{$ensembl2}.= "\t".$keyword_line2;
						#print $ID_LIST{$ensembl2}."\n";
					} else {
						$ID_LIST{$ensembl2}=$keyword_line2;
						#print $ID_LIST{$ensembl2}."\n";
					}
				} else {
					if (($readID2)gt($readID1)){
						#print "this ID bigger than existing"."\n";
						last;
					}
				}
			}
			if ((exists$ID_LIST{$ensembl1})&&(exists$READ_LIST{$readID1})){
				print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$ensembl1}."\n";
			}
		} else {
			if (($readID2)eq($readID1)){
				%ID_LIST = ();
				$ID_LIST{$ensembl2}=$keyword_line2;
				while (<FILEHANDLE2>) {
					$keyword_line2=$_;
					$keyword_line2=~ s/\n//g;
					@LINESPLIT2=split(/\t/,$keyword_line2);
					#if ($split_key eq "on") {
					#	@readID_st=split(/\//,@LINESPLIT2[3]);
					#	$readID2=$readID_st[0];
					#} else {
					#	$readID2=@LINESPLIT2[3];
					#}
					@readID_st=split(/\//,@LINESPLIT2[3]);
					$readID2=$readID_st[0];

					$ensembl2=@LINESPLIT2[16];
					if (($readID2)eq($readID1)){
						if (exists$ID_LIST{$ensembl2}){
							$ID_LIST{$ensembl2} .= "\t".$keyword_line2;
							#print $ID_LIST{$ensembl2}."\n";
						} else {
							$ID_LIST{$ensembl2}=$keyword_line2;
							#print $ID_LIST{$ensembl2}."\n";
						}
					} else {
						if (($readID2)gt($readID1)){
							#print $keyword_line."\n";
							#print "this ID bigger than existing"."\n";
							last;
						}
					}
				}
				if ((exists$ID_LIST{$ensembl1})&&(exists$READ_LIST{$readID1})){
					print OUTFILE $keyword_line."\t".$READ_LIST{$readID1}."\t".$ID_LIST{$ensembl1}."\n";
				}
			} else {
				%ID_LIST = ();
			}
		}
	}
	$pre_readID=$readID1;
}