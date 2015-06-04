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

#check the ID to see whether it is totally same, if yes, generate two lines
#if not, then pairs, put the one to other at the front and to query at the end
#if the readID is totally the same but not on query, means it aligned at the overlapped part of two genes, so record the big range and combine two ensemble together using "+" into one.

#check if there are 5 parameters
die "Usage: perl QF_pair_A1.1_group_to_pair.pl pair_anno_file out_file query.bed read_length alignment_percentage" if @ARGV < 5;

my $BED=$ARGV[0];
my $OUT=$ARGV[1];
my $QUERY_FILE=$ARGV[2];
my $READ_LEN=$ARGV[3];
my $Align_percent=$ARGV[4];
#my $COOR_KEY=$ARGV[4];
#my $QUERY_FILE=$ARGV[5];
open(FILEHANDLE1,$BED)  || die("Can't open the file£º$!");
open(FILEHANDLE2,$QUERY_FILE)  || die("Can't open the file£º$!");
open(OUTFILE, ">$OUT") or die "Cannot write $OUT: $!.\n";

my $keyword_line2=<FILEHANDLE2>;
$keyword_line2=~ s/\n//g;
my @LINESPLIT2=split(/\t/,$keyword_line2); 
my $QUERY=@LINESPLIT2[4];
my $COR1;
my $COR2;
my $COOR;
my $ENSEM;
my $keyword_line=<FILEHANDLE1>;
$keyword_line=~ s/\n//g;
my @LINESPLIT1=split(/\t/,$keyword_line); 
my $readID1=@LINESPLIT1[3];
my $pre_readID=$readID1;
my $pre_line=$keyword_line;
my $ensembl1=@LINESPLIT1[16];
my $match1=@LINESPLIT1[17];
my @pre_LINESPLIT1;
my $pre_readID1;
my $pre_ensembl1;
my $pre_COR1;
my $pre_COR2;
my $pre_match1;
my $PART1;
my $PART2;
my $comb_COR1;
my $comb_COR2;
my $com_ensembl;
my $com_geneID;
my $pre_gene_ID;
my $gene_ID;
my @readID_st1;
my @readID_st2;
my $readID2;
my $pre_ensembl_COR2;
my $pre_ensembl_COR1;
my $ensembl_COR2;
my $ensembl_COR1;
	while (<FILEHANDLE1>) {
		$keyword_line=$_;
		$keyword_line=~ s/\n//g;
		@LINESPLIT1=split(/\t/,$keyword_line);
		$readID1=@LINESPLIT1[3];
		$ensembl1=@LINESPLIT1[16];
		$COR1=@LINESPLIT1[1];
		$COR2=@LINESPLIT1[2];
		$match1=@LINESPLIT1[17];
		if ($match1>=($READ_LEN*$Align_percent/100)){
			#print $readID1."\n";
			#print $pre_readID."\n";
			if (($readID1)ne($pre_readID)){
				@readID_st1=split(/\//,$readID1);
				#print @readID_st."\n";
				$readID1=$readID_st1[0];
				@readID_st2=split(/\//,$pre_readID);
				#print @readID_st."\n";
				$readID2=$readID_st2[0];
				
				if (($readID1)eq($readID2)){
					@pre_LINESPLIT1=split(/\t/,$pre_line);
					$pre_readID1=@pre_LINESPLIT1[3];
					$pre_ensembl1=@pre_LINESPLIT1[16];
					$pre_COR1=@pre_LINESPLIT1[1];
					$pre_COR2=@pre_LINESPLIT1[2];
					$pre_match1=@pre_LINESPLIT1[17];
					if (($ensembl1)eq($QUERY)){
						if ($COR1<$COR2){
							$PART2=$keyword_line;
						} else {
							$PART2=@LINESPLIT1[0]."\t".@LINESPLIT1[2]."\t".@LINESPLIT1[1]."\t".@LINESPLIT1[3]."\t".@LINESPLIT1[4]."\t".@LINESPLIT1[5]."\t".@LINESPLIT1[6]."\t".@LINESPLIT1[7]."\t".@LINESPLIT1[8]."\t".@LINESPLIT1[9]."\t".@LINESPLIT1[10]."\t".@LINESPLIT1[11]."\t".@LINESPLIT1[12]."\t".@LINESPLIT1[13]."\t".@LINESPLIT1[14]."\t".@LINESPLIT1[15]."\t".@LINESPLIT1[16]."\t".@LINESPLIT1[17];
						}
						if ($pre_COR1<$pre_COR2){
							$PART1=$pre_line;	
						} else {
							$PART1=@pre_LINESPLIT1[0]."\t".@pre_LINESPLIT1[2]."\t".@pre_LINESPLIT1[1]."\t".@pre_LINESPLIT1[3]."\t".@pre_LINESPLIT1[4]."\t".@pre_LINESPLIT1[5]."\t".@pre_LINESPLIT1[6]."\t".@pre_LINESPLIT1[7]."\t".@pre_LINESPLIT1[8]."\t".@pre_LINESPLIT1[9]."\t".@pre_LINESPLIT1[10]."\t".@pre_LINESPLIT1[11]."\t".@pre_LINESPLIT1[12]."\t".@pre_LINESPLIT1[13]."\t".@pre_LINESPLIT1[14]."\t".@pre_LINESPLIT1[15]."\t".@pre_LINESPLIT1[16]."\t".@pre_LINESPLIT1[17];
						}
					} else {
						if (($pre_ensembl1)eq($QUERY)){
							if ($COR1<$COR2){
								$PART1=$keyword_line;
							} else {
								$PART1=@LINESPLIT1[0]."\t".@LINESPLIT1[2]."\t".@LINESPLIT1[1]."\t".@LINESPLIT1[3]."\t".@LINESPLIT1[4]."\t".@LINESPLIT1[5]."\t".@LINESPLIT1[6]."\t".@LINESPLIT1[7]."\t".@LINESPLIT1[8]."\t".@LINESPLIT1[9]."\t".@LINESPLIT1[10]."\t".@LINESPLIT1[11]."\t".@LINESPLIT1[12]."\t".@LINESPLIT1[13]."\t".@LINESPLIT1[14]."\t".@LINESPLIT1[15]."\t".@LINESPLIT1[16]."\t".@LINESPLIT1[17];
							}
							if ($pre_COR1<$pre_COR2){
								$PART2=$pre_line;	
							} else {
								$PART2=@pre_LINESPLIT1[0]."\t".@pre_LINESPLIT1[2]."\t".@pre_LINESPLIT1[1]."\t".@pre_LINESPLIT1[3]."\t".@pre_LINESPLIT1[4]."\t".@pre_LINESPLIT1[5]."\t".@pre_LINESPLIT1[6]."\t".@pre_LINESPLIT1[7]."\t".@pre_LINESPLIT1[8]."\t".@pre_LINESPLIT1[9]."\t".@pre_LINESPLIT1[10]."\t".@pre_LINESPLIT1[11]."\t".@pre_LINESPLIT1[12]."\t".@pre_LINESPLIT1[13]."\t".@pre_LINESPLIT1[14]."\t".@pre_LINESPLIT1[15]."\t".@pre_LINESPLIT1[16]."\t".@pre_LINESPLIT1[17];
							}
						}
					}
				
					
					
				} else {
					print OUTFILE $PART1."\t".$PART2."\n";
					$pre_line=$keyword_line;
					$pre_readID=@LINESPLIT1[3];
				}
			} else {
				#print $ensembl1;
				if (($ensembl1)eq($QUERY)){
					$pre_line=$keyword_line;
					$pre_readID=$readID1;
				} else {
					@pre_LINESPLIT1=split(/\t/,$pre_line);
					$pre_readID1=@pre_LINESPLIT1[3];
					$pre_ensembl1=@pre_LINESPLIT1[16];
					$pre_COR1=@pre_LINESPLIT1[1];
					$pre_COR2=@pre_LINESPLIT1[2];
					$pre_match1=@pre_LINESPLIT1[17];
					#print $pre_ensembl1;
					#print $ensembl1;
					if (($pre_ensembl1)ne($QUERY)){
						$pre_gene_ID=@pre_LINESPLIT1[15];
						$pre_ensembl_COR2=@pre_LINESPLIT1[14];
						$pre_ensembl_COR1=@pre_LINESPLIT1[13];
						$gene_ID=@LINESPLIT1[15];
						$ensembl_COR2=@LINESPLIT1[14];
						$ensembl_COR1=@LINESPLIT1[13];
						if ($pre_ensembl_COR1<$ensembl_COR1){
							$comb_COR1=$pre_ensembl_COR1;
						} else {
							$comb_COR1=$ensembl_COR1;
						}
						if ($pre_ensembl_COR2>$ensembl_COR2){
							$comb_COR2=$pre_ensembl_COR2;
						} else {
							$comb_COR2=$ensembl_COR2;
						}
						$com_ensembl=$pre_ensembl1."+".$ensembl1;
						$com_geneID=$pre_gene_ID."+".$gene_ID;
						$pre_line=@pre_LINESPLIT1[0]."\t".@pre_LINESPLIT1[2]."\t".@pre_LINESPLIT1[1]."\t".@pre_LINESPLIT1[3]."\t".@pre_LINESPLIT1[4]."\t".@pre_LINESPLIT1[5]."\t".@pre_LINESPLIT1[6]."\t".@pre_LINESPLIT1[7]."\t".@pre_LINESPLIT1[8]."\t".@pre_LINESPLIT1[9]."\t".@pre_LINESPLIT1[10]."\t".@pre_LINESPLIT1[11]."\t".@pre_LINESPLIT1[12]."\t".@pre_LINESPLIT1[13]."\t".@pre_LINESPLIT1[14]."\t".$com_geneID."\t".$com_ensembl."\t".@pre_LINESPLIT1[17];
							
					}
				}
			}
		} 
	}

print OUTFILE $PART1."\t".$PART2."\n";
