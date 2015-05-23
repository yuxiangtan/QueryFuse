use strict;

#Transform the format into same as split mate to query
#Filter out reads with not full matched reads.

#scan line by line

#check if there are 2 parameters
die "Warning: Usage: perl split_mate_to_other_triple_summary_transform.pl split_mate_to_other_blat.psl_summary_sorted.bed_merge_with_mate_summary.bed read_lenth Align_percent:def_98" if @ARGV < 3;

my $PSL_SUM=$ARGV[0];
my $READ_LEN=$ARGV[1];
my $ALIGN_PER=$ARGV[2];
my $OUT_BED=$PSL_SUM."_transform.bed";

open(FILEHANDLE1,$PSL_SUM)  || die("Can't open the file£º$!");

open(OUTFILE, ">$OUT_BED") or die "Cannot write $OUT_BED: $!.\n";

my $count=1;

while (<FILEHANDLE1>) {
	my $keyword_line=$_;
	$keyword_line=~ s/\n//g;
	if ($count==1){
		#print $keyword_line."\n";
		#print $CHR."\n";
		$count=2;
	}
	my @LINESPLIT1=split(/\t/,$keyword_line); 
	my $readID1=@LINESPLIT1[9];
	my $ensembl1=@LINESPLIT1[13];
	my $match1=@LINESPLIT1[0];
	my $gene_name=@LINESPLIT1[57];
	my $CHR=@LINESPLIT1[54];
	my $strand=@LINESPLIT1[8];
	my $begin_loc=@LINESPLIT1[55];
	my $end_loc=@LINESPLIT1[56];
	my $range1=@LINESPLIT1[15];
	my $range2=@LINESPLIT1[16];
	my $match2=@LINESPLIT1[21];
	my $rangeleft;
	my $rangeright;
	if ($range1<=$range2){
		$rangeleft=$range1;
		$rangeright=$range2;
	} else {
		$rangeleft=$range2;
		$rangeright=$range1;
	}
	my $RANGE_LEFT=$begin_loc+$rangeleft;
	my $RANGE_RIGHT=$begin_loc+$rangeright;
	if (($match1+$match2)>($READ_LEN*$ALIGN_PER/100)){
		if (($match1+$match2)<=$READ_LEN){
			print OUTFILE $CHR."\t".$RANGE_LEFT."\t".$RANGE_RIGHT."\t".$readID1."\t"."NA"."\t".$strand."\t".$RANGE_LEFT."\t".$RANGE_RIGHT."\t"."0,0,0"."\t"."1"."\t".$match1.","."\t"."0,"."\t".$CHR."\t".$begin_loc."\t".$end_loc."\t".$gene_name."\t".$ensembl1."\t".$match1."\t".@LINESPLIT1[21]."\t".@LINESPLIT1[22]."\t".@LINESPLIT1[23]."\t".@LINESPLIT1[24]."\t".@LINESPLIT1[25]."\t".@LINESPLIT1[26]."\t".@LINESPLIT1[27]."\t".@LINESPLIT1[28]."\t".@LINESPLIT1[29]."\t".@LINESPLIT1[30]."\t".@LINESPLIT1[31]."\t".@LINESPLIT1[32]."\t".@LINESPLIT1[33]."\t".@LINESPLIT1[34]."\t".@LINESPLIT1[35]."\t".@LINESPLIT1[36]."\t".@LINESPLIT1[37]."\t".@LINESPLIT1[38]."\t".@LINESPLIT1[39]."\t".@LINESPLIT1[40]."\t".@LINESPLIT1[41]."\t".@LINESPLIT1[42]."\t".@LINESPLIT1[43]."\t".@LINESPLIT1[44]."\t".@LINESPLIT1[45]."\t".@LINESPLIT1[46]."\t".@LINESPLIT1[47]."\t".@LINESPLIT1[48]."\t".@LINESPLIT1[49]."\t".@LINESPLIT1[50]."\t".@LINESPLIT1[51]."\t".@LINESPLIT1[52]."\t".@LINESPLIT1[53]."\t".@LINESPLIT1[54]."\t".@LINESPLIT1[55]."\t".@LINESPLIT1[56]."\t".@LINESPLIT1[57]."\t".@LINESPLIT1[58]."\t".@LINESPLIT1[59]."\n";
		}
	}
}	
	
	
	
	