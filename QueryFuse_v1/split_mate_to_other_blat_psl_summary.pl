use strict;

#work on split_mate_to_other_blat.psl
#for each read, get the top1 matched as true alignment, if more than 1, then keep both.
#The file is sorted by readID already, so just read line by line.

#check if there are 1 parameters
die "Usage: perl split_mate_to_other_blat_psl_summary.pl split_mate_to_other_blat.psl" if @ARGV < 1;

my $PSL=$ARGV[0];

my $OUT_BED=$PSL."_summary.bed";

open(FILEHANDLE1,$PSL)  || die("Can't open the file£º$!");
my $keyword_line=<FILEHANDLE1>;
my @LINESPLIT1=split(/\t/,$keyword_line); 
my $pre_readID=@LINESPLIT1[9];
my $max_match_num=@LINESPLIT1[0];
my $match_max_row=$keyword_line;

open(OUTFILE, ">$OUT_BED") or die "Cannot write $OUT_BED: $!.\n";

while (<FILEHANDLE1>) {
	my $keyword_line=$_;
	my @LINESPLIT1=split(/\t/,$keyword_line);
	my $readID=@LINESPLIT1[9];
	my $match_num=@LINESPLIT1[0];
	if (($readID)eq($pre_readID)){
		if ($match_num>$max_match_num){
				$match_max_row=$keyword_line;
				$max_match_num=$match_num;
			} else {
				if ($match_num==$max_match_num){
					$match_max_row=$match_max_row.$keyword_line;
				}
			}
	} else {
			print OUTFILE $match_max_row;
			$pre_readID=$readID;
			$max_match_num=$match_num;
			$match_max_row=$keyword_line;
	}
}

print OUTFILE $match_max_row;