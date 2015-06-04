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