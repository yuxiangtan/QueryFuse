use Bio::DB::Sam;
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

#check if there are 3 parameters
die "Usage: perl split_bam_file.pl FILE_DIR BAM_FILE PROCESS_DIR" if @ARGV < 3;

my $FILEDIR=$ARGV[0];
my $BAMFILE=$ARGV[1];
my $PROCESSDIR=$ARGV[2];

#make a processing directory
if (!(-d $PROCESSDIR)){
   print 'making'.$PROCESSDIR.'... ';
   unless(mkdir $PROCESSDIR) {
   		die "Unable to create $PROCESSDIR\n";
   }
}

#open bam file
my $bam    = Bio::DB::Bam->open($FILEDIR.$BAMFILE);
my $header = $bam->header;

#open output files and write headers
my $singleton = Bio::DB::Bam->open($PROCESSDIR.'singleton.bam',"w");
my $paired = Bio::DB::Bam->open($PROCESSDIR.'paired.bam',"w");
my $unmapped = Bio::DB::Bam->open($PROCESSDIR.'unmapped.bam',"w");
#to save reads that have secondary alignment location.
my $secondary = Bio::DB::Bam->open($PROCESSDIR.'sec_alignment.bam',"w");

my $status_code = $singleton ->header_write($header);
$status_code = $paired ->header_write($header);
$status_code = $unmapped  ->header_write($header);
$status_code = $secondary  ->header_write($header);

my $bytes;
my $flag;
my @counter= (0,0,0);

#read line by line and output in 3 different files
while (my $align = $bam->read1) {
   $flag=$align->flag;
   #check if fourth flag is set (fragment is unmapped)
   if ($flag & 256){
      #0x256 ==1 secondary alignment
           $bytes = $secondary ->write1($align);
           $counter[3]++;
   } else{
      if ($flag & 4){
         #check if eight flag is set (next fragment is unmapped)
         if ($flag & 8) {
           #0x4 ==1 & 0x8 == 1
           $bytes = $unmapped ->write1($align);
           $counter[2]++;
         }else{
           #0x4 ==1 & 0x8 == 0
           $bytes = $singleton ->write1($align);
           $counter[1]++;
         }
      }else{
        if ($flag & 8) {
           #0x4 == 0 & 0x8 == 1
           $bytes = $singleton ->write1($align);
           $counter[1]++;
         }else{
           #0x4 == 0 & 0x8 == 0
           $bytes = $paired ->write1($align);
           $counter[0]++;
         }
      }
   }
}
#write stats into file
my $out=$PROCESSDIR.'stats.txt';
open(OUTFILE, ">$out") or die "Cannot write $out: $!.\n";
print OUTFILE 'Singletons: ', $counter[1],' ';
print OUTFILE 'Paired: ',$counter[0],' ';
print OUTFILE 'Unmapped: ',$counter[2],' ';
print OUTFILE 'Secondary alignment: ',$counter[3],' ';
close OUTFILE
