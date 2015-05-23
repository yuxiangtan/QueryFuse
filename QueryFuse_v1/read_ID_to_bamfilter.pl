use strict;

#transform ID list into bamtool filter script

#check if there are 3 parameters
die "Usage: perl read_ID_to_bamfilter.pl FILE_DIR read_ID PROCESS_DIR" if @ARGV < 3;

my $FILEDIR=$ARGV[0];
my $IDFILE=$ARGV[1];
my $PROCESSDIR=$ARGV[2];

#make a processing directory
if (!(-d $PROCESSDIR)){
   print 'making'.$PROCESSDIR.'... ';
   unless(mkdir $PROCESSDIR) {
   		die "Unable to create $PROCESSDIR\n";
   }
}


#write filter into file
my $out=$PROCESSDIR.$IDFILE.'filter.txt';
#print "$FILEDIR.$IDFILE";
open (IDLIST, $FILEDIR.$IDFILE) or die "can't open the file";
open(OUTFILE, ">$out") or die "Cannot write $out: $!.\n";
print OUTFILE '{"filters":[';
print OUTFILE "\n";

while( <IDLIST>)
{
    print OUTFILE '{"name":"'."$_".'"},'; 
}

print OUTFILE '{"name":"THIS_IS_NOT_A_READ_TOALIGN"}';
print OUTFILE "\n";
print OUTFILE ']}';
print OUTFILE "\n";

close OUTFILE
