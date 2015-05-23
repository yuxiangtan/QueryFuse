use strict;
#extract .fa file by read ID
#read readID into hash table
#scan through the .fa file line by line
#need to know interval (how many rows betweent two read ID) first, can get by head and sed or grep.

#check if there are 3 parameters
die "Usage: perl extract_fa_by_readID.pl fa_file readID_list row_interval out_fa" if @ARGV < 4;

my $FAFILE=$ARGV[0];
my $READID=$ARGV[1];
my $ROWINT=$ARGV[2];
my $OUTFA=$ARGV[3];

open(FILEHANDLE1,$FAFILE)  || die("Can't open the file£º$!");

open(FILEHANDLE2,$READID)  || die("Can't open the file£º$!");
my %ID_LIST;
while (<FILEHANDLE2>) {
	my $keyword_line=$_;
	$ID_LIST{$keyword_line}=0;
}



open(OUTFILE, ">$OUTFA") or die "Cannot write $OUTFA: $!.\n";

while (<FILEHANDLE1>) {
	my $keyword_line=$_;
	if (substr($keyword_line,0,1)eq">" ) {
		my $keyword_line1=$keyword_line;
		$keyword_line1=~s/^.//;
		$keyword_line1=~s/^\s+//;
		
		#in order to get rid of potential extra /1 or /2 end mark in the ID
		$keyword_line1=~s/\/1$//;
		$keyword_line1=~s/\/2$//;
		
		if(exists$ID_LIST{$keyword_line1}){
			print OUTFILE "$keyword_line";
		       	for (my $count = $ROWINT; $count >= 1; $count--) {
			 	my $line=<FILEHANDLE1>;
			 	print OUTFILE "$line";
			 	}
		}else{
			for (my $count = $ROWINT; $count >= 1; $count--) {
				my $line=<FILEHANDLE1>;
			}
		}
	}else{
		die "Fail: READ ID is not begin with >, or row interval number is wrong \n";
	}
}