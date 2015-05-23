use strict;

#test sorted in perl
#check if there are 1 parameters
die "Usage: perl sort_perl.pl file col_num > outfile" if @ARGV < 2;

my $FILE_IN=$ARGV[0];
my $COL=$ARGV[1];

open(FILEHANDLE1,$FILE_IN)  || die("Can't open the file£º$!");

my @by_uid;
#3 steps, first build the new array col1 with original line and col2 with target to sort col, 2nd sort by col2, 3rd, output only the col1.
@by_uid = map { $_->[0]} sort {$a->[1] cmp $b->[1] } map { [$_,(split (/\t/),$_)[$COL-1]] } <FILEHANDLE1>;

for(@by_uid){print ;}

