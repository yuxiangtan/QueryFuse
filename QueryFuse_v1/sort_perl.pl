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

