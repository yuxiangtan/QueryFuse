#!/bin/bash

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

#This script is just to pre-sort all the files that need to be sort by perl in advance.
#QF_single_end_perl_sort.sh

if [ $# -ne 5 ]
then
  echo ""
    echo "Warning: Usage: QF_single_end_perl_sort.sh UNMAP_ON_QUERY_PSL_1 UNMAP_ON_QUERY_PSL_2 SINGLE_ON_QUERY_BED SINGLE_ON_OTHER_BED QueryFuse_path perl_script"
    echo ""
    echo "UNMAP_ON_QUERY_PSL_1 - unmapped.bam_first_mate_on_query.psl"
    echo "UNMAP_ON_QUERY_PSL_2 - unmapped.bam_second_mate_on_query.psl"
    echo "SINGLE_ON_QUERY_BED - singleton_sorted_to_query.bed"
    echo "SINGLE_ON_OTHER_BED - singleton_sorted_to_other.bed"
    echo "QueryFuse_path - QueryFuse path."
    echo ""

  exit 1
fi

#if [ ! -f $1 ] 
#then
#  echo ""
#  echo "Warning: The file $1 does not exist in QF_psl_to_unmap_process.sh, exit."
#  echo ""
#  exit 1
#fi


#Build the name for files
UNMAP_ON_QUERY_PSL_1_SORT=$1"_sort_pl"
UNMAP_ON_QUERY_PSL_2_SORT=$2"_sort_pl"
SINGLE_ON_QUERY_BED_SORT=$3"_sort_pl"
SINGLE_ON_OTHER_BED_SORT=$4"_sort_pl"
UNMAP_ON_QUERY_PSL_1_SORT_TEMP=$1"_temp"
UNMAP_ON_QUERY_PSL_2_SORT_TEMP=$2"_temp"

#for psl, need to remove the headers first
sed -e '1,5d' $1 > $UNMAP_ON_QUERY_PSL_1_SORT_TEMP
perl $5"/sort_perl.pl" $UNMAP_ON_QUERY_PSL_1_SORT_TEMP 10 > $UNMAP_ON_QUERY_PSL_1_SORT
rm $UNMAP_ON_QUERY_PSL_1_SORT_TEMP

sed -e '1,5d' $2 > $UNMAP_ON_QUERY_PSL_2_SORT_TEMP
perl $5"/sort_perl.pl" $UNMAP_ON_QUERY_PSL_2_SORT_TEMP 10 > $UNMAP_ON_QUERY_PSL_2_SORT
rm $UNMAP_ON_QUERY_PSL_2_SORT_TEMP


perl $5"/sort_perl.pl" $3 4 > $SINGLE_ON_QUERY_BED_SORT
perl $5"/sort_perl.pl" $4 4 > $SINGLE_ON_OTHER_BED_SORT