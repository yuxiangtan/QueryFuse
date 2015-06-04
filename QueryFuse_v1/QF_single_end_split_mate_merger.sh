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

#Gene_specific_fusion_query subfunction: single_end pocessing (7) span mate to other only and (8) span mate to other for to both (but output name overlapped)
#only process one end at a time. To process both, need to run for each end.
#GSFQ_single_end_split_mate_merger.sh (for mate_query or mate_both)
#merge the triple subset of a pair together
#the sub on other is in bed format, need to transform it into qsl like format,
#find the partner of those reads, they should align to the same 

if [ $# -ne 6 ]
then
  echo ""
    echo "Warning: Usage: QF_single_end_split_mate_merger.sh annotation_bed_file SINGLE_ON_QUERY_BED UNMAP_ON_QUERY_PSL Read_length QueryFuse_path align_percent_min"
    echo ""
    echo "annotation_bed_file - location of the SINGLE_ON_QUERY_SPLIT_MATE_QUERY_ID_SUBTRACT_ANNO."
    echo "SINGLE_ON_QUERY_BED - bed file about genes where singletons aligned to, except query."
    echo "UNMAP_ON_QUERY_PSL - psl file of UNMAP_1/2_ON_QUERY_PSL"
    echo "read_length - length of reads."
    echo "QueryFuse_path - QueryFuse path."
    echo "align_percent_min - min percentage of alignment."
    
    echo ""

  exit 1
fi

if [ ! -s $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist in QF_single_end_split_mate_merger.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $2 ] 
then
  echo ""
  echo "Warning: The file $2 does not exist in QF_single_end_split_mate_merger.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $3 ] 
then
  echo ""
  echo "Warning: The file $3 does not exist in QF_single_end_split_mate_merger.sh, exit."
  echo ""
  exit 1
fi

if [ ! -x $5 ] 
then
  echo ""
  echo "Warning: Script-path $5 does not exist in QF_single_end_split_mate_merger.sh, exit."
  echo ""
  exit 1
fi

SORTED_ANNO_BED=$1"_sort"
read_length=$4
QueryFuse_path=$5
align_percent_min=$6
#echo "sort the output by read_ID"
UNMAP_ON_QUERY_PSL_SORT=$3"_sort_pl"
perl $QueryFuse_path"/sort_perl.pl" $1 4 > $SORTED_ANNO_BED
SORTED_SINGLE=$2"_sort_pl"
#perl $QueryFuse_path"/sort_perl.pl" $2 4 > $SORTED_SINGLE

#echo "find its partner and match them"
perl $QueryFuse_path"/split_mate_to_query_grouped_into_triple_summary.pl" $SORTED_ANNO_BED $SORTED_SINGLE $UNMAP_ON_QUERY_PSL_SORT $read_length $align_percent_min
