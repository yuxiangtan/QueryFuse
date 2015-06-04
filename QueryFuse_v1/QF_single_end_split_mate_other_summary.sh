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

#only process one end at a time. To process both, need to run for each end.
#QF_single_end_split_mate_other_summary.sh
#summary the split_mate_to_other_blat_psl
#for each read, get what genes are the top three.
#find the partner of those reads, they should align to the same 

if [ $# -ne 7 ]
then
  echo ""
    echo "Warning: Usage: QF_single_end_split_mate_other_summary.sh SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL SINGLE_ON_OTHER_BED UNMAP_ON_QUERY_PSL read_length Script_path LOG_ERR Align_percent"
    echo ""
    echo "SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL - location of the SINGLE_ON_QUERY_SPLIT_MATE_OTHER_BLAT_PSL."
    echo "SINGLE_ON_OTHER_BED - bed file about genes where singletons aligned to, except query."
    echo "UNMAP_ON_QUERY_PSL - psl file of UNMAP_1/2_ON_QUERY_PSL"
    echo "read_length - length of reads"
    echo "Script_path - folder of QueryFuse"
    echo "LOG_ERR"
    echo "Align_percent"
    echo ""

  exit 1
fi

if [ ! -s $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist, exit."
  echo ""
  exit 1
fi

if [ ! -s $2 ] 
then
  echo ""
  echo "Warning: The file $2 does not exist, exit."
  echo ""
  exit 1
fi

if [ ! -s $3 ] 
then
  echo ""
  echo "Warning: The file $3 does not exist, exit."
  echo ""
  exit 1
fi

if [ ! -x $5 ] 
then
  echo ""
  echo "Warning: Script-path $5 does not exist, exit."
  echo ""
  exit 1
fi
read_length=$4
Script_path=$5
LOG_ERR=$6
Align_percent=$7
#TEMP_GENERATION function to generate temp file and also check existence. If exist, regenerate. Must have two arguments
TEMP_GEN () {
   RAN_NUM=`date +%s%N`
   TEMP_OUT=${1}${2}${RAN_NUM}
    while [ -f $TEMP_OUT ]; do  
        RAN_NUM=`date +%s%N` 
        TEMP_OUT=${1}${2}${RAN_NUM}  
    done
    echo $TEMP_OUT
}

#get the top3 for each reads.
perl $Script_path"/split_mate_to_other_blat_psl_summary.pl" $1

SUMMARY_BED=$1"_summary.bed"
SORTED_SUMMARY_BED=$1"_summary_sorted.bed"
SINGLE_ON_OTHER_BED_SORT=$2

#echo "sort the output by read_ID in QF_single_end_split_mate_other_summary.sh"
#UNMAP_ON_QUERY_PSL_TEMP=`TEMP_GEN $3 "TEMP"`
UNMAP_ON_QUERY_PSL_SORT=$3"_sort_pl"
perl $Script_path"/sort_perl.pl" $SUMMARY_BED 10 > $SORTED_SUMMARY_BED
#sed -e '1,5d' $3 > $UNMAP_ON_QUERY_PSL_TEMP
#perl $Script_path"/sort_perl.pl" $UNMAP_ON_QUERY_PSL_TEMP 10 > $UNMAP_ON_QUERY_PSL_SORT
#rm $UNMAP_ON_QUERY_PSL_TEMP
#perl sort_perl.pl $SINGLE_ON_OTHER_BED 4 > $SINGLE_ON_OTHER_BED_SORT

#echo "find its partner and match them in QF_single_end_split_mate_other_summary.sh"
#function is to figure out which has highest match,
perl $Script_path"/split_mate_to_other_grouped_into_triple_summary.pl" $SORTED_SUMMARY_BED $SINGLE_ON_OTHER_BED_SORT $UNMAP_ON_QUERY_PSL_SORT
SORTED_SUMMARY_BED_MERGE=$1"_summary_sorted.bed_merge_with_mate_summary.bed"
#group into pairs and also filter our pair with too long or too much mismatch
perl $Script_path"/split_mate_to_other_triple_summary_transform.pl" $SORTED_SUMMARY_BED_MERGE $read_length $Align_percent