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


#grep SAM by ID
#read ID into R and get the row number matched
#use awk to get the rows

if [ $# -ne 5 ]
then
  echo ""
    echo "Usage: grep_reads_from_SAM_by_ID.sh ID_list SAM_file output_name R_script QF_path"
  echo ""

  exit 1
fi

if [ ! -s $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist in grep_reads_from_SAM_by_ID.sh, exit"
  echo ""
  exit 1
fi

if [ ! -s $2 ] 
then
  echo ""
  echo "Warning: The file $2 does not exist in grep_reads_from_SAM_by_ID.sh, exit"
  echo ""
  exit 1
fi

#get the ID from SAM
TEMP_SAM_ID=$2"temp_ID.txt"
cut -f1 $2 > $TEMP_SAM_ID

#get the row number of matched rows
TEMP_ROW_ID=$3"temp_ID.txt"
$4 $5"match_row_ID.R" file.bed=$TEMP_SAM_ID file.list=$1 file.out=$TEMP_ROW_ID

if [ -s $TEMP_ROW_ID ] 
    then
        awk 'FNR==NR{a[$1];next}(FNR in a){print}' $TEMP_ROW_ID $2 > $3
    else
        echo "Warning: match_row_ID.R failed to produce result in grep_reads_from_SAM_by_ID.sh"
fi

rm $TEMP_SAM_ID
rm $TEMP_ROW_ID 


