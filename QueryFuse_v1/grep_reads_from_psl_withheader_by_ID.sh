#!/bin/bash

#grep psl by ID
#this function will sed the header
#read ID into R and get the row number matched
#use awk to get the rows

if [ $# -ne 5 ]
then
  echo ""
    echo "Usage: grep_reads_from_psl_withheader_by_ID.sh ID_list psl_file output_name R_script QueryFuse_path"
    echo "R_script - path of a specific R version."
    echo "QueryFuse_path - QueryFuse path."
    
  echo ""

  exit 1
fi

if [ ! -s $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist in grep_reads_from_psl_withheader_by_ID.sh, exit"
  echo ""
  exit 1
fi

if [ ! -s $2 ] 
then
  echo ""
  echo "Warning: The file $2 does not exist in grep_reads_from_psl_withheader_by_ID.sh, exit"
  echo ""
  exit 1
fi

TEMP_psl=$2"temp_psl"
sed '1,5d' $2 > $TEMP_psl

#get the ID from bed
TEMP_BED_ID=$2"temp_ID.txt"
cut -f10 $TEMP_psl | cut -d'/' -f1 > $TEMP_BED_ID

#get the row number of matched rows
TEMP_ROW_ID=$3"temp_ID.txt"
$4 $5"match_row_ID.R" file.bed=$TEMP_BED_ID file.list=$1 file.out=$TEMP_ROW_ID

if [ -s $TEMP_ROW_ID ] 
    then
        awk 'FNR==NR{a[$1];next}(FNR in a){print}' $TEMP_ROW_ID $TEMP_psl > $3
    else
        echo "Warning: match_row_ID.R failed to produce result in grep_reads_from_psl_withheader_by_ID.sh."
fi

rm $TEMP_BED_ID
rm $TEMP_ROW_ID 


