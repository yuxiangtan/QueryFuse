#!/bin/bash

#grep bed by ID
#read ID into R and get the row number matched
#use awk to get the rows

if [ $# -ne 6 ]
then
  echo ""
    echo "Warning: Usage: grep_reads_from_Bed_by_ID.sh ID_list bed_file output_name R_script Script_path LOG_ERR"
  echo ""

  exit 1
fi

if [ ! -f $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist in grep_reads_from_Bed_by_ID.sh, exit."
  echo ""
  exit 1
fi

if [ ! -f $2 ] 
then
  echo ""
  echo "Warning: The file $2 does not exist in grep_reads_from_Bed_by_ID.sh, exit."
  echo ""
  exit 1
fi

if [ ! -x $5 ] 
then
  echo ""
  echo "Warning: Script-path $5 does not exist in grep_reads_from_Bed_by_ID.sh, exit."
  echo ""
  exit 1
fi

TEMP_GEN () {
   RAN_NUM=`date +%s%N`
   TEMP_OUT=${1}${2}${RAN_NUM}
    while [ -f $TEMP_OUT ]; do  
        RAN_NUM=`date +%s%N` 
        TEMP_OUT=${1}${2}${RAN_NUM}  
    done
    echo $TEMP_OUT
}


#get the ID from bed
RAN_NUM=`date +%s%N`
TEMP_BED_ID=`TEMP_GEN $2 "temp_ID"`
cut -f4 $2 | cut -d'/' -f1 > $TEMP_BED_ID
LOG_ERR=$6

#get the row number of matched rows
TEMP_ROW_ID=`TEMP_GEN $3 "temp_ID"`
MATCH_ROW_ID_R=$5"/match_row_ID.R"
if [ -s $TEMP_BED_ID ];
	then
		if [ -s $1 ];
                    then
                        $4 $MATCH_ROW_ID_R file.bed=$TEMP_BED_ID file.list=$1 file.out=$TEMP_ROW_ID
                    else
                        echo "Warning: "$1" is empty in grep_reads_from_Bed_by_ID.sh, exit." >> $LOG_ERR
                        echo "Warning: "$1" is empty in grep_reads_from_Bed_by_ID.sh, exit."
                        rm $TEMP_BED_ID
                        exit 1
                fi
    	else
		echo "Warning: "$2" is empty in grep_reads_from_Bed_by_ID.sh, exit." >> $LOG_ERR
                echo "Warning: "$2" is empty in grep_reads_from_Bed_by_ID.sh, exit."
                exit 1
fi

if [ -s $TEMP_ROW_ID ] 
then
awk 'FNR==NR{a[$1];next}(FNR in a){print}' $TEMP_ROW_ID $2 > $3
fi

rm $TEMP_BED_ID
rm $TEMP_ROW_ID 


