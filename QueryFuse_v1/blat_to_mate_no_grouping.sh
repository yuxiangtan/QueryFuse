#!/bin/bash

#read in the read ID list, for each ID, grep its mate's gene name
#get the location of the gene from whole_gene_list.bed
#fastFromBed to get the seq of gene
#get the location from the substract.bed list
#get the subtract seq
#blat the subtract to the gene


if [ $# -ne 11 ]
then
  echo ""
    echo "Warning: Usage: blat_to_mate_no_grouping.sh read_list mate_bed whole_gene_list.bed reference.fa subtract.bed target_reads.fa output_file Script_path LOG_ERR size_other Align_percent"
    echo ""
    echo "Read_list - A file with one read ID each line."
    echo "mate_bed - A bed file of mate alignment."
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "reference.fa - for example, hg19.fa"
    echo "subtract.bed - A bed file containing the target reads with the location of unaligned subtract"
    echo "target_reads.fa - A fa file containing the sequence of all the target reads"
    echo "outputfile"
    echo "Script_path - folder of QueryFuse"
    echo "LOG_ERR"
    echo "blat step_size to other genes"
    echo "Align_percent"
  echo ""

  exit 1
fi

#Read_list
if [ ! -f $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

#mate_bed
if [ ! -f $2 ] 
then
  echo ""
  echo "Warning: The file $2 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

#whole_gene_list
if [ ! -f $3 ] 
then
  echo ""
  echo "Warning: The file $3 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

#reference.fa
if [ ! -f $4 ] 
then
  echo ""
  echo "Warning: The file $4 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

#subtract.bed
if [ ! -f $5 ] 
then
  echo ""
  echo "Warning: The file $5 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

#target_reads.fa
if [ ! -f $6 ] 
then
  echo ""
  echo "Warning: The file $6 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

#outputfile
if [ -f $7 ] 
then
  rm $7
fi

#Script_path
if [ ! -x $8 ] 
then
  echo ""
  echo "Warning: Script-path $8 does not exist in blat_to_mate_no_grouping.sh, exit."
  echo ""
  exit 1
fi

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


#parameter
RAN_NUM=`date +%s%N`
READ_LIST=$1
MATE_BED=$2
GENE_BED=$3
REF_FA=$4
SUBTRACT_BED=$5
TARGET_FA=$6
outputfile=$7
Script_path=$8
LOG_ERR=${9}
size_other=${10}
Align_percent=${11}
TEMP_BED=`TEMP_GEN $READ_LIST "_"`
TEMP_BED1=`TEMP_GEN $READ_LIST "TEMP1"`
TEMP_BED_NAME=`TEMP_GEN $READ_LIST "TEMP_BED_NAME"`
TEMP_GENE=`TEMP_GEN $READ_LIST "TEMP_GENE"`
#TEMP_ROW_ID=`TEMP_GEN $READ_LIST "TEMP_ROW_ID"`
#TEMP_NAME=`TEMP_GEN $READ_LIST "TEMP_NAME"`
TEMP_GENE_BED=`TEMP_GEN $READ_LIST "TEMP_GENE_BED"`
TEMP_GENE_BED_1=`TEMP_GEN $READ_LIST "TEMP_GENE_BED1"`
TEMP_GENE_FA=`TEMP_GEN $READ_LIST "TEMP_GENE_FA"`
TEMP_SUBTRACT_BED=`TEMP_GEN $READ_LIST "TEMP_SUBTRACT_BED"`
TEMP_SUBTRACT_FA=`TEMP_GEN $READ_LIST "TEMP_SUBTRACT_FA"`
#TEMP_PSL=`TEMP_GEN $READ_LIST "TEMP_psl"`
TEMP_NAME_FA=`TEMP_GEN $READ_LIST "TEMP_NAME_FA"`
TEMP_FA_ROW_DIS=`TEMP_GEN $READ_LIST "TEMP_FA_ROW_DIS"`
rep_match=$((1024*11/size_other))


echo `date`
echo "get mate bed based on ID list in blat_to_mate_no_grouping.sh"
$Script_path"/grep_reads_from_Bed_by_ID.sh" $READ_LIST $MATE_BED $TEMP_BED Rscript $Script_path $LOG_ERR

echo `date`
echo "sort the bed file in blat_to_mate_no_grouping.sh"
if [ -s $TEMP_BED ];
    then
    sort -k17 $TEMP_BED > $TEMP_BED1
    else
        echo $TEMP_BED" is empty in blat_to_mate_no_grouping.sh, exit." >> $LOG_ERR
        echo $TEMP_BED" is empty in blat_to_mate_no_grouping.sh, exit."
        exit 1
fi

echo `date`
echo "get gene names in blat_to_mate_no_grouping.sh"
cut -f17 $TEMP_BED1 | sort -u > $TEMP_GENE
cut -f4 $TEMP_BED1 | cut -d'/' -f1 > $TEMP_BED_NAME


$Script_path"/grep_reads_from_Bed_by_ID.sh" $TEMP_BED_NAME $SUBTRACT_BED $TEMP_SUBTRACT_BED Rscript $Script_path $LOG_ERR
$Script_path"/grep_reads_from_whole_gene_list_Bed_by_ID.sh" $TEMP_GENE $GENE_BED $TEMP_GENE_BED_1 Rscript $Script_path $LOG_ERR
if [ -s $TEMP_GENE_BED_1 ];
    then
    cut -f1,2,3,5 $TEMP_GENE_BED_1 > $TEMP_GENE_BED 
    else
        echo "Warning: "$TEMP_GENE_BED_1" is empty in blat_to_mate_no_grouping.sh, exit." >> $LOG_ERR
        echo "Warning: "$TEMP_GENE_BED_1" is empty in blat_to_mate_no_grouping.sh, exit."
        rm TEMP_BED
        rm TEMP_BED1
        rm TEMP_BED_NAME
        rm TEMP_GENE
        exit 1
fi

echo `date`
echo "get gene FA in blat_to_mate_no_grouping.sh"
fastaFromBed -fi $REF_FA -bed $TEMP_GENE_BED -fo $TEMP_GENE_FA -name


echo `date`
echo "get ROW_DIS in blat_to_mate_no_grouping.sh"
grep "^>" -n -m2 $TARGET_FA | cut -d":" -f1 > $TEMP_FA_ROW_DIS
NUM1=`sed -n "1p" $TEMP_FA_ROW_DIS`
NUM2=`sed -n "2p" $TEMP_FA_ROW_DIS`
ROW_DIS=$[ NUM2 - NUM1 -1 ]
#echo $ROW_DIS

echo `date`
echo "get fa file with target subtract seq in blat_to_mate_no_grouping.sh"

if [ -s $TEMP_BED_NAME ];
	then
        if [ -s $TARGET_FA ];
            then
            perl $Script_path"/extract_fa_by_readID.pl" $TARGET_FA $TEMP_BED_NAME $ROW_DIS $TEMP_NAME_FA
            else
                echo "Warning: "$TARGET_FA" is empty in blat_to_mate_no_grouping.sh, exit." >> $LOG_ERR
                echo "Warning: "$TARGET_FA" is empty in blat_to_mate_no_grouping.sh, exit."
                rm TEMP_BED
                rm TEMP_BED1
                rm TEMP_BED_NAME
                rm TEMP_GENE
                rm TEMP_GENE_BED
                rm TEMP_GENE_FA
                rm TEMP_FA_ROW_DIS
                exit 1
        fi
        else
            echo $TEMP_BED_NAME" is empty in blat_to_mate_no_grouping.sh, exit." >> $LOG_ERR
            echo $TEMP_BED_NAME" is empty in blat_to_mate_no_grouping.sh, exit."
            rm TEMP_BED
            rm TEMP_BED1
            rm TEMP_GENE
            rm TEMP_GENE_BED
            rm TEMP_GENE_FA
            rm TEMP_FA_ROW_DIS
            exit 1
fi

if [ -s $TEMP_NAME_FA".fai" ];
	then
		rm $TEMP_NAME_FA".fai"
fi
	fastaFromBed -fi $TEMP_NAME_FA -bed $TEMP_SUBTRACT_BED -fo $TEMP_SUBTRACT_FA -name

echo `date`
echo "blat to mate fa in blat_to_mate_no_grouping.sh"
echo "blat -stepSize=$size_other -repMatch=$rep_match"
blat -stepSize=$size_other -repMatch=$rep_match -minIdentity=$Align_percent -noHead $TEMP_GENE_FA $TEMP_SUBTRACT_FA $outputfile

rm $TEMP_BED
rm $TEMP_BED1
rm $TEMP_BED_NAME
rm $TEMP_GENE
rm $TEMP_GENE_BED
rm $TEMP_GENE_BED_1
rm $TEMP_GENE_FA
rm $TEMP_SUBTRACT_BED
rm $TEMP_SUBTRACT_FA
rm $TEMP_FA_ROW_DIS
rm $TEMP_NAME_FA
