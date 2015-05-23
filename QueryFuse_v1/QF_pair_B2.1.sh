#!/bin/bash

#Gene_specific_fusion_query subfunction: pocessing pair scenario B2.1
#GSFQ_pair_B2.1.sh
if [ $# -ne 6 ]
then
  echo ""
    echo "Usage: QF_pair_B2.1.sh PAIR_SORT UP_ROW ROW_NUM whole_gene_list.bed"
    echo ""
    echo "file_prefix - PATH of the intermediate output folder"
    echo "UP_ROW - row number of upper bound"
    echo "ROW_NUM - The total number of rows."
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "QF_path"
    echo "LOG_ERROR"
    echo ""
    exit 1
fi
file_prefix=$1
UP_ROW_DUP_SPLIT=$2
ROW_NUM_DUP_SPLIT=$3
GENE_BED=$4
QF_path=${5}
LOG_ERROR=${6}

PAIR_SORT=$file_prefix"paired_sorted"
PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"

PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_low.psl" 
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_good.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_NAME="paired_sorted_to_query_pair_filter_ID_duplicate_on_query_split_ID_good.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_FILTER=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_good.txtfilter.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_BAM=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_good.bam"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_ANNO=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_good_gene_anno.bed"

#If match <76, splitting read candidates, get the information of mates, and realign the subtract to the mate gene region (b2.1)
#get IDs 
if [ -n "$UP_ROW_DUP_SPLIT" ];
	then
		if [ "$UP_ROW_DUP_SPLIT" != "NA" ];
                then
                    if [ "$UP_ROW_DUP_SPLIT" -gt 1 ];
		    then
			#echo `date`
			echo 'extract information for b2.1 splitting candidates in QF_pair_B2.1.sh'
			UP_ROW_DUP_SPLIT_1=$[ UP_ROW_DUP_SPLIT -1 ]
			sed "1,${UP_ROW_DUP_SPLIT_1}p" -n $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW
			cut -f10 $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD
                    else
                        echo "Warning: upper bound number is 0 in QF_pair_B2.1.sh, exit"
                        echo "Warning: upper bound number is 0 in QF_pair_B2.1.sh, exit" >> $LOG_ERROR
                        exit 1
                    fi
                fi
	else
		echo "Warning: upper bound number is not available in QF_pair_B2.1.sh, exit"
		echo "Warning: upper bound number is not available in QF_pair_B2.1.sh, exit" >> $LOG_ERROR
                exit 1
fi
#get filter script
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD ];
	then
		echo "transform read ID into bam filter script in QF_pair_B2.1.sh"
		perl $QF_path"read_ID_to_bamfilter.pl" $file_prefix $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_NAME $file_prefix
		#echo ""
	else
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD" is not exist in QF_pair_B2.1.sh, exit"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD" is not exist in QF_pair_B2.1.sh, exit" >> $LOG_ERROR
                exit 1

fi

#extract alignment information to bed format
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_FILTER ];
        then
                echo `date`
                echo "get bam for good split reads (b2.1) in QF_pair_B2.1.sh"
                bamtools filter -in $PAIR_TO_QUERY_BAM -out $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_BAM -script $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_FILTER
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_FILTER" is not exist in QF_pair_B2.1.sh, exit"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_FILTER" is not exist in QF_pair_B2.1.sh, exit" >> $LOG_ERROR
                exit 1
fi

#get gene information by intersectBed
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_BAM ];
        then
                echo `date`
                echo "get gene information for good split reads (b2.1)"
                intersectBed -abam $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_BAM -b $4 -bed -wo > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_ANNO
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_BAM" is not exist in QF_pair_B2.1.sh"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_GOOD_BAM" is not exist in QF_pair_B2.1.sh" >> $LOG_ERROR
fi
#For each ID, get the end info, filter it out, transform to fa format (not doing this yet, optional)
#get the seq of the sub-fraction, realign to query gene. (not doing this yet, optional)
#summary of the result by reading bed file.
#
#
#
#
