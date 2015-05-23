#!/bin/bash

#Gene_specific_fusion_query subfunction: pocessing pair scenario B2.2
#GSFQ_pair_B2.2.sh
if [ $# -ne 9 ]
then
  echo ""
    echo "Usage: QF_pair_B2.2.sh PAIR_SORT UP_ROW ROW_NUM whole_gene_list.bed QUERY_BED read_length QF_path LOG_ERROR Align_percent"
    echo ""
    echo "PAIR_SORT - PATH of PAIR_SORT"
    echo "UP_ROW - row number of upper bound"
    echo "ROW_NUM - The total number of rows."
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "QUERY_BED"
    echo "read_length - length of reads"
    echo "QF_path"
    echo "LOG_ERROR"
    echo "Align_percent"
    echo ""
    exit 1
fi
PAIR_SORT=$1
UP_ROW_DUP_SPLIT=$2
ROW_NUM_DUP_SPLIT=$3
GENE_BED=$4
QUERY_BED=$5
READ_LEN=$6
QF_path=${7}
LOG_ERROR=${8}
Align_percent=${9}


PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_SAM=$PAIR_SORT"_to_query.sam"
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_MATCH=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_match.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_ROW=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_row.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_up.psl" 
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP_GOOD=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_up_good.psl" 
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_NO_HEADER_SAM=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable_no_header.sam"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_FILTER=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.txtfilter.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_SAM=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.sam"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.bam"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable_gene_anno.bed"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO_PAIR=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable_gene_anno_pair.bed"


#pair B2.2
if [ -n "$UP_ROW_DUP_SPLIT" ];
	then
		if [ "$UP_ROW_DUP_SPLIT" != "NA" ];
			then
			#echo `date`
			echo "UP_ROW_DUP_SPLIT: "$UP_ROW_DUP_SPLIT
			echo 'get unreliable span candidate in QF_pair_B2.2.sh'
			sed "${UP_ROW_DUP_SPLIT},${ROW_NUM_DUP_SPLIT}p" -n $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP
			#to reduce false positive support reads, control the maximum mismatch to at most 1/100
			python $QF_path"eliminate_mismatch_read_psl.py" $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP $READ_LEN $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP_GOOD $Align_percent
			cut -f10 $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP_GOOD > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE
                #No else needed because if this number equal to NA, then it should exited in QF_pair_B2.sh 
                fi
	else
		echo "Warning: Upper bound is not available in QF_pair_B2.2.sh, exit."
		echo "Warning: Upper bound is not available in QF_pair_B2.2.sh, exit." >> $LOG_ERROR
                exit 1
fi

#get filter script, not used anymore, bamtools is too slow.
#if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE ];
#	then
#		perl read_ID_to_bamfilter.pl ./ $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE ./
#		echo $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE" is not exist"
#	else
#		echo $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE" is not exist"
#fi
#extract alignment information
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE ];
    then
            #echo `date`
            #echo "get bam for unreliabe span reads (b2.2) in QF_pair_B2.2.sh"
            #bamtools filter -in $PAIR_TO_QUERY_BAM -out $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM -script $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_FILTER
            echo "grep_reads from sam in QF_pair_B2.2.sh"
            $QF_path"grep_reads_from_SAM_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE $PAIR_TO_QUERY_SAM $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_NO_HEADER_SAM Rscript $QF_path
            #echo `date`
            echo "add header in QF_pair_B2.2.sh"
            samtools view -H $PAIR_TO_QUERY_BAM | cat > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_SAM
            cat $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_NO_HEADER_SAM >> $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_SAM
            echo "sam to bam in QF_pair_B2.2.sh"
            samtools view -S $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_SAM -b > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM
            #echo `date`
            rm $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_NO_HEADER_SAM
    else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE" is not exist in QF_pair_B2.2.sh, exit."
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE" is not exist in QF_pair_B2.2.sh, exit." >> $LOG_ERROR
            exit 1
fi


#get gene information by intersectBed
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM ];
        then
                #echo `date`
                echo "get gene information for unreliabe span reads in QF_pair_B2.2.sh"
                intersectBed -abam $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM -b $4 -bed -wo > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM" is not exist in QF_pair_B2.2.sh, exit."
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM" is not exist in QF_pair_B2.2.sh, exit."  >> $LOG_ERROR
                exit 1
fi

#summary of the result by reading bed file.(by compare pair, one on query, one is ont.
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO_PAIR_TEMP=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable_gene_anno_pair.bed_temp"

if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO ];
        then
                perl $QF_path"QF_pair_A1.1_group_to_pair.pl" $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO_PAIR_TEMP $QUERY_BED $READ_LEN $Align_percent
                sort -u $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO_PAIR_TEMP > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO_PAIR
                rm $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO_PAIR_TEMP
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO" is not exist in QF_pair_B2.2.sh"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO" is not exist in QF_pair_B2.2.sh"  >> $LOG_ERROR
                
fi

