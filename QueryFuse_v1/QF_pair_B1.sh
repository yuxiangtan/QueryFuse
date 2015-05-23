#!/bin/bash

#Gene_specific_fusion_query subfunction: pocessing pair scenario B1
#GSFQ_pair_B1.sh
if [ $# -ne 8 ]
then
  echo ""
    echo "Usage: QF_pair_B1.sh file_prefix READ_LEN QUERY_FA whole_gene_list.bed tophat_genome_reference LOG_ERROR QueryFuse_path size_query"
    echo ""
    echo "file_prefix -  The folder that all the defualt can use (for intermedia files)."
    echo "READ_LEN - the length of reads"
    echo "QUERY_FA - PATH of QUERY_FA"
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "tophat_genome_reference - indexed genome_reference."
    echo "LOG_ERROR - path to the error log file"
    echo "QueryFuse_path - QueryFuse path."
    echo "size_query - step_size for blat to query"
    echo ""
    exit 1
fi


#name the parameters
file_prefix=$1
GENE_BED=$4
READ_LEN=$2
GENOME_REF=$5
QUERY_FA=$3
LOG_ERROR=$6
QF_path=$7
size_query=${8}

PAIR_SORT=$file_prefix"paired_sorted"
PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"
PAIR_TO_QUERY_FILTER_ID_DUPLI_FA=$PAIR_SORT"_to_query_pair_filter_ID_duplicate.fa"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_TEMP_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_temp.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_TEMP_MATCH=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_temp_match.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_TEMP_ROW=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_temp_row.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_TEMP_PSL_LOW=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_temp_low.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_SPLIT=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_split.txt"

#get the ID with two output
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL ];
        then
        echo `date`
        echo "preprocessing Pair B1 scenario in QF_pair_B1.sh"
        sed '1,5d' $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL | cut -f10 |uniq -d > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID
        $QF_path"grep_reads_from_psl_withheader_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_PSL Rscript $QF_path
        #grep  -f $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_PSL
        else
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL" is not generated in QF_pair_B1.sh, exit"
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL" is not generated in QF_pair_B1.sh, exit" >> $LOG_ERROR
                exit 1
fi


echo "finished pre-processing of pair scenario b1 in QF_pair_B1.sh"
echo "===================="
#echo $PAIR_TO_QUERY_FILTER_ID_DUPLI_FA
$QF_path"QF_pair_unreliable_split.sh" $file_prefix $READ_LEN $QUERY_FA $GENE_BED $GENOME_REF $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_PSL $PAIR_TO_QUERY_FILTER_ID_DUPLI_FA $QF_path $LOG_ERROR $size_query

echo "finished processing pair scenario b1 in QF_pair_B1.sh"
echo "===================="
