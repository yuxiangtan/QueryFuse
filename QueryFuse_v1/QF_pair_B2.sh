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


#Gene_specific_fusion_query subfunction: pocessing pair scenario B2
#GSFQ_pair_B2.sh
if [ $# -ne 10 ]
then
  echo ""
    echo "Usage: QF_pair_B2.sh File_prefix BAM_file_folder_prefix result_folder_prefix whole_gene_list.bed human_genome.fa read_length LOG_ERROR QueryFuse_path Align_percent MIN_SCORE"
    echo ""
    echo "File_prefix - The folder that all the defualt can use (for intermedia files)."
    echo "BAM_file_folder_prefix - that has all three needed input bam files."
    echo "result_folder_prefix - folder location for final result."
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "human_genome.fa - fa file of the whole human genome, such as hg19.fa"
    echo "read_length - length of reads"
    echo "LOG_ERROR - path to the error log file"
    echo "QueryFuse_path - QueryFuse path."
    echo "Align_percent"
    echo "MIN_SCORE min_alignment length"
		
    echo ""
    exit 1


fi

file_prefix=$1
bam_fd=$2
outresult_fd=$3
GENE_BED=$4
#in fact, HG_FA is in need
HG_FA=$5
READ_LEN=$6
LOG_ERROR=${7}
QF_path=${8}
Align_percent=${9}
MIN_SCORE=${10}


PAIR_BAM=$bam_fd"paired.bam"
PAIR_SORT=$file_prefix"paired_sorted"
PAIR_SORT_BAM=$bam_fd"paired_sorted.bam"
QUERY_BED=$outresult_fd"query_gene.bed"
QUERY_FA=$outresult_fd"query_gene.fa"
PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"
PAIR_TO_QUERY_FILTER_ID_DUPLI=$PAIR_SORT"_to_query_pair_filter_ID_duplicate.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_FA=$PAIR_SORT"_to_query_pair_filter_ID_duplicate.fa"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID.psl"


#pair b2
#get the ID with only one output from pair b
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL ];
        then
            echo "pre-process scenario B2 in QF_pair_B2.sh"
            echo `date`
            sed '1,5d' $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL | cut -f10 |uniq -u > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID
            $QF_path"grep_reads_from_psl_withheader_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL Rscript $QF_path
            #grep  -f $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL
        else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL" is not exist in QF_pair_B2.sh, exit"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL" is not exist in QF_pair_B2.sh, exit" >> $LOG_ERROR
            exit 1
fi

PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_MATCH=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_match.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_ROW=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_row.txt"
#PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL_UP=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID_temp_up.psl" 
#PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.txt"
#PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_FILTER=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.txtfilter.txt"
#PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_BAM=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable.bam"
#PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPAN_ID_UNRELIABLE_ANNO=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_span_ID_unreliable_gene_anno.bed"
#If match >=76, unreliable spanning read candidates. Get the information of mates:
#get IDs 
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL ];
	then
	    echo "get unreliable spanning read candidates in scenario B2 in QF_pair_B2.sh"
            sort -g $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL
	else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL" is not exist in QF_pair_B2.sh, exit"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL" is not exist in QF_pair_B2.sh, exit" >> $LOG_ERROR
            exit 1
fi

if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL ];
	then
            cut -f1 $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_MATCH
            Rscript $QF_path"boundary_row_ID.R" file.in=$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_MATCH lower_bound=0 upper_bound=$[READ_LEN-MIN_SCORE] file.out=$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_ROW
            UP_ROW_DUP_SPLIT=`sed '2p' -n $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_ROW`
            ROW_NUM_DUP_SPLIT=`wc -l $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_MATCH | cut -d' ' -f1`
            echo "the upper boundary is "$UP_ROW_DUP_SPLIT" and should be matched by the next number"
	else
            UP_ROW_DUP_SPLIT="NA"
            ROW_NUM_DUP_SPLIT="NA"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL" is not exist in QF_pair_B2.sh, exit"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_TEMP_PSL" is not exist in QF_pair_B2.sh, exit" >> $LOG_ERROR
            exit 1
fi

echo "finished pre-processing of pair scenario b2"
echo "===================="
$QF_path"QF_pair_B2.2.sh" $PAIR_SORT $UP_ROW_DUP_SPLIT $ROW_NUM_DUP_SPLIT $GENE_BED $QUERY_BED $READ_LEN $QF_path $LOG_ERROR $Align_percent

echo "finished processing of pair scenario b2.2 in QF_pair_B2.sh"
echo "===================="

$QF_path"QF_pair_B2.1.sh" $file_prefix $UP_ROW_DUP_SPLIT $ROW_NUM_DUP_SPLIT $GENE_BED $QF_path $LOG_ERROR
echo "finished processing of pair scenario b2.1 in QF_pair_B2.sh"
echo "===================="