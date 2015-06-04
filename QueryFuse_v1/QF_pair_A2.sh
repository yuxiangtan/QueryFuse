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

#Gene_specific_fusion_query subfunction: pocessing pair scenario A2
#GSFQ_pair_A2.sh
if [ $# -ne 8 ]
then
  echo ""
    echo "Usage: QF_pair_A2.sh file_prefix READ_LEN QUERY_FA whole_gene_list.bed tophat_genome_reference LOG_ERROR QueryFuse_path size_query"
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
PAIR_TO_QUERY_FILTER_ID_UNIQ=$PAIR_SORT"_to_query_filter_ID_uniq.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_FA=$PAIR_SORT"_to_query_filter_ID_uniq.fa"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_ID=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_split_ID.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_split_ID.psl"

#Separate splitting candidates.(assume max duplicate level is 2)
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL ];
        then
            echo `date`
            echo 'Separate splitting candidates in QF_pair_A2.sh'
            sed '1,5d' $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL | cut -f10 |uniq -d > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_ID
            $QF_path"grep_reads_from_psl_withheader_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_ID $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_PSL Rscript $QF_path
            #grep  -f $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_ID $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_PSL
        else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL " is not exist in QF_pair_A2.sh, exit"             
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL " is not exist in QF_pair_A2.sh, exit" >> $LOG_ERROR
            exit
fi


echo "finished pre-processing of pair scenario a2 in QF_pair_A2.sh"
echo "===================="
$QF_path"QF_pair_unreliable_split.sh" $file_prefix $READ_LEN $QUERY_FA $GENE_BED $GENOME_REF $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_SPLIT_PSL $PAIR_TO_QUERY_FILTER_ID_UNIQ_FA $QF_path $LOG_ERROR $size_query
echo "finished processing pair scenario a2 in QF_pair_A2.sh"
echo "===================="
