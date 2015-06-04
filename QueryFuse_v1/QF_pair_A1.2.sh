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

#Gene_specific_fusion_query subfunction: pocessing pair scenario A1.2
#GSFQ_pair_A1.2.sh
if [ $# -ne 6 ]
then
  echo ""
    echo "Usage: QF_pair_A1.2.sh PAIR_SORT UP_ROW ROW_NUM whole_gene_list.bed QF_path LOG_ERROR"
    echo ""
    echo "PAIR_SORT - PATH of PAIR_SORT"
    echo "UP_ROW - row number of upper bound"
    echo "ROW_NUM - The total number of rows."
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "QF_path"
    echo "LOG_ERROR"
    echo ""
    exit 1
fi


#name the parameters
PAIR_SORT=$1
UP_ROW=$2
ROW_NUM=$3
GENE_BED=$4
QF_path=${5}
LOG_ERROR=${6}


PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_SAM=$PAIR_SORT"_to_query.sam"
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"
PAIR_TO_QUERY_FILTER_ID_UNIQ=$PAIR_SORT"_to_query_filter_ID_uniq.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_FA=$PAIR_SORT"_to_query_filter_ID_uniq.fa"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_ID.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_MATCH=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_match"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_ROW=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_row"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_LOW=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_low.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_UNRELIABLE_ID.txt"
#PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_FILTER=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_UNRELIABLE_ID.txtfilter.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_NO_HEADER_SAM=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_UNRELIABLE_ID_extract_no_head.sam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_SAM=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_UNRELIABLE_ID_extract.sam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_UNRELIABLE_ID_extract.bam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_ANNO_BED=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_UNRELIABLE_gene_anno.bed"


#If match <read_length-23, they are unreliable split/span reads
if [ -n "$UP_ROW" ];
	then
		if [ "$UP_ROW" -gt 1 ];
			then
			#echo `date`
			echo 'get less reliable split/span reads'
			UP_ROW_1=$[ UP_ROW -1 ]
			sed "1,${UP_ROW_1}p" -n $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_LOW
			cut -f10 $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_LOW > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE
		fi
	else
		echo "Warning: upper bound number is not available in QF_pair_A1.2.sh, exit."
		echo "Warning: upper bound number is not available in QF_pair_A1.2.sh, exit." >> $LOG_ERROR
                exit 1
fi


#if there are reads fulfil the requirement
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE ];
	then
            #extract alignment information from bam file using ID
            #use ID list to generate bamtools filter script
            #perl read_ID_to_bamfilter.pl ./ $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE ./
            #echo ""
            #
            #bamtools filter -in $PAIR_TO_QUERY_BAM -out $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM -script $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_FILTER
            #echo `date`
            #echo "bam to sam in QF_pair_A1.2.sh"
            #samtools view $PAIR_TO_QUERY_BAM > $PAIR_TO_QUERY_SAM
            #echo `date`
            echo "grep_reads from sam in QF_pair_A1.2.sh"
            $QF_path"grep_reads_from_SAM_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE $PAIR_TO_QUERY_SAM $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_NO_HEADER_SAM Rscript $QF_path
            #echo `date`
            #echo "add header in QF_pair_A1.2.sh"
            samtools view -H $PAIR_TO_QUERY_BAM | cat > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_SAM
            cat $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_NO_HEADER_SAM >> $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_SAM
            #echo "sam to bam in QF_pair_A1.2.sh"
            samtools view -S $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_SAM -b > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM
            #echo `date`
            rm $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_NO_HEADER_SAM
        else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE" is not generated in QF_pair_A1.2.sh"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE" is not generated in QF_pair_A1.2.sh" >> $LOG_ERROR
            exit 1
fi

if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM ];
	then
            echo "get gene information by intersectBed in QF_pair_A1.2.sh"
            #get gene information by intersectBed
             intersectBed -abam $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM -b $GENE_BED -bed -wo > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_ANNO_BED
            #summary of the result by reading bed file.
            #
            #
	else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM" is not generated in QF_pair_A1.2.sh"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_UNRELIABLE_EXTRACT_BAM" is not generated in QF_pair_A1.2.sh" >> $LOG_ERROR
fi
