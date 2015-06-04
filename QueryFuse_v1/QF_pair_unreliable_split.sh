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

#Gene_specific_fusion_query subfunction: pocessing unreliable split reads (A2 and B1)
#GSFQ_pair_unreliable_split.sh
if [ $# -ne 10 ];
then
  echo ""
    echo "Usage: QF_pair_unreliable_split.sh file_prefix READ_LEN QUERY_FA whole_gene_list.bed tophat_genome_reference QUERY_SPLIT_PSL QUERY_SPLIT_CANDIDATE_ID QUERY_FILTER_ID_FA QF_path LOG_ERROR size_query"
    echo ""
    echo "file_prefix -  The folder that all the defualt can use (for intermedia files)."
    echo "READ_LEN - the length of reads"
    echo "QUERY_FA - PATH of QUERY_FA"
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "tophat_genome_reference - indexed genome_reference."
    echo "QUERY_SPLIT_PSL - PATH of QUERY_SPLIT_PSL"
    echo "QUERY_FILTER_ID_FA - PATH of QUERY_FILTER_ID_FA"
    echo "QF_path"
    echo "LOG_ERROR"
    echo "size_query - step_size for blat to query"
    echo ""
	exit 1
fi

file_prefix=$1
READ_LEN=$2
QUERY_FA=$3
whole_gene_list=$4
genome_ref=$5
QUERY_SPLIT_PSL=$6
QUERY_FILTER_ID_FA=$7
QF_path=${8}
LOG_ERROR=${9}
size_query=${10}
rep_match=$((1024*11/size_query))

PAIR_SORT=$file_prefix"paired_sorted"
PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"
QUERY_SPLIT_PSL_NAME=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_split_ID.psl"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA=$QUERY_SPLIT_PSL".fa"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_PSL=$QUERY_SPLIT_PSL".psl"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_CANDIDATE_ID=$QUERY_SPLIT_PSL".candiate_ID.txt"



#For each ID, grep the pairs from psl, If match <76 splitting read candidates. Max match>=90, record 1 or 2 end, then the other end must <76. If not, filter out.
if [ -s $QUERY_SPLIT_PSL ];
        then
                #echo `date`
                echo "get split candiates and their seq into fa file in QF_pair_unreliable_split.sh"
                #in this small script, grep is used, it may make the program slow. 
                $QF_path"psl_to_ID_extract_split_major_on_query.sh" $QUERY_SPLIT_PSL $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_CANDIDATE_ID $QUERY_FILTER_ID_FA $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA $READ_LEN
        else
                echo "Warning: "$QUERY_SPLIT_PSL" is not exist in QF_pair_unreliable_split.sh, exit"
                echo "Warning: "$QUERY_SPLIT_PSL" is not exist in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
                exit 1
fi

#Blat again to make sure which end is split
if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA ];
        then
                echo `date`
                echo 'Blat to get splitting end from candidates in QF_pair_unreliable_split.sh'
                blat -stepSize=$size_query -repMatch=$rep_match -noHead $QUERY_FA $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_PSL
                echo `date`
                echo 'Blat finished in QF_pair_unreliable_split.sh'
                
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA" is not generated in QF_pair_unreliable_split.sh, exit"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA" is not generated in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
                exit 1
fi

PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL=$QUERY_SPLIT_PSL"_sorted_temp.psl"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_MATCH=$QUERY_SPLIT_PSL"_sorted_temp_match.txt"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_ROW=$QUERY_SPLIT_PSL"_sorted_temp_row.txt"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW=$QUERY_SPLIT_PSL"_sorted_temp_row_low.psl"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT=$QUERY_SPLIT_PSL"_sorted_split_ID.txt"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT_NAME=$QUERY_SPLIT_PSL_NAME"_sorted_split_ID.txt"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT_FILTER=$QUERY_SPLIT_PSL"_sorted_split_ID.txtfilter.txt"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_BED=$QUERY_SPLIT_PSL"_sorted_temp_row_low_subtract.bed"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA_FAI=$QUERY_SPLIT_PSL".fa.fai"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_FA=$QUERY_SPLIT_PSL"_sorted_temp_row_low_subtract.fa"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT=$QUERY_SPLIT_PSL"_sorted_temp_row_low_subtract_tophat_out"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT_BAM=$QUERY_SPLIT_PSL"_sorted_temp_row_low_subtract_tophat_out/accepted_hits.bam"
PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_ANNO=$QUERY_SPLIT_PSL"_sorted_temp_row_low_subtract_anno.bed"

if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_PSL ]
	then
            sort -g $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_PSL > $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL

	else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_PSL" is not generated in QF_pair_unreliable_split.sh, exit"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_PSL" is not generated in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
            exit 1
fi

if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL ]
	then
            cut -f1 $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_MATCH
            Rscript $QF_path"boundary_row_ID.R" file.in=$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_MATCH lower_bound=0 upper_bound=$[READ_LEN-23] file.out=$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_ROW
            UP_ROW_SPLIT=`sed '2p' -n $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_ROW`
            ROW_NUM_SPLIT=`wc -l $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_MATCH | cut -d' ' -f1`
	
	else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL" is not generated in QF_pair_unreliable_split.sh, exit"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL" is not generated in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
            exit 1
fi

if [ -n "$UP_ROW_SPLIT" ];
	then
		if [ "$UP_ROW_SPLIT" -gt 1 ];
			then
			#echo `date`
			#echo 'extract subtract information for splitting candidates'
			UP_ROW_SPLIT_1=$[ UP_ROW_SPLIT -1 ]
			sed "1,${UP_ROW_SPLIT_1}p" -n $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW
			cut -f10 $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW > $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT
		fi
	else
		echo "Warning: Upper bound row number is not available in QF_pair_unreliable_split.sh, exit"
		echo "Warning: Upper bound row number is not available in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
                exit 1
fi

#if there are reads fulfil the requirement
if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT ];
	then
            #echo ""
            #use ID list and psl to get subtract.bed
            echo "get subtract_bed in QF_pair_unreliable_split.sh"
            Rscript $QF_path"subtract_psl.R" file.in=$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW file.out=$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_BED read_length=$READ_LEN

	else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT" is not generated in QF_pair_unreliable_split.sh, exit"
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_SPLIT" is not generated in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
            exit 1
fi 


#get the subtract fa file 
if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA_FAI ];
	then
		rm $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA_FAI
fi		

if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA -a -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_BED ];
then
	fastaFromBed -fi $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA -bed $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_BED -fo $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_FA -name
else
	echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA" and/or "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_BED" is not exist in QF_pair_unreliable_split.sh, exit"
        echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_FA" and/or "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_BED" is not exist in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
        exit 1
fi

#get gene information by Tophat: 
if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_FA ];
	then
		#echo `date`
		echo 'Get splitting end location by tophat in QF_pair_unreliable_split.sh'
		tophat \
                --segment-length 12 \
		-o $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT \
		$genome_ref \
		$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_FA
	else
	echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_FA" is not generated in QF_pair_unreliable_split.sh, exit"
	echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_FA" is not generated in QF_pair_unreliable_split.sh, exit" >> $LOG_ERROR
        exit 1
fi

#get the alignment location from Tophat output:
if [ -s $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT_BAM ];
	then
		intersectBed -abam $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT_BAM -b $whole_gene_list -bed -wo > $PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_ANNO
	else
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT_BAM" is not generated in QF_pair_unreliable_split.sh"
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_ON_QUERY_SPLIT_ID_TEMP_PSL_LOW_SUBTRACT_TOPHAT_OUT_BAM" is not generated in QF_pair_unreliable_split.sh" >> $LOG_ERROR
fi
#summary of the result by reading bed file.
#
#
#
