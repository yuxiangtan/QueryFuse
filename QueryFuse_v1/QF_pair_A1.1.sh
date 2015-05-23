#!/bin/bash

#Gene_specific_fusion_query subfunction: pocessing pair scenario A1.1
#GSFQ_pair_A1.1.sh
if [ $# -ne 9 ]
then
  echo ""
    echo "Usage: QF_pair_A1.1.sh PAIR_SORT UP_ROW ROW_NUM whole_gene_list.bed QUERY_BED read_length QF_path LOG_ERROR Align_percent"
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


#name the parameters
PAIR_SORT=$1
UP_ROW=$2
ROW_NUM=$3
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
PAIR_TO_QUERY_FILTER_ID_UNIQ=$PAIR_SORT"_to_query_filter_ID_uniq.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_FA=$PAIR_SORT"_to_query_filter_ID_uniq.fa"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_ID.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_MATCH=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_match"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_ROW=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_row"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_up.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP_GOOD=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_up_good.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_ID.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_FILTER=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_ID.txtfilter.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_NO_HEADER_SAM=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_ID_extract_no_header.sam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_ID_extract.sam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_BAM=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_ID_extract.bam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_gene_anno.bed"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED_PAIRED=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed"

if [ -n "$UP_ROW" ];
	then
		if [ "$UP_ROW" != "NA" ];
			then
			echo `date`
			echo "upper boundary is: "$UP_ROW
			echo 'get reliable span candidate in QF_pair_A1.1.sh'
			sed "${UP_ROW},${ROW_NUM}p" -n $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP
			#to reduce false positive support reads, control the maximum mismatch to at most 1- alignment percentage
			python $QF_path"eliminate_mismatch_read_psl.py" $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP $READ_LEN $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP_GOOD $Align_percent
			cut -f10 $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP_GOOD > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD
		fi
	else
		echo "Warning: Upper bound row number is not available in QF_pair_A1.1.sh, exit."
		echo "Warning: Upper bound row number is not available in QF_pair_A1.1.sh, exit." >> $LOG_ERROR
                exit 1
fi

#if there are reads fulfil the requirement
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD ];
    then
#extract alignment information from bam file using ID
#echo "use ID list to filter bam"
#perl read_ID_to_bamfilter.pl ./ $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD ./
#echo ""
#
#bamtools filter -in $PAIR_TO_QUERY_BAM -out $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_BAM -script $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_FILTER
        #echo `date`
        #echo "bam to sam in QF_pair_A1.1.sh"
        #sam is generated in QF_pair_A.sh already
        #samtools view $PAIR_TO_QUERY_BAM > $PAIR_TO_QUERY_SAM
        #echo `date`
        echo "grep_reads from sam in QF_pair_A1.1.sh"
#echo "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD\n"
#echo "$PAIR_TO_QUERY_SAM\n"
#echo "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_NO_HEADER_SAM\n"
        $QF_path"grep_reads_from_SAM_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD $PAIR_TO_QUERY_SAM $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_NO_HEADER_SAM Rscript $QF_path
        #echo `date`
        #echo "add header in QF_pair_A1.1.sh"
        samtools view -H $PAIR_TO_QUERY_BAM | cat > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM
        cat $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_NO_HEADER_SAM >> $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM
        #echo "sam to bam in QF_pair_A1.1.sh"
    else
        echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD" is not generated in QF_pair_A1.1.sh, exit."
	echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD" is not generated in QF_pair_A1.1.sh, exit." >> $LOG_ERROR
        exit 1
fi

if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM ];
    then
        samtools view -S $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM -b > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_BAM
#echo `date`
        rm $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_NO_HEADER_SAM
        echo "get gene information by intersectBed in QF_pair_A1.1.sh"
        intersectBed -abam $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_BAM -b $GENE_BED -bed -wo > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED
#summary of the result by reading bed file.
#
#
    else
	echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM" is not generated in QF_pair_A1.1.sh, exit."
	echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_SAM" is not generated in QF_pair_A1.1.sh, exit." >> $LOG_ERROR
        exit 1
fi

PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED_PAIRED_TEMP=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_good_gene_anno_paired.bed_TEMP"
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED ];
        then
                perl $QF_path"QF_pair_A1.1_group_to_pair.pl" $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED_PAIRED_TEMP $QUERY_BED $READ_LEN $Align_percent
                sort -u $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED_PAIRED_TEMP > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED_PAIRED
                rm $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED_PAIRED_TEMP
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED" is not generated in QF_pair_A1.1.sh"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID_GOOD_EXTRACT_ANNO_BED" is not generated in QF_pair_A1.1.sh"  >> $LOG_ERROR
fi
