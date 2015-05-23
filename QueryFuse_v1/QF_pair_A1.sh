#!/bin/bash

#Gene_specific_fusion_query subfunction: pocessing pair scenario A1 (pre-process, a1.1, a1.2)
#GSFQ_pair_A1.sh
if [ $# -ne 9 ]
then
  echo ""
    echo "Usage: QF_pair_A1.sh File_prefix BAM_file_folder_prefix result_folder_prefix whole_gene_list.bed human_genome.fa read_length LOG_ERROR QueryFuse_path Align_percent"
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
		echo ""
                exit 1
fi

if [ ! -d $1 ] 
then
  echo ""
  echo "Warning: The directory $1 does not generated in QF_pair_A1.sh, exit. in QF_pair_A1.sh, exit."
  echo ""
  exit 1
fi

if [ ! -x $3 ] 
then
	mkdir "$3"
fi

if [ ! -s $4 ] 
then
  echo ""
  echo "Warning: The file $4 does not exist in QF_pair_A1.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $5 ] 
then
  echo ""
  echo "Warning: The file $5 does not exist in QF_pair_A1.sh, exit."
  echo ""
  exit 1
fi

#Build the name for files
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


PAIR_BAM=$bam_fd"paired.bam"
PAIR_SORT=$file_prefix"paired_sorted"
PAIR_SORT_BAM=$bam_fd"paired_sorted.bam"
QUERY_BED=$outresult_fd"query_gene.bed"
QUERY_FA=$outresult_fd"query_gene.fa"
PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
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


#Separate spanning candidates. (just use grep,because number of span candidates should be really few)
sed '1,5d' $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL | cut -f10 |uniq -u > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID -a -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL ]
	then
		echo "Separate spanning candidates in QF_pair_A1.sh"
		$QF_path"grep_reads_from_psl_withheader_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL Rscript $QF_path
		#grep  -f $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL

	else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID" and/or "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL" is not generated in QF_pair_A1.sh, exit."
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID" and/or "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL" is not generated in QF_pair_A1.sh, exit." >> LOG_ERROR
                exit 1
fi


#If match >=read_length-23, they are reliable spanning read candidates. Get the information of mates:
#If match <read_length-23, they are unreliable split/span reads
#get IDs 
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL ]
	then
            #the psl file has no header already.
            sort -g $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL
            #echo "separate reads from a1 to a2 in QF_pair_A1.sh"
	else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL" is not generated in QF_pair_A1.sh, exit."
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL" is not generated in QF_pair_A1.sh, exit." >> LOG_ERROR
            exit 1
fi
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL ];
	then
		#echo "check row boundary in QF_pair_A1.sh"
                cut -f1 $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL > $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_MATCH
                Rscript $QF_path"boundary_row_ID.R" file.in=$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_MATCH lower_bound=0 upper_bound=$[READ_LEN-23] file.out=$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_ROW
                UP_ROW=`sed '2p' -n $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_ROW`
                ROW_NUM=`wc -l $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_MATCH | cut -d' ' -f1`
                #echo "the upper boundary is "$UP_ROW" should be same as the next number"
        else
                UP_ROW="NA"
                ROW_NUM="NA"   
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL" is not generated in QF_pair_A1.sh, exit."
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL" is not generated in QF_pair_A1.sh, exit." >> LOG_ERROR
                exit 1
fi

#if [ ! $UP_ROW ] 
#then
#  echo ""
#  echo "upper bound of QF_pair_A1.sh step does not available, exit this step."
#  echo ""
#  exit 1
#fi

echo "finished pre-processing of pair scenario a1"
echo "===================="
$QF_path"QF_pair_A1.1.sh" $PAIR_SORT $UP_ROW $ROW_NUM $GENE_BED $QUERY_BED $READ_LEN $QF_path $LOG_ERROR $Align_percent


echo "finished processing of pair scenario a1.1"
echo "===================="
$QF_path"QF_pair_A1.2.sh" $PAIR_SORT $UP_ROW $ROW_NUM $GENE_BED $QF_path $LOG_ERROR
echo "finished processing of pair scenario a1.2"
echo "===================="
