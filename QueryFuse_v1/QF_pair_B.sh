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

#Gene_specific_fusion_query subfunction: pocessing pair scenario B (pre-process, b1, b2)
#GSFQ_pair_B.sh
if [ $# -ne 11 ]
then
  echo ""
    echo "Usage: QF_pair_B.sh File_prefix BAM_file_folder_prefix result_folder_prefix whole_gene_list.bed human_genome.fa read_length tophat_genome_reference LOG_ERROR QF_path Align_percent size_query "
    echo ""
    echo "File_prefix - The folder that all the defualt can use (for intermedia files)."
    echo "BAM_file_folder_prefix - that has all three needed input bam files."
    echo "result_folder_prefix - folder location for final result."
    echo "whole_gene_list.bed - A bed file containing the whole gene list in the reference genome"
    echo "human_genome.fa - fa file of the whole human genome, such as hg19.fa"
    echo "read_length - length of reads"
    echo "tophat_genome_reference - indexed genome_reference."
    echo "LOG_ERROR - path to the error log file "
    echo "QueryFuse_path - QueryFuse path."
    echo "Align_percent"
    echo "size_query - step_size for blat to query"
    echo ""
    exit 1
fi

if [ ! -d $1 ] 
then
  echo ""
  echo "Warning: The directory $1 does not exist in QF_pair_B.sh, exit."
  echo ""
  exit 1
fi

if [ ! -x $3 ] 
then
	mkdir "$3"
fi

if [ ! -s $5 ] 
then
  echo ""
  echo "Warning: The file $5 does not exist in QF_pair_B.sh, exit."
  echo ""
  exit 1
fi

if [ ! -s $4 ] 
then
  echo ""
  echo "Warning: The file $4 does not exist in QF_pair_B.sh, exit."
  echo ""
  exit 1
fi


#name the parameters
file_prefix=$1
bam_fd=$2
outresult_fd=$3
GENE_BED=$4
HG_FA=$5
READ_LEN=$6
GENOME_REF=$7
LOG_ERROR=$8
QF_path=${9}
Align_percent=${10}
size_query=${11}
rep_match=$((1024*11/size_query))


#Build the name for files
PAIR_BAM=$bam_fd"paired.bam"
PAIR_SORT=$file_prefix"paired_sorted"
PAIR_SORT_BAM=$bam_fd"paired_sorted.bam"
QUERY_BED=$outresult_fd"query_gene.bed"
QUERY_FA=$outresult_fd"query_gene.fa"
PAIR_TO_QUERY_BAM=$PAIR_SORT"_to_query.bam"
PAIR_TO_QUERY_SAM=$PAIR_SORT"_to_query.sam" 
PAIR_TO_QUERY_BED=$PAIR_SORT"_to_query.bed"
PAIR_TO_QUERY_FILTER_BED=$PAIR_SORT"_to_query_filter.bed"



#pre-process pair B section
PAIR_TO_QUERY_FILTER_ID_DUPLI=$PAIR_SORT"_to_query_pair_filter_ID_duplicate.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_FA=$PAIR_SORT"_to_query_pair_filter_ID_duplicate.fa"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID.txt"
PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_SPLIT_ID_PSL=$PAIR_SORT"_to_query_pair_filter_ID_duplicate_on_query_split_ID.psl"

echo "pre-process scenario B in QF_pair_B.sh"

#get read ID both ends left
if [ -s $PAIR_TO_QUERY_FILTER_BED ];
	then
		echo `date`
		echo "get read ID both ends left in QF_pair_B.sh"
		cut -f4 $PAIR_TO_QUERY_FILTER_BED | cut -f1 | sort -u | cut -d'/' -f1 | uniq -d > $PAIR_TO_QUERY_FILTER_ID_DUPLI
	else
		echo "Warning: "$PAIR_TO_QUERY_FILTER_BED" is not exist in QF_pair_B.sh, exit"	
		echo "Warning: "$PAIR_TO_QUERY_FILTER_BED" is not exist in QF_pair_B.sh, exit"  >> $LOG_ERROR
                exit 1

fi

#get the duplicate ID reads and convert to fasta
if [ -s $PAIR_TO_QUERY_FILTER_ID_DUPLI ];
	then
		echo `date`
		echo "get the duplicate ID reads and convert to fasta in QF_pair_B.sh"
		PAIR_TO_QUERY_FILTER_ID_DUPLI_SAM=$PAIR_SORT"_to_query_filter_ID_dupl.sam"
                #PAIR_TO_QUERY_SAM is from QF_pair_A.sh
		$QF_path"grep_reads_from_SAM_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_DUPLI $PAIR_TO_QUERY_SAM $PAIR_TO_QUERY_FILTER_ID_DUPLI_SAM Rscript $QF_path
		#change the sam to fa
		cut -f'1,10' $PAIR_TO_QUERY_FILTER_ID_DUPLI_SAM | sort -u | awk '{OFS="\t"; print ">"$1"\n"$2}' > $PAIR_TO_QUERY_FILTER_ID_DUPLI_FA
		#samtools view $PAIR_TO_QUERY_BAM | grep -f $PAIR_TO_QUERY_FILTER_ID_DUPLI | cut -f'1,2,10' | sort -u | awk '{OFS="\t"; print ">"$FILE_DIR"\n"$GENE_ID"\n"}' > $PAIR_TO_QUERY_FILTER_ID_DUPLI_FA
	else
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI" is not exist in QF_pair_B.sh, exit"	
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI" is not exist in QF_pair_B.sh, exit"  >> $LOG_ERROR
                exit 1
fi


#blat the duplicate to query
if [ -s $QUERY_FA -a -s $PAIR_TO_QUERY_FILTER_ID_DUPLI_FA ];
	then
		echo `date`
		echo "blat the duplicate to query in QF_pair_B.sh"		
		blat -stepSize=$size_query -repMatch=$rep_match -minIdentity=$Align_percent $QUERY_FA $PAIR_TO_QUERY_FILTER_ID_DUPLI_FA $PAIR_TO_QUERY_FILTER_ID_DUPLI_ON_QUERY_PSL
	else
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_FA" is not exsit in QF_pair_B.sh, exit"
		echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_DUPLI_FA" is not exsit in QF_pair_B.sh, exit"  >> $LOG_ERROR
                exit 1
fi

echo "finished pre-processing of pair scenario b"
echo "===================="
$QF_path"QF_pair_B2.sh" $file_prefix $bam_fd $outresult_fd $GENE_BED $HG_FA $READ_LEN $LOG_ERROR $QF_path $Align_percent


$QF_path"QF_pair_B1.sh" $file_prefix $READ_LEN $QUERY_FA $GENE_BED $GENOME_REF $LOG_ERROR $QF_path $size_query
