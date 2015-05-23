#!/bin/bash

#Gene_specific_fusion_query subfunction: pocessing pair scenario A (pre-process, a1, a2)
#GSFQ_pair_A.sh


if [ $# -ne 11 ]
then
  echo ""
    echo "Usage: QF_pair_A.sh File_prefix BAM_file_folder_prefix result_folder_prefix whole_gene_list.bed human_genome.fa read_length tophat_genome_reference LOG_ERROR QF_path Align_percent size_query"
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


#File_prefix
if [ ! -d $1 ] 
then
  echo ""
  echo "Warning: The directory $1 does not exist in QF_pair_A.sh, exit."
  echo ""
  exit 1
fi
#result_folder_prefix
if [ ! -x $3 ] 
then
	mkdir "$3"
fi

#whole_gene_list
if [ ! -s $4 ] 
then
  echo ""
  echo "Warning: The file $4 does not exist in QF_pair_A.sh, exit."
  echo ""
  exit 1
fi
#human_genome
if [ ! -s $5 ] 
then
  echo ""
  echo "Warning: The file $5 does not exist in QF_pair_A.sh, exit."
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
PAIR_TO_QUERY_FILTER_ID_UNIQ=$PAIR_SORT"_to_query_filter_ID_uniq.txt"
PAIR_TO_QUERY_FILTER_ID_SAM=$PAIR_SORT"_to_query_filter_ID.sam"
PAIR_TO_QUERY_FILTER_ID_UNIQ_FA=$PAIR_SORT"_to_query_filter_ID_uniq.fa"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_ID=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_ID.txt"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp.psl"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_MATCH=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_match"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_ROW=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_row"
PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL_SPAN_TEMP_PSL_UP=$PAIR_SORT"_to_query_filter_ID_uniq_on_query_span_temp_up.psl"

echo "pre-process scenario A"
echo `date`
echo "get read ID has only one end left"
cut -f4  $PAIR_TO_QUERY_FILTER_BED | cut -f1 | sort -u | cut -d'/' -f1 | uniq -u > $PAIR_TO_QUERY_FILTER_ID_UNIQ
#get the unique ID reads and convert to fasta
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ ];
    then
            echo "convert unique ID reads to fasta file in QF_pair_A.sh"
            #samtools view $PAIR_TO_QUERY_BAM > $PAIR_TO_QUERY_SAM #moved to pair_end_process because many steps will call this file.
            $QF_path"grep_reads_from_SAM_by_ID.sh" $PAIR_TO_QUERY_FILTER_ID_UNIQ $PAIR_TO_QUERY_SAM $PAIR_TO_QUERY_FILTER_ID_SAM Rscript $QF_path
            cut -f'1,10' $PAIR_TO_QUERY_FILTER_ID_SAM | sort -u | awk '{OFS="\t"; print ">"$1"\n"$2}' > $PAIR_TO_QUERY_FILTER_ID_UNIQ_FA
    else
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ" is not generated in QF_pair_A.sh, exit" >> $LOG_ERROR
            echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ" is not generated in QF_pair_A.sh, exit"
            exit 1
fi

#blat on query again for unique ID: (each ID can have one or two output, one means spanning, two means splitting)
if [ -s $PAIR_TO_QUERY_FILTER_ID_UNIQ_FA ];
        then
                echo "blat on query again in QF_pair_A.sh"
                blat -stepSize=$size_query -repMatch=$rep_match -minIdentity=$Align_percent $QUERY_FA $PAIR_TO_QUERY_FILTER_ID_UNIQ_FA $PAIR_TO_QUERY_FILTER_ID_UNIQ_ON_QUERY_PSL
        else
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_FA" is not generated in QF_pair_A.sh, exit"
                echo "Warning: "$PAIR_TO_QUERY_FILTER_ID_UNIQ_FA" is not generated in QF_pair_A.sh, exit" >> $LOG_ERROR
                exit 1
fi

echo "finished pre-processing of pair scenario a"
echo "===================="
#echo "file_prefix"$file_prefix
#echo "bam_fd"$bam_fd
#echo "outresult_fd"$outresult_fd
#echo "GENE_BED"$GENE_BED
#echo "HG_FA"$HG_FA
#echo "READ_LEN"$READ_LEN
#echo "Scc"$SCC
#echo "LOG_ERROR"$LOG_ERROR
#echo "R_script"$R_script
#echo "QF_path"$QF_path
#echo "Perl_script"$Perl_script
#echo "Python_script"$Python_script
#echo "Align_percent"$Align_percent
$QF_path"QF_pair_A1.sh" $file_prefix $bam_fd $outresult_fd $GENE_BED $HG_FA $READ_LEN $LOG_ERROR $QF_path $Align_percent

#echo $GENOME_REF
$QF_path"QF_pair_A2.sh" $file_prefix $READ_LEN $QUERY_FA $GENE_BED $GENOME_REF $LOG_ERROR $QF_path $size_query
echo "finished processing of pair scenario a2"
echo "===================="