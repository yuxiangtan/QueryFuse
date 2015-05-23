#!/bin/bash

#Gene_specific_fusion_query subfunction: psl_to_unmap pocessing (extract candidates)
#GSFQ_psl_to_unmap_process.sh

if [ $# -ne 3 ]
then
  echo ""
    echo "Warning: Usage: QF_psl_to_unmap_process.sh unmap_query_psl READ_LEN R_script QF_path python_script"
    echo ""
    echo "unmap_query_psl - location of the unmapped_first/second_mate_over_query.psl."
    echo "READ_LEN - length of reads"
    echo "QF_path"
    echo ""

  exit 1
fi

if [ ! -s $1 ] 
then
  echo ""
  echo "Warning: The file $1 does not exist in QF_psl_to_unmap_process.sh, exit."
  echo ""
  exit 1
fi


#Build the name for files
READ_LEN=$2
QF_path=$3
UNMAP_ON_QUERY_SPLIT=$1"_split"
UNMAP_ON_QUERY_SPAN_MATE_QUERY=$1"_span_mate_query"
UNMAP_ON_QUERY_SPAN_MATE_NOT_QUERY=$1"_span_mate_not_query"
UNMAP_ON_QUERY_SPLIT_ID=$1"_split_ID.txt"
UNMAP_ON_QUERY_SPLIT_ID_SUBTRACT_BED=$1"_split_ID_subtract.bed"
UNMAP_ON_QUERY_SPAN_MATE_QUERY_ID=$1"_span_mate_query_ID.txt"
UNMAP_ON_QUERY_SPAN_MATE_NOT_QUERY_ID=$1"_span_mate_not_query_ID.txt"
UNMAP_ON_QUERY_SPLIT_ID_SUBTRACT_ID=$1"_split_ID_subtract_ID.txt"


echo "separate unmap_query psl into candidate groups in QF_psl_to_unmap_process.sh"

python $QF_path"psl_separator_by_aligned_len_no_filter.py" -i $1 -a $UNMAP_ON_QUERY_SPAN_MATE_QUERY -b $UNMAP_ON_QUERY_SPLIT -c $UNMAP_ON_QUERY_SPAN_MATE_NOT_QUERY -n 12 -N $[READ_LEN-23] -j 5
Rscript $QF_path"subtract_psl.R" file.in=$UNMAP_ON_QUERY_SPLIT read_length=$READ_LEN file.out=$UNMAP_ON_QUERY_SPLIT_ID_SUBTRACT_BED

#in case some of the ID has "/" in it, so use the second cut to make the format uniform.
cut -f10 $UNMAP_ON_QUERY_SPAN_MATE_QUERY | cut -f1 -d"/" > $UNMAP_ON_QUERY_SPAN_MATE_QUERY_ID
cut -f10 $UNMAP_ON_QUERY_SPAN_MATE_NOT_QUERY | cut -f1 -d"/" > $UNMAP_ON_QUERY_SPAN_MATE_NOT_QUERY_ID

cut -f1 $UNMAP_ON_QUERY_SPLIT_ID_SUBTRACT_BED | cut -f1 -d"/" > $UNMAP_ON_QUERY_SPLIT_ID_SUBTRACT_ID
