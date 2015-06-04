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

#extract reads in pair-b1 and a1 scenario.
#scan the sorted psl file(already sorted by read ID), compare max and min.
#should be read pair both aligned to query. One is fully(>90%) aligned and the other is partially(<len-23) aligned.  
#maybe not in the efficient enough way, but generally should have not that many candidates.

if [ $# -ne 5 ]
then
  echo ""
    echo "Usage: psl_to_ID_extract_split_major_on_query.sh psl_file output_name fa_file out_put_fa_file read_length"
  echo ""

  exit 1
fi

if [ ! -s $1 ] 
then
  echo ""
  echo "The file $1 does not exist"
  echo ""
  exit 1
fi

if [ ! -d logs ];
then
   mkdir logs
fi


if [ -f $2 ];
then
   rm $2
fi

if [ -f $4 ];
then
   rm $4
fi

if [ ! -s $3 ] 
then
  echo ""
  echo "The file $3 does not exist"
  echo ""
  exit 1
fi


## Start run for each file

READ_LEN=$5
PRE_ID=""
MIN_MATCH=0
MAX_MATCH=0

temp_file5=$1"temp_file5"
temp_file6=$1"temp_file6"

while read myfile	##read each line in the param file
do
	ID="$(echo $myfile | cut -d' ' -f10)"
	MATCH=$(echo $myfile | cut -d' ' -f1)
	if [ "$ID"x != "$PRE_ID"x ];
		then
			if [ $MAX_MATCH -ge $[READ_LEN *9/10] -a $MIN_MATCH -lt $[READ_LEN -23] ];
				then

			echo -e $PRE_ID >> $2
			grep -C1 $PRE_ID $3 >$temp_file5
			sed -n "/$PRE_ID/,$ p" $temp_file5 > $temp_file6
 			sed -i '1s/\($\)/\/1/g' $temp_file6
		  sed -i '3s/\($\)/\/2/g' $temp_file6
		  sed '' $temp_file6 >> $4
  		fi
			
			MIN_MATCH=$MATCH
			MAX_MATCH=$MATCH
			PRE_ID=$ID
		else
			if [ $MAX_MATCH -lt $MATCH ];
				then
					MAX_MATCH=$MATCH
			fi
			if [ $MIN_MATCH -gt $MATCH ];
				then
					MIN_MATCH=$MATCH
			fi
	fi
done < $1 ##very important! the $2 means the input file for the while loop.

if [ $MAX_MATCH -ge $[READ_LEN *9/10] -a $MIN_MATCH -lt $[READ_LEN -23] ];
		then

			echo -e $PRE_ID >> $2
			grep -C1 $PRE_ID $3 >$temp_file5
			sed -n "/$PRE_ID/,$ p" $temp_file5 > $temp_file6
 			sed -i '1s/\($\)/\/1/g' $temp_file6
		  sed -i '3s/\($\)/\/2/g' $temp_file6
		  sed '' $temp_file6 > $4
fi

if [ -s $temp_file5 ];
	then
		rm $temp_file5
fi

if [ -s $temp_file6 ];
	then
		rm $temp_file6
fi