#!/bin/bash



I="$1"
N="$2"
OFFSET="0"
if [ "$#" == "3" ]; then
    OFFSET="$3"
elif [ ! "$#" == "2" ]; then
    echo "Error: Incorrect number of arguments" >&2
    echo "" >&2
    echo "Typical syntax" >&2
    echo "" >&2
    echo "  select_structure.sh i n [offset] < traj.raw > structure_10.raw" >&2
    echo "  i=the structure you want from your file (indexing begins at 1)" >&2
    echo "  n=the number of lines per structure (INCLUDING BLANK LINE SEPARATORS)" >&2
    echo "  offset (optional) is the number of lines at the beginning of the file to skip" >&2
    echo "" >&2
    echo "example:" >&2
    echo "  select_structure.sh 10 4001 < traj.raw > structure_10.raw" >&2
    exit 1
fi


LINE_START="$((1+OFFSET+(I-1)*N))"

#echo "N=$N"
#echo "LINE_START=$LINE_START"
tail -n "+$LINE_START" | head -n "$N"
