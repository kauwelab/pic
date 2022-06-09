#! /bin/bash
dataSource=$1
destination=$2
par=$3
mpx=$4
outFile="${destination}${par}aboveR_${mpx}minPx.tsv"

for d in $dataSource/*
do
compName="`basename $d`"
gene1="`basename $d | cut -d'_' -f 1`"
gene2="`basename $d | cut -d'_' -f 2`"
file1="$d/$gene1$x$gene2" 
file2="$d/$gene2$x$gene1"
python mi_smart_filters.py -f1 "$file1" -f2 "$file2" -p1 $gene1 -p2 $gene2 -s $outFile -r $par -p $mpx -xi -xd
done 
