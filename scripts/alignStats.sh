#!/bin/bash

if [ $# != 1]; then
	echo "USAGE: alignStats.sh <qsub directory>"
	exit 1;
else
	
	cd $1
	for i in `ls *out`; do echo -e $(ls ${i} | cut -f1 -d".")"\t"$(head -1 ${i} | cut -f2 -d":" | sed 's/^ //')"\t"$(head -2 $i | tail -1 | cut -f2 -d":" | sed 's/^ //' | cut -f1 -d" ")"\t"$(head -2 $i | tail -1 | cut -f2 -d":" | sed 's/^ //' | cut -f2 -d" " | sed 's/[()]//g')"\t"$(head -3 $i | tail -1 | cut -f2 -d":" | sed 's/^ //' | cut -f1 -d" ")"\t"$(head -3 $i | tail -1 | cut -f2 -d":" | sed 's/^ //' | cut -f2 -d" " | sed 's/[()]//g'); done > alignStats.txt

fi
