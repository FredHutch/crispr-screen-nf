#!/bin/bash

if [ $# != 5 ]; then
        echo "USAGE: mageck.sh <counts_table> <treatment cols comma-sep> <control cols comma-sep> <output_dir> <sgRNA controls (NTC)>"
        exit 1;
else
	counts=$1
	treat=$2
	control=$3
	outDir=$4
	sgCtrl=$5
	outBase=$(basename $counts .txt)

	ml texlive/20200624
	cd $outDir
	#mkdir -p $outBase
	#cd $outBase
	/home/pchanana/miniconda3/bin/mageck test -k ${counts} -t ${treat} -c ${control} -n ${outBase}_median --pdf-report --normcounts-to-file
	/home/pchanana/miniconda3/bin/mageck test -k ${counts} -t ${treat} -c ${control} -n ${outBase}_control --pdf-report --normcounts-to-file --norm-method control --control-sgrna ${sgCtrl} 
	#/home/pchanana/miniconda3/bin/mageck test -k ${counts} -t ${treat} -c ${control} -n ${outBase}_total --pdf-report --normcounts-to-file --norm-method total

fi
