#!/bin/bash

##usage: ./outlier.sh inputdata.tar output.tar

E_NOARGS=85
E_NOFILE=86
E_EXISTFILE=87

if [ -z "$2" ]
	then
		echo "Usage `basename $0` inputfile outfile"
		exit $E_NOARGS
fi

if [ ! -f "$1" ]
	then
		echo "$1 not found"
		exit $E_NOFILE
fi

if [ -f "$2" ]
	then
		echo "$2 already exists"
		exit $E_EXISTSFILE
fi

#inputfolder=$(pwd)/input
#outputfolder=$(pwd)/output
if [ ! -d "/input/data" ]
	then
	mkdir -p /input/data
fi 
if [ ! -d "/output/res" ]
	then
	mkdir -p /output/res
fi 
tar -C /input/data -xf $1 

Rscript /prog/get_outlier.R

##compress output result into tar file 
tar -cf $2 -C /output/res . 

##clean up
rm -rf /input/data
rm -rf /output/res 


