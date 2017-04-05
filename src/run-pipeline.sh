#!/bin/bash

# $Id$

##
## Setup and run analysis pipeline
##
## Expected usage:
##  First called with OPERATION = init
##  Create tarball from initialized directory structure and pass as input to
##  all other OPERATIONs
##
##


#
# 1. Setup directory structure and create links to code and data
#

function usage
{
  echo "Usage: $0 OPERATION -i <input-tarball>"
  echo "             -s <SM-output-file> -e <expt-design-file> -r <analysis_directory>"
  echo "             -c <common-code-dir> [-d <common-data-dir>]"
  echo "             [-p <parameters-file>] [-t <type>]  [-m <data>]  [-v]"
  echo "   OPERATION is one of init, parseSM, normalize, RNAcorr"
  echo "   <input-tarball> is a tar file with the initialized directory structure"
  echo "       if specfied, -s, -r, -c and -d are not required -- allowed only when OPERATION is not init"
  echo "   <SM-output-file> is Spectrum Mill output in ssv format"
  echo "   <expt-design-file> should contain Sample.ID, Experiment, Channel columns + optional annotation"
  echo "   <analysis_directory> is the directory name where the analysis is performed"
  echo "       if path is specified, only the basename is used; also used as tarfile name"
  echo "   <type> is proteome [default] or phosphoproteome"
  echo "   <data> is unimodal [default] (or other subtype)"
  echo "   -v enabels running cna plots (default: do not run cna)"
  echo "   Parameters in [ ... ] are optional. Rest are required."
  echo "   Use -h to print this message."
}


# set defaults
# work with absolute paths for file/dir links/names
sm_file=
expt_file=
input_tar=
analysis_dir=
code_dir=
data_dir=
param_file=
prefix="proteome"
data="unimodal"
run_cna="FALSE"
data_source="default"
log_file="RUN-LOG.txt"
all_args=$@

op=$1
shift

# check $op is a supported operation
case $op in
  init|parseSM|norm|RNAcorr) # OK
    ;;
  *)    usage
        exit 1
esac

# read in arguments
while [ "$1" != "" ]; do
    case $1 in
	-s )     shift
	         sm_file=`readlink -f $1`
		 ;;
	-e )     shift
	         expt_file=`readlink -f $1`
		 ;;
	-i )     shift
	         input_tar=`readlink -f $1`
		 ;;
	-r )     shift
	         analysis_dir=`readlink -f $1`
		 ;;
	-c )     shift
	         code_dir=`readlink -f $1`
		 ;;
	-d )     shift
	         data_dir=`readlink -f $1`
		 ;;
	-p )     shift
	         param_file=`readlink -f $1`
	   ;;
	-t )     shift
	         prefix=$1
		 ;;
	-m )     shift
	         data=$1
		 ;;
	-v )     run_cna="TRUE"
	         ;;
        -h )     usage
	         exit
		 ;;
	* )      usage
	         exit 1
    esac
    shift
done


if [ $op = "init" ]
then
  if [ "$sm_file" = "" -o "$expt_file" = "" -o "$analysis_dir" = "" -o "$code_dir" = "" ]
  then
    usage
    exit 1
  fi
else
  if [ "$input_tar" = "" ]
  then
    usage
    exit 1
  fi
fi



### OPERATION init
if [ $op = "init" ]
then
  
  ## create directory where all files are put
  if [ -d $analysis_dir ]
  then
    echo "Directory $analysis_dir exists. To avoid confusion, specify a new directory"
    exit
  fi
  analysis_dir=`basename $analysis_dir`   # ignore path
  mkdir $analysis_dir
  cd $analysis_dir
  
  ## operations below in $analysis_dir
  
  # create run log
  cat <<EOT >> $log_file
PGDAC data analysis pipeline
Command line: $0 $all_args
Analysis directory: $analysis_dir
Input SM data file: $sm_file
Code directory: $code_dir
Input data directory: $data_dir

Data directory contents:
`ls -l $data_dir`
EOT
  
  
  echo "Creating directories and initializing ..." >> $log_file
  
  ## data directory
  mkdir data
  # copy input data to keep everything together
  cp $sm_file data/$prefix-SMout.ssv
  cp $expt_file data/exptdesign.csv
  if [ "$data_dir" != "" ]   # if $data_dir is specified, copy contents
  then
    for f in `ls $data_dir`
    do
      cp $data_dir/$f data/$f
    done
  fi
  
  
  ## create appropriate directories for analysis
  sub_dirs="normalization parsed-data rna-seq"
  
  ## copy config and concat parameters -- this could change from run to run
  if [ "$param_file" = "" ]
  then
    cp $code_dir/config.r config.r
  else
    cat $code_dir/config.r $param_file > config.r
  fi
  

  for d in $sub_dirs
  do
    mkdir $d
    (cd $d; cp ../config.r config.r)
  done
  
  
  ## copy appropriate code in preparation for running pipeline components
  for f in parseMSinput.r; do cp $code_dir/$f parsed-data/$f; done
  for f in create-cls.r filter.r normalize.r; do cp $code_dir/$f normalization/$f; done
  for f in rna-seq.r rna-seq-correlation.r; do cp $code_dir/$f rna-seq/$f; done
  # for f in cna-analysis.r generate-cna-plots.r run-cna-analysis.r run-cna-analysis-all.r runR.sh; do cp $code_dir/$f cna/$f; done
  # for f in netgestalt-input-data.r; do cp $code_dir/$f netgestalt/$f; done
  # for f in mut-analysis.r sigmut-genes.r mutated-gene-markers.r; do cp $code_dir/$f mutation/$f; done
  # for f in assoc-analysis.r; do cp $code_dir/$f association/$f; done
  # for f in correlation.r run-correlation.r runR.sh; do cp $code_dir/$f correlation/$f; done
  # 
  # # link the following files since they are used by a LSF array job script
  # # (and is easier to manage as files in the current directory)
  # for f in gene-location.csv chr-length.csv; do (cd cna; cp ../data/$f $f); done
  # for f in gene-location.csv; do (cd correlation; cp ../data/$f $f); done

  cd ..    # go back to parent of $analysis_dir
else
  # extract tarball in current directory and set $analysis_dir
  tar -x -f $input_tar
  analysis_dir=`tar -t -f $input_tar | head -1 | sed -e 's/\/.*//'`
fi   # OPERATION init

#
# 2. Run the various components
#


### OPERATION parseSM
if [ $op = "parseSM" ]
then
  
  ## data preprocessing (parsed-data)
  echo "Data preprocessing ..." >> $log_file
  (cd "$analysis_dir/parsed-data";
   R CMD BATCH --vanilla "--args $prefix $data" parseMSinput.r)

fi   # OPERATION parseSM



### OPERATION norm
if [ $op = "norm" ]
then

  ## normalization (normalization)
  echo "Normalization ..." >> $log_file
  (cd "$analysis_dir/normalization";
   R CMD BATCH --vanilla "--args $prefix $data" normalize.r;
   R CMD BATCH --vanilla "--args $prefix $data" filter.r)
  
  # (cd normalization; R CMD BATCH --vanilla "--args $prefix $data" create-cls.r)
  
fi   # OPERATION norm



### OPERATION RNAcorr
if [ $op = "RNAcorr" ]
then

  ## rna-seq
  echo "RNAseq correlation ..." >> $log_file
  (cd "$analysis_dir/rna-seq";
   R CMD BATCH --vanilla "--args $prefix $data" rna-seq.r;
   R CMD BATCH --vanilla "--args $prefix $data" rna-seq-correlation.r)
  
fi    # OPERATION RNAcorr

# ## association
# echo "Association analysis ..." >> $log_file
# (cd association;
#  R CMD BATCH --vanilla "--args $prefix $data" assoc-analysis.r)
# 
# 
# # ## mutated gene marker analysis
# # echo "Mutated gene marker analysis ..." >> $log_file
# # (cd mutation;
# #  R CMD BATCH --vanilla "--args $prefix $data" mutated-gene-markers.r)
# 
# 
# # ## mutation analysis (generate tables needed by netgestalt)
# # if [ $data_source != "whim" ]
# # then
# #   echo "Mutation analysis ..." >> $log_file
# #   (cd mutation;
# #    R CMD BATCH --vanilla "--args $prefix $data" mut-analysis.r)
# # fi
#   
# 
# ## netgestalt
# echo "Generating tables for NetGestalt ..." >> $log_file
# (cd netgestalt;
#  R CMD BATCH --vanilla "--args $prefix $data" netgestalt-input-data.r)
#   
#   
# ## correlation (needs tables from netgestalt)
# echo "Correlation analysis ..." >> $log_file
# (cd correlation;
#  R CMD BATCH --vanilla "--args $prefix $data" run-correlation.r)
# 
#   
# # ## additional mutation analysis (uses global proteome table from netgestalt)
# # if [ $data_source != "whim" ]
# # then
# #   echo "More mutation analysis ..." >> $log_file
# #   (cd mutation;
# #    R CMD BATCH --vanilla "--args $prefix $data" sigmut-genes.r)
# # fi
#   
# 
# if [ $run_cna = "TRUE" ]
# then
# ## cna (plots only -- relies on tables created by netgestalt)
# echo "Generating CNA plots ..." >> $log_file
# (cd cna;
#  R CMD BATCH --vanilla "--args $prefix $data" run-cna-analysis-all.r)
# fi


### Output tarball with $analysis_dir
# at this point, cwd should be parent of $analysis_dir
tar -c -f $op-output.tar $analysis_dir
# cp $analysis_dir.tar /results/.   # for testing

