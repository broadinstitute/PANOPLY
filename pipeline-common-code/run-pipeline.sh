#!/bin/bash

# $Id$

##
## Setup and run analysis pipeline
##


#
# 1. Setup directory structure and create links to code and data
#

function usage
{
  echo "Usage: $0 -s <SM-output-file> -n <norm-data> "
  echo "             -r <directory> -c <common-code-dir> -d <common-data-dir>"
  echo "             [-l <previous-pipeline-dir>] [-t <type>]  [-m <data>]  [-u | -x]  [-v]"
  echo "       <SM-output-file> is Spectrum Mill output in ssv format"
  echo "       <norm-data> is normalized input data table"
  echo "       <type> is proteome [default] or phosphoproteome"
  echo "       <data> is unimodal [default] (or other subtype)"
  echo "       <previous-pipeline-dir> is a directory where pipeline was previously run"
  echo "           if -l is specified, data-with-intensity and 2comp-normalization are"
  echo "           symlinked to previous run"
  echo "       -x enables use of SLURM (default: do not use SLURM)"
  echo "       -u enables use of UGER (default: do not use UGER)"
  echo "           if both -s and -u are specified, the option specified last will take precedence"
  echo "       -v enabels running cna plots (default: do not run cna)"
  echo "       Parameters in [ ... ] are optional. Rest are required."
  echo "       Use -h to print this message."
}


# set defaults
# work with absolute paths for file/dir links/names
sm_file=
norm_data=
analysis_dir=
code_dir=
data_dir=
prev_dir=
prefix="proteome"
data=
use_cluster="none"
cluster_queue=
scratch_dir=
run_cna="FALSE"
data_source="default"
log_file="RUN-LOG.txt"
all_args=$@


while [ "$1" != "" ]; do
    case $1 in
	-s )     shift
	         sm_file=`readlink -f $1`
		 ;;
	-n )     shift
	         norm_data=`readlink -f $1`
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
	-l )     shift
	         prev_dir=`readlink -f $1`
		 ;;
	-t )     shift
	         prefix=$1
		 ;;
	-m )     shift
	         data=$1
		 ;;
	-x )     use_cluster="slurm"
		 ;;
	-u )     use_cluster="uger"
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


if [[ ( "$sm_file" = "" && "$norm_data" = "" ) || "$analysis_dir" = "" || "$code_dir" = "" || "$data_dir" = "" ]]
then
  usage
  exit 1
fi

if [[ "$use_cluster" = "slurm" ]]
then
  cluster_queue="${CLUSTER_QUEUE:?Set CLUSTER_QUEUE environment variable when using SLURM}"
  scratch_dir="${SCRATCH_DIR:?Set SCRATCH_DIR environment variable to point to scratch file system}"
fi


## create directory where all files are put
if [ -d $analysis_dir ]
then
  echo "Directory $analysis_dir exists. To avoid confusion, specify a new directory"
  exit
fi
mkdir $analysis_dir
cd $analysis_dir

## operations below in $analysis_dir

# create run log
cat <<EOT >> $log_file
CPTAC/PGDAC data analysis pipeline
Command line: $0 $all_args
Analysis directory: $analysis_dir
Normalized data file: $norm_data
Common code directory: $code_dir
Common data directory: $data_dir
Previous data version: $prev_dir

Data directory contents:
`ls -l $data_dir`

EOT


echo "Creating directories and links ..." >> $log_file

## data directory
mkdir data
for f in `ls $data_dir`
do
  ln -s $data_dir/$f data/$f
done


## create appropriate directories for analysis
sub_dirs="2comp-normalization netgestalt cna rna-seq mutation association correlation"
if [ "$sm_file" != "" ]
then
  sub_dirs="data-with-intensity $sub_dirs"
fi


## copy preamble -- this could change from run to run
cp $code_dir/preamble.r preamble.r
echo "compute.cluster.type <- '$use_cluster'" >> preamble.r
echo "cluster.queue <- '$cluster_queue'" >> preamble.r   # will be used only for SLURM
echo "scratch.fs <- sprintf ('%s/pid-%s/', '$scratch_dir', Sys.getpid())" >> preamble.r   # will be used only for SLURM
for d in $sub_dirs
do
  mkdir $d
  (cd $d; ln -s ../preamble.r preamble.r)
done


## link to previous datasets
if [ "$prev_dir" != "" ]
then
  for nf in `ls -1 $prev_dir/2comp-normalization/*.gct | sed '/modal-/d'`
  do
    ln -s $nf 2comp-normalization/`basename $nf`
  done
  for f in create-cls.r; do ln -s $code_dir/$f 2comp-normalization/$f; done
else
  if [ "$norm_data" != "" ]
  then
    ln -s $norm_data 2comp-normalization/$prefix-ratio-norm.gct
    for f in create-cls.r filter.r; do ln -s $code_dir/$f 2comp-normalization/$f; done
  else  # $sm_data != ""
    for f in parseMSinput.r; do ln -s $code_dir/$f data-with-intensity/$f; done
    for f in create-cls.r filter.r normalize.r; do ln -s $code_dir/$f 2comp-normalization/$f; done
  fi
fi



## link appropriate code in preparation for running pipeline components
for f in rna-seq.r rna-seq-correlation.r; do ln -s $code_dir/$f rna-seq/$f; done
for f in cna-analysis.r generate-cna-plots.r run-cna-analysis.r run-cna-analysis-all.r runR-$use_cluster.sh; do ln -s $code_dir/$f cna/$f; done
for f in netgestalt-input-data.r; do ln -s $code_dir/$f netgestalt/$f; done
for f in mut-analysis.r sigmut-genes.r mutated-gene-markers.r; do ln -s $code_dir/$f mutation/$f; done
for f in assoc-analysis.r; do ln -s $code_dir/$f association/$f; done
for f in correlation.r run-correlation.r runR-$use_cluster.sh; do ln -s $code_dir/$f correlation/$f; done

# link the following files since they are used by a LSF array job script
# (and is easier to manage as files in the current directory)
for f in gene-location.csv chr-length.csv; do (cd cna; ln -s ../data/$f $f); done
for f in gene-location.csv; do (cd correlation; ln -s ../data/$f $f); done


#
# 2. Run the various components
#    (some components are run synchronously with -K to address data dependencies)
#


# decide whether to use GridEngine or not
cmd_sync=""
cmd_async=""
jname="job$BASHPID"
if [ $use_cluster = "uger" ]
then
  cmd_sync="qsub -sync y -q short -o uger.out -cwd -j y -N $jname"
  cmd_async="qsub -q short -o uger.out -cwd -j y -N $jname"
  if [ $prefix = "phosphoproteome" ]
  then
    cmd_sync="$cmd_sync -l h_vmem=12g $code_dir/runR-uger.sh"
    cmd_async="$cmd_async -l h_vmem=12g $code_dir/runR-uger.sh"
  else
    cmd_sync="$cmd_sync -l h_vmem=8g $code_dir/runR-uger.sh"
    cmd_async="$cmd_async -l h_vmem=8g $code_dir/runR-uger.sh"
   fi
fi

if [ $use_cluster = "slurm" ]
then
  # request 16G memory in all cases
  cmd_sync="salloc -p $cluster_queue -t 5:0:0 --mem=12G -J $jname /usr/bin/srun $code_dir/runR-slurm.sh"
  cmd_async="sbatch -p $cluster_queue -t 5:0:0 --mem=12G -o slurm.out -J $jname $code_dir/runR-slurm.sh"
fi



## data preprocessing (data-with-intensity) -- do not run if linked to previous version
if [ "$prev_dir" = "" -a "$norm_data" = "" ]
then
  echo "Data preprocessing ..." >> $log_file
  (cd data-with-intensity;
   $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" parseMSinput.r)
fi


## normalization (2comp-normalization)
if [ "$prev_dir" = "" -a "$norm_data" = "" ]
then
  echo "Normalization & Filtering ..." >> $log_file
  (cd 2comp-normalization;
   $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" normalize.r;
   $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" filter.r)
else
  if [ "$norm_data" != "" ]
  then
    echo "Filtering ..." >> $log_file
    (cd 2comp-normalization;
     $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" filter.r)
  fi
fi
# run only create-cls.r if linked to previous version (also if normalized data provided)
(cd 2comp-normalization; 
 $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" create-cls.r)


## RNA-seq
if [ $data_source != "whim" ]
then
  echo "RNA-seq correlation ..." >> $log_file
  (cd rna-seq;
    $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" rna-seq.r;
    $cmd_async R CMD BATCH --vanilla "--args $prefix $data" rna-seq-correlation.r)
fi


## association
echo "Association analysis ..." >> $log_file
(cd association;
 $cmd_async R CMD BATCH --vanilla "--args $prefix $data" assoc-analysis.r)


# ## mutated gene marker analysis
# echo "Mutated gene marker analysis ..." >> $log_file
# (cd mutation;
#  $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" mutated-gene-markers.r)


# ## mutation analysis (generate tables needed by netgestalt)
# if [ $data_source != "whim" ]
# then
#   echo "Mutation analysis ..." >> $log_file
#   (cd mutation;
#    $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" mut-analysis.r)
# fi
  

## netgestalt
echo "Generating tables for NetGestalt ..." >> $log_file
(cd netgestalt;
 $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" netgestalt-input-data.r;
 ln -s rnaseq-matrix.csv mrna-matrix.csv)
  
  
## correlation (needs tables from netgestalt)
echo "Correlation analysis ..." >> $log_file
(cd correlation;
 $cmd_async R CMD BATCH --vanilla "--args $prefix $data" run-correlation.r)

  
# ## additional mutation analysis (uses global proteome table from netgestalt)
# if [ $data_source != "whim" ]
# then
#   echo "More mutation analysis ..." >> $log_file
#   (cd mutation;
#    $cmd_sync R CMD BATCH --vanilla "--args $prefix $data" sigmut-genes.r)
# fi
  

if [ $run_cna = "TRUE" ]
then
## cna (plots only -- relies on tables created by netgestalt)
echo "Generating CNA plots ..." >> $log_file
(cd cna;
 $cmd_async R CMD BATCH --vanilla "--args $prefix $data" run-cna-analysis-all.r)
fi


