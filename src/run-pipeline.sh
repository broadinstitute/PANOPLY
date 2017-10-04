#!/bin/bash

### Broad Institute PGDAC PIPELINE
### Script for invoking various pipeline components

##
## Setup and run analysis pipeline
##
## Expected usage:
##  First called with OPERATION = inputSM or inputNorm
##   - this will create normalized and filtered datasets ready for analysis
##   - output: filtered data and tarball
##   - experiment design/class labels are embedded in the normalized/filtered data table (gct3)
##  Next, call each analysis module OPERATION = <analysis>, as needed
##   - using either input tarball or (config.r, filtered.data) will create
##   - appropriate analysis directory, copy code and perform analysis
##  Each call creates a tarball from initialized directory structure (or input config/data files)
##   to pass as input to subsequent OPERATIONs
##


## Utility functions

function usage {
  echo "Usage: $0 OPERATION -i <input-tarball>"
  echo "             -s <SM-output-file>  -n <normalized-data>  -f <filtered-data>"
  echo "             -e <expt-design-file> -r <analysis_directory>"
  echo "             -c <common-code-dir> -d <common-data-dir>"
  echo "             -p <parameters-file> -t <type>  -m <data>"
  echo "             -rna <rna-expression> -cna <copy-number-data>"
  echo "             -g <groups> -pe <#processors>"
  echo "   OPERATION is one of inputSM, inputNorm, RNAcorr"
  echo "   <input-tarball> is a tar file with the initialized directory structure"
  echo "       if specfied, -s, -r, -c and -d are not required -- allowed only when OPERATION is not init"
  echo "   <SM-output-file> is Spectrum Mill output in ssv format"
  echo "   <normalized-data> is normalized input data table in gct2 or gct3 format"
  echo "   <expt-design-file> should contain Sample.ID, Experiment, Channel columns + optional annotation"
  echo "   <analysis_directory> is the directory name where the analysis is performed"
  echo "       if path is specified, only the basename is used; also used as tarfile name"
  echo "   <type> is proteome [default] or phosphoproteome"
  echo "   <data> is unimodal [default] (or other subtype)"
  echo "   <rna-expression> GCT file with RNA expression data"
  echo "   <copy-number-data> GCT file with CNA data"
  echo "   <groups> is a expt-design-like-file used for subsetting CNA analysis"
  echo "   <pe> is the number of processors (jobs) to use for CNA/correlation analysis"
  echo "   Input Requirements:"
  echo "     OPERATION inputSM requires (-s, -e, -r, -c, -d)"
  echo "     OPERATION inputNorm requires (-n, -r, -c)"
  echo "     OPERATION RNAcorr requires (-i, -c, -rna) OR (-f, -r, -c, -rna)"
  echo "     OPERATION harmonize requires (-i, -c, -d, -rna, -cna) OR (-r, -f, -c, -d, -rna, -cna)"
  echo "     OPERATION CNAsetup requires (-i, -c), optional (-g, -pe)"
  echo "   Use -h to print this message."
}


function runlog {
  # create run log
  cat <<EOT >> $1
PGDAC data analysis pipeline
Command line: $0 $all_args
Analysis directory: $analysis_dir
Input SM data file: $sm_file
Input nomalized data: $norm_data
Input filtered data: $filt_data
Input tar file: $input_tar
Code directory: $code_dir
Common data directory: $common_data
Parameters file: $param_file
Data prefix: $prefix
Data type: $data
Data source: $data_source
EOT
}


function createAnalysisDir {
  # create analysis directory, if needed
  analysis_dir=`basename $1`   # ignore path
  if [ ! -d $analysis_dir ]; then
    mkdir $analysis_dir
  fi
}


function initDataDir {
  # create and initialize data directory (in analysis_dir)
  if [ ! -d $data_dir ]; then
    mkdir $data_dir
  fi
  
  if [ "$common_data" != "" ]; then   # if $common_data is specified, copy contents
    for f in `ls $common_data`; do
      if [ ! -f $data_dir/$f ]; then
        cp $common_data/$f $data_dir/$f
      fi
    done  
  fi
  # copy expt_file if specified
  if [ "$expt_file" != "" ]; then
    cp $expt_file $data_dir/$expt_design_file
  fi
}


function createConfig {
  # copy config and concat parameters
  if [ "$param_file" = "" ]; then
    cp $code_dir/config.r config.r
  else
    cat $code_dir/config.r $param_file > config.r
  fi
}


function createSubdirs {
  # create subdirectories (if they do not exisit) and copy config file
  for d in $1; do
    if [ ! -d $d ]; then
      mkdir $d
    fi
    (cd $d; cp ../config.r config.r)
  done
}


function analysisInit {
  ## initialize directories and check data files for specified analysis
  analysis=$1
  if [ ! -d $norm_dir ]; then
    echo "Directory not found: $analysis_dir/$norm_dir ... abort"
    exit 1
  fi
  if [ "$filt_data" = "" ]; then
    if [ ! -f $norm_dir/$filtered_output ]; then
      echo "Filtered dataset not found: $norm_dir/$filtered_output ... abort"
      exit 1
    fi
  else
    cp $filt_data $norm_dir/$filtered_output
  fi
  
  case $analysis in
    RNAcorr )   if [ "$rna_data" = "" ]; then
                  echo "RNA expression data not found ... abort"
                  exit 1
                else
                  createSubdirs $rna_dir
                  cp $rna_data $data_dir/$rna_data_file
                fi ;;
    harmonize ) if [ "$rna_data" = "" -o "$cna_data" = "" ]; then
                  echo "RNA expression or CNA data not found ... abort"
                  exit 1
                else
                  createSubdirs $harmonize_dir
                  cp $rna_data $data_dir/$rna_data_file
                  cp $cna_data $data_dir/$cna_data_file
                fi ;;
    CNAsetup ) if [ ! -f $harmonize_dir/rna-matrix.csv -o ! -f $harmonize_dir/cna-matrix.csv -o ! -f $harmonize_dir/$prefix-matrix.csv ]; then
                  echo "Harmonized data not found -- run harmonize first ... abort"
                  exit 1
                else
                  createSubdirs $cna_dir
                  # add to config.r if specified
                  if [ "$pe" != "" ]; then
                    echo "pe.max <- $pe" >> $cna_dir/config.r
                  fi
                  if [ "$groups" != "" ]; then
                    echo "cna.subgroups <- $groups" >> $cna_dir/config.r
                  fi
                fi ;;
    * )
  esac
}



## set defaults
#  work with absolute paths for file/dir links/names
sm_file=
norm_data=
filt_data=
expt_file=
input_tar=
analysis_dir=
code_dir=
common_data=
param_file=
data=
rna_data=
cna_data=
groups=
pe=
prefix="proteome"
data_source="default"
log_file="RUN-LOG.txt"
all_args=$@

op=$1
shift

## check $op is a supported operation
case $op in
  inputSM|inputNorm|RNAcorr|harmonize|CNAsetup) # OK
    ;;
  *)    echo "ERROR: Unknown OPERATION $op"; usage
        exit 1
esac

## read in arguments
while [ "$1" != "" ]; do
  case $1 in
	-c )     shift; code_dir=`readlink -f $1` ;;
	-d )     shift; common_data=`readlink -f $1` ;;
	-e )     shift; expt_file=`readlink -f $1` ;;
	-f )     shift; filt_data=`readlink -f $1` ;;
	-g )     shift; groups=$1 ;;
	-i )     shift; input_tar=`readlink -f $1` ;;
	-m )     shift; data=$1 ;;
	-n )     shift; norm_data=`readlink -f $1` ;;
	-p )     shift; param_file=`readlink -f $1` ;;
	-r )     shift; analysis_dir=`readlink -f $1` ;;
	-s )     shift; sm_file=`readlink -f $1` ;;
	-t )     shift; prefix=$1 ;;
	-pe )    shift; pe=$1 ;;
	-rna )   shift; rna_data=`readlink -f $1` ;;
	-cna )   shift; cna_data=`readlink -f $1` ;;
  -h )     usage; exit ;;
	* )      usage; exit 1
  esac
  shift
done


## check appropriate parameters have been provided
case $op in 
  inputSM )   if [[ "$sm_file" = "" || "$expt_file" = "" || "$analysis_dir" = ""   \
                    || "$code_dir" = "" || "$common_data" = "" ]]
              then
                usage
                exit 1
              fi ;;
  inputNorm ) if [[ "$norm_data" = "" || "$analysis_dir" = "" || "$code_dir" = "" ]]
              then
                usage
                exit 1
              fi ;;
  RNAcorr )   if [[ ("$input_tar" = "")  &&  ("$filt_data" = "" || "$analysis_dir" = "")   \
                    || "$code_dir" = "" || "$rna_data" = "" ]]
              then
                usage
                exit 1
              fi ;;
  harmonize ) if [[ ("$input_tar" = "")  &&  ("$filt_data" = "" || "$analysis_dir" = "")   \
                    || "$code_dir" = "" || "$common_data" = "" || "$rna_data" = "" || "$cna_data" = "" ]] 
              then
                usage
                exit 1
              fi ;;
  CNAsetup )  if [[ ("$input_tar" = "") || "$code_dir" = "" ]]
              then
                usage
                exit 1
              fi ;;
  * )        usage           # unknown operation
             exit 1
esac
  

## some definitions to more easily coordinate with config.r
data_dir="data"
parse_dir="parsed-data"
norm_dir="normalized-data"
harmonize_dir="harmonized-data"
rna_dir="rna"
cna_dir="cna"
if [ "$data" = "" ]; then
  subset_str=""
else 
  subset_str="-$data"
fi
expt_design_file="exptdesign.csv"
normalized_output="$prefix-ratio-norm.gct"
filtered_output="$prefix-ratio-norm-NArm$subset_str.gct"
rna_data_file="rna-data.gct"
cna_data_file="cna-data.gct"


## INITIALIZATION 
## Directory setup and/or extract tarball
if [ $op = "inputSM" -o $op = "inputNorm" -o "$input_tar" = "" ]
then
  ### input tar file not specified, or ignored (when $op=inputSM or inputNorm)
  ## create directory where all files are put
  createAnalysisDir $analysis_dir
  cd $analysis_dir
  
  ## operations below in $analysis_dir
  runlog $log_file
  ## data directory
  initDataDir
  ## copy config and concat parameters -- this could change from run to run
  createConfig
  createSubdirs $norm_dir   # always needed for any OPERATION
else   
  ### tar file specified
  ## extract tarball in current directory and set $analysis_dir
  tar -x -f $input_tar
  analysis_dir=`tar -t -f $input_tar | head -1 | sed -e 's/\/.*//'`
  cd $analysis_dir
  ## data directory (in case previously not specified)
  initDataDir
fi


## At this point, CWD is $analysis_dir

## OPERATIONs
## Perform initialization and processing for each $op type:
#   - copy input data to keep everything together; and
#   - copy appropriate code in preparation for running pipeline components
#   - run various code/components
case $op in 
#   inputSM: input is a SpectrumMill ssv file that should be parsed, normalized and filtered
    inputSM )   createSubdirs $parse_dir
                cp $sm_file $data_dir/$prefix-SMout.ssv
                cp $expt_file $data_dir/$expt_design_file
                for f in parseMSinput.r; do cp $code_dir/$f $parse_dir/$f; done
                for f in create-cls.r filter.r normalize.r; do cp $code_dir/$f $norm_dir/$f; done
               
                ## data preprocessing (parsed-data)
                (cd $parse_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" parseMSinput.r)
 
                ## normalization (normalization)
                (cd $norm_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" normalize.r)
                 
                ## filtering and cls file generation
                (cd $norm_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" filter.r;
                 R CMD BATCH --vanilla "--args $prefix $data" create-cls.r)
            ;;
  
#   inputNorm: input is a normalized gct (2 or 3) file that should just be filtered
    inputNorm ) cp $norm_data $norm_dir/$normalized_output
                for f in create-cls.r filter.r; do cp $code_dir/$f $norm_dir/$f; done
                
                ## filtering and cls file generation
                (cd $norm_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" filter.r;
                 R CMD BATCH --vanilla "--args $prefix $data" create-cls.r)
            ;;
#   RNAcorr: RNA-seq (or microarray) expression correlation with proteome
    RNAcorr )   analysisInit "RNAcorr"
                for f in rna-seq.r rna-seq-correlation.r; do cp $code_dir/$f $rna_dir/$f; done
                (cd $rna_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" rna-seq.r;
                 R CMD BATCH --vanilla "--args $prefix $data" rna-seq-correlation.r)
             ;;
#   harmonize: Harmonize RNA, CNA and proteome data to create gene-centric tables with common
#              rows (genes) and columns (samples)
    harmonize ) analysisInit "harmonize"
                for f in harmonize.r; do cp $code_dir/$f $harmonize_dir/$f; done
                (cd $harmonize_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" harmonize.r)
             ;;
#   CNAsetup: setup directories and code for running CNA analysis
#             input must be tar file obtained after harmonize
    CNAsetup )  analysisInit "CNAsetup"
                for f in cna-analysis.r cna-analysis-setup.r generate-cna-plots.r; do cp $code_dir/$f $cna_dir/$f; done
                (cd $cna_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" cna-analysis-setup.r)
                # copy required outputs to job wd (parent of $analysis_dir)
                cp $cna_dir/jids.txt $cna_dir/subgroups.txt $cna_dir/file_table.tsv ../.
                for f in `cat $cna_dir/file_table.tsv`; do cp $cna_dir/$f ../.; done
             ;;
#   Unknown operation
    * )         echo "ERROR: Unknown OPERATION $op"; exit 1 
esac

cd ..    # go back to parent of $analysis_dir (for tarball creation)

### Output tarball with $analysis_dir
# at this point, cwd should be parent of $analysis_dir
tar -c -f $op-output.tar $analysis_dir
# cp $analysis_dir.tar /results/.   # for testing


# ## copy appropriate code in preparation for running pipeline analysis components
# for f in cna-analysis.r generate-cna-plots.r run-cna-analysis.r run-cna-analysis-all.r runR.sh; do cp $code_dir/$f cna/$f; done
# for f in netgestalt-input-data.r; do cp $code_dir/$f netgestalt/$f; done
# for f in assoc-analysis.r; do cp $code_dir/$f association/$f; done
# for f in correlation.r run-correlation.r runR.sh; do cp $code_dir/$f correlation/$f; done
# for f in mut-analysis.r sigmut-genes.r mutated-gene-markers.r; do cp $code_dir/$f mutation/$f; done
# 
# # link the following files since they are used by a LSF array job script
# # (and is easier to manage as files in the current directory)
# for f in gene-location.csv chr-length.csv; do (cd cna; cp ../data/$f $f); done
# for f in gene-location.csv; do (cd correlation; cp ../data/$f $f); done


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



