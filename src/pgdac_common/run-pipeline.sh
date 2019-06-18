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
  echo "Usage: $0 OPERATION -i <input-tarball> -o <output-tarball>"
  echo "             -s <SM-output-file>  -a <parsed-data>"
  echo "             -n <normalized-data>  -f <filtered-data>"
  echo "             -e <expt-design-file> -r <analysis_directory>"
  echo "             -c <common-code-dir> -d <common-data-dir>"
  echo "             -p <parameters-file> -t <type>  -m <data>"
  echo "             -rna <rna-expression> -cna <copy-number-data>"
  echo "             -g <groups> -pe <#processors>"
  echo "             -CMAPgroup <CMAP-group-name>  -CMAPtype <CMAP-data-type>"
  echo "             -CMAPnperm <CMAP # permutations> -CMAPpmt <CMAP-permutation-scores-directory>"
  echo "             -CMAPscr <CMAP-subset-scores-directory> -CMAPcfg <CMAP-config-file>" 
  echo "   OPERATIONs allowed are listed below"
  echo "   <input-tarball> is a tar file with the initialized directory structure"
  echo "       if specfied, -s, -r, -c and -d are not required -- allowed only when OPERATION is not input*"
  echo "   <output-tarball> is the output tar file (defaults to OPERATION-output.tar)"
  echo "   <SM-output-file> is Spectrum Mill output in ssv format"
  echo "   <parsed-data> is parsed (SM or other) data file, before normalization/filtering"
  echo "   <normalized-data> is normalized input data table in gct2 or gct3 format"
  echo "   <expt-design-file> should contain Sample.ID, Experiment, Channel columns + optional annotation"
  echo "   <analysis_directory> is the directory name where the analysis is performed"
  echo "       if path is specified, only the basename is used; also used as tarfile name"
  echo "   <type> is proteome [default] or phosphoproteome"
  echo "   <data> is QC.pass [default] (or other subtype)"
  echo "   <rna-expression> GCT file with RNA expression data"
  echo "   <copy-number-data> GCT file with CNA data"
  echo "   <groups> is a expt-design-like-file used for subsetting CNA analysis"
  echo "   <#processors> is the number of processors (jobs) to use for CNA/correlation analysis"
  echo "   <CMAP-group-name> is group name for which CNA analysis was performed; defaults to 'all'"
  echo "   <CMAP-data-type> is data type for CMAP analysis; 'pome' or 'mrna'"
  echo "   <CMAP # permutations> is number of CMAP permutations to derive FDR"
  echo "   <CMAP-subset-scores-directory> is directory containing CMAP scores for data subsets run"
  echo "   <CMAP-permutation-scores-directory> is directory containing CMAP scores for permutation"
  echo "   <CMAP-config-file> is file with CMAP analysis parameters to over-ride defaults"
  echo "   Input Requirements:"
  echo "     OPERATION inputSM requires (-s, -e, -r, -c, -d)"
  echo "     OPERATION inputNorm requires (-n, -r, -c)"
  echo "     OPERATION normalize requires (-a, -r, -c) or (-i, -c); -i output from inputSM "
  echo "     OPERATION RNAcorr requires (-i, -c, -rna) OR (-f, -r, -c, -rna)"
  echo "     OPERATION harmonize requires (-i, -c, -d, -rna, -cna) OR (-r, -f, -c, -d, -rna, -cna)"
  echo "     OPERATION sampleQC requires (-i, -c); -i is tar output from harmonize"
  echo "     OPERATION CNAsetup requires (-i, -c), optional (-g, -pe); -i is tar output from harmonize"
  echo "     OPERATION CNAcorr requires (-i,); -i is tar output from CNAsetup"
  echo "     OPERATION CMAPsetup reqires (-i, -c); optional (-CMAPgroup, -CMAPtype, -CMAPcfg, -CMAPnperm); -i is tar output from CNAcorr"
  echo "     OPERATION CMAPconn reqires (-i, -CMAPscr); optional (-CMAPgroup, -CMAPtype, -CMAPcfg, -CMAPnperm, -CMAPpmt); -i is tar output from CMAPsetup"
  echo "     OPERATION assoc requires (-i, -c), optional (-g); or (-f, -r, -c, -g); -i is tar output from normalize/harmonize"
  echo "     OPERATION cluster required (-i, -c) OR (-f, -r, -c), optional (-g); -i is tar output from normalize/harmonize"
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
  # copy expt_file if specified (and keep original -- expt_design_file may get overwritten)
  if [ "$expt_file" != "" ]; then
    cp $expt_file $data_dir/$expt_design_file
    cp $expt_file $data_dir/.
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
  
  # for OPERATION = normalize
  if [ "$analysis" = "normalize" ]; then
    createSubdirs $parse_dir
    createSubdirs $norm_dir
    if [ ! -f "$parse_dir/$parsed_output" ]; then
      if [ "$parsed_data" = "" ]; then 
        echo "Parsed data not found ... abort"
        exit 1
      else
        cp $parsed_data $parse_dir/$parsed_output
      fi
    fi
    return
  fi
  
  # all other OPERATIONs
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
    CNAsetup )  if [ ! -f $harmonize_dir/rna-matrix.csv -o ! -f $harmonize_dir/cna-matrix.csv -o ! -f $harmonize_dir/$prefix-matrix.csv ]; then
                  echo "Harmonized data not found -- run harmonize first ... abort"
                  exit 1
                else
                  createSubdirs $cna_dir
                  # add to config.r if specified
                  if [ "$pe" != "" ]; then
                    echo "pe.max <- $pe" >> $cna_dir/config.r
                  fi
                  if [ "$groups" != "" ]; then
                    echo "cna.subgroups <- \"$groups\"" >> $cna_dir/config.r
                    cp $groups $cna_dir/.
                  fi
                fi ;;
    CNAcorr )   if [ ! -f $cna_dir/subgroups.txt -o ! -f $cna_dir/file_table.tsv -o ! -f $data_dir/gene-location.csv -o ! -f $data_dir/chr-length.csv ]; then
                  echo "Required initialization files not found -- run CNAsetup/harmonize first ... abort"
                  exit 1
                else
                  createSubdirs $cna_dir
                fi ;;
    CMAPsetup ) if [ ! -f $cna_dir/$cmap_prefix-matrix.csv -o ! -f $cna_dir/$cmap_prefix-vs-cna-sigevents.csv -o ! -f $cna_dir/$cmap_prefix-vs-cna-pval.csv -o ! -f $data_dir/cmap-knockdown-genes-list.txt ]; then
                  echo "Required initialization files not found -- run CNAcorr first ... abort"
                  exit 1
                else
                  createSubdirs $cmap_dir
                fi ;;
    CMAPconn ) if [ ! -d $cmap_dir ]; then
                  echo "CMAP directory not found -- run CMAPsetup first ... abort"
                  exit 1
                else
                  createSubdirs $cmap_dir
                fi ;;
    sampleQC )  if [ ! -f $harmonize_dir/rna-matrix.csv -o ! -f $harmonize_dir/cna-matrix.csv -o ! -f $harmonize_dir/$prefix-matrix.csv ]; then
                  echo "Harmonized data not found -- run harmonize first ... abort"
                  exit 1
                else
                  createSubdirs $qc_dir
                fi ;;
    assoc )     createSubdirs $assoc_dir
                # add to config.r if specified
                if [ "$groups" != "" ]; then
                    echo "assoc.subgroups <- \"$groups\"" >> $assoc_dir/config.r
                    cp $groups $assoc_dir/.
                fi ;;
    cluster )   createSubdirs $cluster_dir
                # add to config.r if specified 
                if [ "$groups" != "" ]; then
                    # class vectors to determine cluster enrichment
                    echo "cluster.enrichment.subgroups <- \"$groups\"" >> $cluster_dir/config.r
                    cp $groups $cluster_dir/.
                fi ;;
    * )
  esac
}



## set defaults
#  work with absolute paths for file/dir links/names
sm_file=
norm_data=
parsed_data=
filt_data=
expt_file=
input_tar=
output_tar=
analysis_dir=
code_dir=
common_data=
param_file=
data=
rna_data=
cna_data=
groups=
pe=
cmap_group="all"
cmap_type="pome"
cmap_scores=
cmap_nperm="0"
cmap_permutation="tmp"
cmap_config_file=
prefix="proteome"
data_source="default"
log_file="RUN-LOG.txt"
all_args=$@

op=$1
shift

## check $op is a supported operation
case $op in
  inputSM|inputNorm|normalize|RNAcorr|harmonize|CNAsetup|CNAcorr|CMAPsetup|CMAPconn|sampleQC|assoc|cluster) # OK
    ;;
  *)    echo "ERROR: Unknown OPERATION $op"; usage
        exit 1
esac

## read in arguments
while [ "$1" != "" ]; do
  case $1 in
	-a )     shift; parsed_data=`readlink -f $1` ;;
	-c )     shift; code_dir=`readlink -f $1` ;;
	-d )     shift; common_data=`readlink -f $1` ;;
	-e )     shift; expt_file=`readlink -f $1` ;;
	-f )     shift; filt_data=`readlink -f $1` ;;
	-g )     shift; groups=`readlink -f $1` ;;
	-i )     shift; input_tar=`readlink -f $1` ;;
	-o )     shift; output_tar=`readlink -f $1` ;;
	-m )     shift; data=$1 ;;
	-n )     shift; norm_data=`readlink -f $1` ;;
	-p )     shift; param_file=`readlink -f $1` ;;
	-r )     shift; analysis_dir=`readlink -f $1` ;;
	-s )     shift; sm_file=`readlink -f $1` ;;
	-t )     shift; prefix=$1 ;;
	-pe )    shift; pe=$1 ;;
	-rna )   shift; rna_data=`readlink -f $1` ;;
	-cna )   shift; cna_data=`readlink -f $1` ;;
	-CMAPgroup )
	         shift; cmap_group=$1 ;;
	-CMAPtype )
	         shift; cmap_type=$1 ;;
	-CMAPscr )
	         shift; cmap_scores=$1 ;;
	-CMAPnperm )
           shift; cmap_nperm=$1 ;;
  -CMAPpmt )
           shift; cmap_permutation=$1 ;;
	-CMAPcfg )
	         shift; cmap_config_file=`readlink -f $1` ;;
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
  normalize ) if [[ ("$input_tar" = "")  &&  ("$parsed_data" = "" || "$analysis_dir" = "") || "$code_dir" = "" ]]
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
  CNAcorr )   if [[ ("$input_tar" = "") ]]
              then
                usage
                exit 1
              fi ;;
  CMAPsetup ) if [[ ("$input_tar" = "") || "$code_dir" = "" ]]
              then
                usage
                exit 1
              fi ;;
  CMAPconn ) if [[ ("$input_tar" = "") || "$cmap_scores" = "" || ($cmap_nperm -gt 0 && "$cmap_permutation" = "") ]]
              then
                usage
                exit 1
              fi ;;
  sampleQC )  if [[ ("$input_tar" = "") || "$code_dir" = "" ]]
              then
                usage
                exit 1
              fi ;;
  assoc )     if [[ ("$input_tar" = "")  &&  ("$filt_data" = "" || "$analysis_dir" = "" || "$groups" = "")   \
                    || "$code_dir" = "" ]]
              then
                usage
                exit 1
              fi ;;
  cluster )   if [[ ("$input_tar" = "")  &&  ("$filt_data" = "" || "$analysis_dir" = "")   \
                    || "$code_dir" = "" ]]
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
cmap_dir="cmap"
cmap_prefix="$cmap_group-$cmap_type"
qc_dir="sample-qc"
assoc_dir="association"
cluster_dir="clustering"
if [ "$data" = "" ]; then
  subset_str=""
else 
  subset_str="-$data"
fi
expt_design_file="exptdesign.csv"
parsed_output="$prefix-ratio.gct"
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

                ## data preprocessing (parsed-data)
                (cd $parse_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" parseMSinput.r)
             ;;
  
#   inputNorm: input is a normalized gct (2 or 3) file that should just be filtered
    inputNorm ) cp $norm_data $norm_dir/$normalized_output
                for f in create-cls.r filter.r; do cp $code_dir/$f $norm_dir/$f; done
                
                ## filtering and cls file generation
                (cd $norm_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" filter.r;
                 R CMD BATCH --vanilla "--args $prefix $data" create-cls.r)
            ;;
#   normalize: start with parsed data (SM or other) and normalize/filter
    normalize ) analysisInit "normalize"
                for f in create-cls.r filter.r normalize.r; do cp $code_dir/$f $norm_dir/$f; done
               
                ## normalization (normalization)
                (cd $norm_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" normalize.r)
                 
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
                cp $cna_dir/subgroups.txt $cna_dir/file_table.tsv ../.
                for f in `cat $cna_dir/file_table.tsv`; do cp $cna_dir/$f ../.; done
             ;;
#   CNAcorr: run CNA analysis
#            input must be tar file obtained after CNAsetup
    CNAcorr )   analysisInit "CNAcorr"
                # this operation is only called from the command line for sequential execution
                # FireCloud module uses scatter/gather for parallelization, and does not call this operation
                for f in gene-location.csv chr-length.csv; do cp $data_dir/$f $cna_dir/$f; done
                (cd $cna_dir;
                 # read subgroups.txt into array
                 groups=`cat subgroups.txt`
                 g=($groups)
                 # iterate through groups and execute cna-analysis.r for each
                 for i in `seq 1 "${#g[@]}"`; do
                   line_i=`sed -n "${i}p" file_table.tsv`
                   l=($line_i)
                   group_prefix=${g[$(($i-1))]}
                   if [ ! -d ${group_prefix}-output ]; then
                     mkdir ${group_prefix}-output
                   fi
                   Rscript cna-analysis.r 1 1 $group_prefix "${l[0]}" "${l[1]}" "${l[2]}"
                   Rscript cna-analysis.r 0 1 $group_prefix NULL NULL NULL
                 done)
             ;;
#   CMAPsetup: run initialization for CMAP analysis -- sets up directories and creats input files (genesets)
#            input must be tar file obtained after CNAcorr
    CMAPsetup ) analysisInit "CMAPsetup"
                for f in cmap-input.r connectivity.r; do cp $code_dir/cmap-analysis/$f $cmap_dir/$f; done
                (cd $cmap_dir;
                 # generate cmap input genesets
                 Rscript cmap-input.r "../$cna_dir" "../$data_dir/cmap-knockdown-genes-list.txt" $cmap_group $cmap_type $cmap_nperm $cmap_config_file)
                # copy output (geneset file) to job wd (parent of $analysis_dir)
                cp $cmap_dir/$cmap_group-cmap-$cmap_type-updn-genes.gmt ../cmap-trans-genesets.gmt
                cp $cmap_dir/$cmap_group-cmap-$cmap_type-permuted-genes-*.gmt ../. 2>/dev/null || :    # permuted files may not always exist
             ;;
#   CMAPconn: calculate connectivity scores for CMAP analysis -- also gather ssGSEA scores from scatter (FireCloud) and assemble
#            input must be tar file obtained after CMAPsetup
    CMAPconn )  analysisInit "CMAPconn"
                mv ../$cmap_scores $cmap_dir/.
                if [ $cmap_nperm -gt 0 ]; then 
                  mv ../$cmap_permutation $cmap_dir/.
                fi
                (cd $cmap_dir;
                 # combine subset scores, and calculate connectivity scores
                 Rscript connectivity.r $cmap_scores $cmap_group $cmap_type $cmap_nperm $cmap_permutation $cmap_config_file)
             ;;
#   sampleQC: sample-level correlation between p'ome, RNA and CNA, with some clustering/fanplots
#             input must be tar file obtained after harmonize
    sampleQC )  analysisInit "sampleQC"
                for f in sample-qc.r; do cp $code_dir/$f $qc_dir/$f; done
                (cd $qc_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" sample-qc.r)
             ;;
#   assoc: association analysis for cls's in GCT or supplied in input
    assoc )     analysisInit "assoc"
                for f in assoc-analysis.r; do cp $code_dir/$f $assoc_dir/$f; done
                (cd $assoc_dir;
                 R CMD BATCH --vanilla "--args $prefix $data" assoc-analysis.r)
             ;;
#   cluster: perform consensus kmeans clustering
    cluster )   analysisInit "cluster"
                cp $code_dir/pgdac_kmeans_consensus.R $cluster_dir/
                cp $code_dir/consensus_clustering.R $cluster_dir/
                cp $code_dir/assoc-analysis.r $cluster_dir/
                cp $code_dir/postprocess.R $cluster_dir/
                (cd $cluster_dir;
                 
                 # extract 'label'  from 'analysis_dir'
                 label=`echo $analysis_dir | sed -E 's,.*/(.*).*,\1,'`
                 tmpdir='.' ## temp-folder
                  
                 ## extract sd threshold from 'config.r'   
                 sdclust=`cat config.r | grep -e '^clustering.sd.threshold' | awk -F' ' '{print $3}'`
                 
                 # run kmeans clustering and best cluster selection
                 # parmaters for minimal and maximal cluster numbers as well as 
                 # number of bootstrap iterations are fixed
                 Rscript pgdac_kmeans_consensus.R -i "${analysis_dir}" -u 2 -v 10 -b 1000 -s $sdclust -l $label -t $prefix -n $norm_dir -c $cluster_dir -d $tmpdir -z $code_dir
                 
                 # run association analysis on clusters to determine markers
                 R CMD BATCH --vanilla "--args $prefix $data" postprocess.R;
                 R CMD BATCH --vanilla "--args $prefix $data" assoc-analysis.r
             );;
#   Unknown operation
    * )         echo "ERROR: Unknown OPERATION $op"; exit 1 
esac

cd ..    # go back to parent of $analysis_dir (for tarball creation)

### Output tarball with $analysis_dir
# at this point, cwd should be parent of $analysis_dir
if [ "$output_tar" = "" ]; then
  tar -c -f $op-output.tar $analysis_dir
else
  tar -c -f $output_tar $analysis_dir
fi
# cp $analysis_dir.tar /results/.   # for testing


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

