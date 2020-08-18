#!/bin/bash
while getopts ":c:o:r:a:s:l:p:" opt; do
    case $opt in
        c) cons_clust_tar="$OPTARG";;
        o) ssgsea_ome="$OPTARG";;
        r) ssgsea_rna="$OPTARG";;
	      a) analysis_dir="$OPTARG";;
        s) ssgsea_assoc="$OPTARG";;
        l) ssgsea_clust="$OPTARG";;
        p) ptmsea="$OPTARG";;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

scatter_processing()
{
  dir_name=$1
  array=($(ls $dir_name/*.tar))
  for index in "${!array[@]}"
  do
    mkdir -p $dir_name/$index && tar xf ${array[index]} -C $dir_name/$index;
    group="$( echo -e "$( grep "input gct" $dir_name/$index/*parameters.txt | \
      cut -d ':' -f 2 | tr -d '[:space:]')" | \
      rev | cut -d '.' -f2- | rev )"
    mkdir -p $dir_name/"ssgsea-$group";
    cp -r $dir_name/$index/* $dir_name/"ssgsea-$group"/.;
    rm -rf $dir_name/$index;
  done
}

dir_create()
{
  cd $src;
  mkdir -p $summ_path
  mkdir -p cons_clust_tar && tar xf $cons_clust_tar -C cons_clust_tar --strip-components 1
  mkdir -p ssgsea_ome && tar xf $ssgsea_ome -C ssgsea_ome
  mkdir -p ssgsea_rna && tar xf $ssgsea_rna -C ssgsea_rna
  scatter_processing $src/$ssgsea_assoc
  scatter_processing $src/$ssgsea_clust
  if [[ ! -z $ptmsea ]]; then
    mkdir -p ptmsea && tar xf $ptmsea -C ptmsea
  fi
}


collect()
{
  cd $src;
  mkdir -p $summ_path/rna; 
  cp cons_clust_tar/rna/*.pdf $summ_path/rna/.;
  mkdir -p $summ_path/cna; 
  cp cons_clust_tar/cna/*.png $summ_path/cna/.;
  mkdir -p $summ_path/sample-qc; 
  cp cons_clust_tar/sample-qc/*.pdf $summ_path/sample-qc/.;
  mkdir -p $summ_path/association; 
  cp cons_clust_tar/association/*.pdf $summ_path/association/.;
  mkdir -p $summ_path/clustering; 
  cp cons_clust_tar/clustering/*.pdf $summ_path/clustering/.; 
  cp cons_clust_tar/clustering/*.png $summ_path/clustering/.;
  
  mkdir -p $full_path;
  cp -r cons_clust_tar/* $full_path/.;
  rm -rf cons_clust_tar;
  cp -r ssgsea_ome $full_path/ssgsea_ome/;
  cp -r ssgsea_rna $full_path/ssgsea_rna/;
  cp -r $ssgsea_assoc $full_path/$ssgsea_assoc/;
  cp -r $ssgsea_clust $full_path/$ssgsea_clust/;

  mkdir -p $summ_path/ssgsea_ome;
  mkdir -p $summ_path/ssgsea_rna;
  mkdir -p $summ_path/$ssgsea_assoc;
  mkdir -p $summ_path/$ssgsea_clust;
  cp -r ssgsea_ome/* $summ_path/ssgsea_ome/.;
  cp -r ssgsea_rna/* $summ_path/ssgsea_rna/.;
  cp -r $ssgsea_assoc/* $summ_path/$ssgsea_assoc/.;
  cp -r $ssgsea_clust/* $summ_path/$ssgsea_clust/.;

  rm -rf ssgsea_ome ssgsea_rna $ssgsea_assoc $ssgsea_clust;

  if [[ ! -z $ptmsea ]]; then
    cp -r ptmsea $full_path/ptmsea/;
    mkdir -p $summ_path/ptmsea;
    cp -r ptmsea $summ_path/ptmsea/.;
    rm -rf ptmsea;
  fi
}

src=`pwd`

summ="$analysis_dir-summary-results"
full="$analysis_dir-full-results"

summ_path=$src/$summ
full_path=$src/$full
dir_create
collect
cd $src;
tar -cvf panoply_main_summary.tar $summ
tar -cvf panoply_main_full.tar $full
