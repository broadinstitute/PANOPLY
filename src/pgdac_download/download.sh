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
  mkdir -p $dir
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
  mkdir -p $dir/rna; 
  cp cons_clust_tar/rna/*.pdf $dir/rna/.;
  mkdir -p $dir/cna; 
  cp cons_clust_tar/cna/*.png $dir/cna/.;
  mkdir -p $dir/sample-qc; 
  cp cons_clust_tar/sample-qc/*.pdf $dir/sample-qc/.;
  mkdir -p $dir/association; 
  cp cons_clust_tar/association/*.pdf $dir/association/.;
  mkdir -p $dir/clustering; 
  cp cons_clust_tar/clustering/*.pdf $dir/clustering/.; 
  cp cons_clust_tar/clustering/*.png $dir/clustering/.;
  
  mkdir -p $full/cons_clust_tar;
  cp -r cons_clust_tar/* $full/cons_clust_tar/.;
  rm -rf cons_clust_tar;
  cp -r ssgsea_ome $full/ssgsea_ome/;
  cp -r ssgsea_rna $full/ssgsea_rna/;
  cp -r $ssgsea_assoc $full/$ssgsea_assoc/;
  cp -r $ssgsea_clust $full/$ssgsea_clust/;

  mkdir -p $dir/ssgsea_ome;
  mkdir -p $dir/ssgsea_rna;
  mkdir -p $dir/$ssgsea_assoc;
  mkdir -p $dir/$ssgsea_clust;
  cp -r ssgsea_ome/* $dir/ssgsea_ome/.;
  cp -r ssgsea_rna/* $dir/ssgsea_rna/.;
  cp -r $ssgsea_assoc/* $dir/$ssgsea_assoc/.;
  cp -r $ssgsea_clust/* $dir/$ssgsea_clust/.;

  rm -rf ssgsea_ome ssgsea_rna $ssgsea_assoc $ssgsea_clust;

  if [[ ! -z $ptmsea ]]; then
    cp -r ptmsea $full/ptmsea/;
    mkdir -p $dir/ptmsea;
    cp -r ptmsea $dir/ptmsea/.;
    rm -rf ptmsea;
  fi
}

src=`pwd`
dir="$src/$analysis_dir-results"
full="$src/$analysis_dir-full-output"
dir_create
collect
cd $src;
tar -czvf pgdac_main_summary.tar.gz $analysis_dir-results
tar -czvf pgdac_main_full.tar.gz $analysis_dir-full-output
