#!/bin/bash
set -e # exit upon error condition
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

while getopts ":t:o:r:a:s:p:n:f:m:e:c:" opt; do
    case $opt in
        t) association_tar="$OPTARG";;
        o) ssgsea_ome="$OPTARG";;
        r) ssgsea_rna="$OPTARG";;
        a) analysis_dir="$OPTARG";;
        s) ssgsea_assoc="$OPTARG";;
        p) ptmsea="$OPTARG";;
        n) nmf_res="$OPTARG";;
        f) nmf_figs="$OPTARG";;
        m) nmf_ssgsea="$OPTARG";;
        e) omicsev_tar="$OPTARG";;
        c) cosmo_tar="$OPTARG";;
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
  mkdir -p association_tar && tar xf $association_tar -C association_tar --strip-components 1
  mkdir -p omicsev_tar && tar xf $omicsev_tar -C omicsev_tar
  mkdir -p cosmo_tar && tar xf $cosmo_tar -C cosmo_tar
  mkdir -p ssgsea_ome && tar xf $ssgsea_ome -C ssgsea_ome
  if [[ ! -z $ssgsea_rna ]]; then
    mkdir -p ssgsea_rna && tar xf $ssgsea_rna -C ssgsea_rna
  fi
  scatter_processing $src/$ssgsea_assoc
  if [[ ! -z $ptmsea ]]; then
    mkdir -p ptmsea && tar xf $ptmsea -C ptmsea
  fi
  if [ ! -z $nmf_res ] && [ ! -z $nmf_figs ]; then
    mkdir -p so_nmf
    mkdir -p so_nmf/nmf_results && tar xf $nmf_res -C so_nmf/nmf_results
    mkdir -p so_nmf/nmf_figures && tar xf $nmf_figs -C so_nmf/nmf_figures
    mkdir -p so_nmf/nmf_ssgsea && tar xf $nmf_ssgsea -C so_nmf/nmf_ssgsea
  fi
}


collect()
{
  cd $src;
  mkdir -p $summ_path/rna; 
  cp association_tar/rna/*.pdf $summ_path/rna/.;
  mkdir -p $summ_path/cna; 
  cp association_tar/cna/*.png $summ_path/cna/.;
  mkdir -p $summ_path/sample-qc; 
  cp association_tar/sample-qc/*.pdf $summ_path/sample-qc/.;
  mkdir -p $summ_path/association; 
  cp association_tar/association/*.pdf $summ_path/association/.;
  
  mkdir -p $full_path;
  cp -r association_tar/* $full_path/.;
  rm -rf association_tar;
  cp -r omicsev_tar/* $full_path/.;
  rm -rf omicsev_tar;
  cp -r cosmo_tar/* $full_path/.;
  rm -rf cosmo_tar;
  cp -r ssgsea_ome $full_path/ssgsea_ome/;
  cp -r $ssgsea_assoc $full_path/$ssgsea_assoc/;
  if [[ ! -z $ssgsea_rna ]]; then
    cp -r ssgsea_rna $full_path/ssgsea_rna/;
  fi

  mkdir -p $summ_path/ssgsea_ome;
  mkdir -p $summ_path/$ssgsea_assoc;
  cp -r ssgsea_ome/* $summ_path/ssgsea_ome/.;
  cp -r $ssgsea_assoc/* $summ_path/$ssgsea_assoc/.;
  rm -rf ssgsea_ome $ssgsea_assoc;
  
  if [[ ! -z $ssgsea_rna ]]; then
    mkdir -p $summ_path/ssgsea_rna;
    cp -r ssgsea_rna/* $summ_path/ssgsea_rna/.;
    rm -rf ssgsea_rna;
  fi

  if [[ ! -z $ptmsea ]]; then
    cp -r ptmsea $full_path/ptmsea/;
    mkdir -p $summ_path/ptmsea;
    cp -r ptmsea $summ_path/ptmsea/.;
    rm -rf ptmsea;
  fi

  if [ ! -z $nmf_res ] && [ ! -z $nmf_figs ]; then
    cp -r so_nmf $full_path/so_nmf/;
    mkdir -p $summ_path/so_nmf;
    cp -r so_nmf/nmf_figures/ $summ_path/so_nmf; # copy only the figures to summary-path
    rm -rf so_nmf;
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
