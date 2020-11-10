#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

## Release a new version of PANOPLY
##  (must be invoked from the release directory)

panoply=`pwd`/..
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}Error:${reg}"
warn="${grn}Warning:${reg}" ## notification

## documentation location (assumes wiki repo is available in path)
doc_dir="$panoply/../PANOPLY.wiki"


display_usage() {
  echo "usage: release.sh -v <VER> [-f] [-l] [-u] [-h]"
  echo "       -v <VER>  specify version string (required)"
  echo "       -f  force update of exisiting production workspaces"
  echo "       -l  rebuild panoply_libs docker from scratch"
  echo "       -u  rebuild panoply_utils docker from scratch"
  echo "       -h  print usage"
  exit
}

while getopts "v:fluh" opt; do
    case $opt in
        v) VER="$OPTARG";;
        f) f_flag="true";;
        l) l_flag="true";;
        u) u_flag="true";;
        h) display_usage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

if [[ -z $VER ]]; then
  echo -e "$err Release version not specified. Exiting.."
  exit 1
fi

# ensure that VER has only characters, numbers, _ and -
if ! [[ "$VER" =~ ^[a-zA-Z0-9_\-]+$ ]]; then
  echo -e "$err Invalid version. Must be alphanumeric with dash/unscore"
  exit 1
fi

# default development and release namespaces
DEV=broadcptacdev
REL=broadcptac
proj=broad-firecloud-cptac
prod="PANOPLY_Production_Pipelines_v$VER"
prod_all="PANOPLY_Production_Modules_v$VER"

# Before beginning, check production workspaces do not exist
if ! fissfc space_exists -w $prod_all -p $proj -q; then
  if [[ $f_flag == "true" ]]; then
    echo -e "$warn Workspace $prod_all exists. Existing methods will be replaced."
  else
    echo -e "$err Workspace $prod_all exists. Delete workspace or use -f option."
    exit 1
  fi
fi
if ! fissfc space_exists -w $prod -p $proj -q; then
  if [[ $f_flag == "true" ]]; then
    echo -e "$warn Workspace $prod exists. Existing methods will be replaced."
  else
    echo -e "$err Workspace $prod exists. Delete workspace or use -f option."
    exit 1
  fi
fi

# Check if documentation (from Wiki) has been checked out
if [ ! -d $doc_dir ]; then
  echo -e "$err Documentation not found in $doc_dir. Checkout wiki from GitHub."
  exit 1
fi


## ** Get latest code version in dev branch, create release branch
cd $panoply
git checkout dev
git pull
# create new branch for release
rel_id="release-$VER"
git checkout -b $rel_id dev
git push origin $rel_id


## ** Rebuild/update all docker images in $DEV
# Build / update panoply_libs
cd $panoply/hydrant
if [[ $l_flag == "true" ]]; then
  ./setup.sh -t panoply_libs -n $DEV -y -b -u 
fi

# Build panoply_utils and panoply_common
if [[ $u_flag == "true" ]]; then
  ./setup.sh -t panoply_utils -n $DEV -y -b -u 
fi
./setup.sh -t panoply_common -n $DEV -y -b -u -x

# Rebuild all other task dockers
modules=( $( ls -d $panoply/hydrant/tasks/panoply_* | xargs -n 1 basename |  \
  sed '/panoply_libs/d' | sed '/panoply_utils/d' | sed '/panoply_common/d') )
for mod in "${modules[@]}"
do
  ./setup.sh -t $mod -n $DEV -y -b -u -x
  ./setup.sh -t $mod -z    # cleanup
done


## ** Create production workspaces and release dockers;
##    copy/update WDLs and install methods on Terra
## ** Also populate $prod and $prod_all workspaces with methods
## ** Configure panoply_main panoply_unified_workflow in $prod
cd $panoply/release
./setup-release.sh -p $DEV -N $REL -T $VER -r $proj -w $prod_all -y $prod 


## ** Update notebook docker
cd $panoply/panda
./build-notebook-docker.sh -n $REL -g $VER -u


## ** Commit and push version branch to GitHub (do not merge with dev)
git commit -a -m "Release $VER"
git push origin $rel_id


echo -e "${red}### \n### COPY PANOPLY-startup-notebook to Production Workspaces \n###${reg}"
