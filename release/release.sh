#!/bin/bash

## Release a new version of PANOPLY
##  (must be invoked from the release directory)

panoply=`pwd`/..
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}error.${reg}"
not="${grn}----->${reg}" ## notification

display_usage() {
  echo "usage: release.sh -v release-version"
  exit
}

while getopts "v:h" opt; do
    case $opt in
        v) VER="$OPTARG";;
        h) display_usage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

if [[ -z $VER ]]; then
  echo -e "$err Release version not specified. Exiting.."
  exit 1
fi

# default development and release namespaces
DEV=broadcptacdev
REL=broadcptac
proj=broad-firecloud-cptac
prod="PANOPLY_Production_Pipelines_v$VER"
prod_all="PANOPLY_Production_v$VER"

# Before beginning, check production workspaces do not exist
if ! fissfc space_exists -w $prod_all -p $proj -q; then
  echo "Workspace $prod_all exists. Delete workspace and re-run."
  exit 1
fi
if ! fissfc space_exists -w $prod -p $proj -q; then
  echo "Workspace $prod_all exists. Delete workspace and re-run."
  exit 1
fi


## ** Get latest code version in dev branch, create release branch
# cd $panoply
# git checkout dev
# git pull
# # create new branch for release
# rel_id="release-$VER"
# git checkout -b $rel_id dev
# git push origin $rel_id


## ** Rebuild/update all docker images in $DEV
# Build / update panoply_libs
cd $panoply/hydrant
# ./setup.sh -t panoply_libs -n $DEV -y -b -u 

# Build panoply_utils and panoply_common
# ./setup.sh -t panoply_utils -n $DEV -y -b -u 
./setup.sh -t panoply_common -n $DEV -y -b -u -x

# Rebuild all other task dockers
modules=( $( ls -d $panoply/hydrant/tasks/panoply_* | xargs -n 1 basename |  \
  sed '/panoply_libs/d' | sed '/panoply_utils/d' | sed '/panoply_common/d') )
for mod in "${modules[@]}"
do
  ./setup.sh -t $mod -n $DEV -y -b -u -x
done


## ** Create production workspaces and release dockers;
##    copy/update WDLs and install methods on Terra
## ** Also populate $prod and $prod_all workspaces with methods
cd $panoply/release
./setup-release.sh -p $DEV -N $REL -T $VER -r $proj -w $prod_all -y $prod 


## ** Update notebook docker
cd $panoply/panda
./build-notebook-docker.sh -n $REL -g $VER -u


## ** Configure pipeline/workflows listed below in $prod
wkflows=( panoply_main panoply_unified_pipeline )
for wkflow in "${wkflows[@]}"
do
  #----configure workflows
done
