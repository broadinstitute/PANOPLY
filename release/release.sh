#!/bin/bash

## Release a new version of PANOPLY

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

while getopts ":vh" opt; do
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

## ** Get latest code version in dev branch, create release branch
cd $panoply
git checkout dev
git pull
# create new branch for release
rel_id="release-$VER"
git checkout -b $rel_id dev
git push origin $rel_id


## ** Rebuild all docker images in $DEV
# Build / update panoply_libs
cd $panoply/hydrant
./setup.sh -t panoply_libs -n $DEV -y -b -u 

# Build panoply_utils and panoply_common
./setup.sh -t panoply_utils -n $DEV -y -b -u 
./setup.sh -t panoply_common -n $DEV -y -b -u -x

# Rebuild all other task dockers
./update.sh -t panoply_common -n $DEV


## ** Create release dockers, copy/update WDLs and install methods on Terra
cd $panoply/release
./setup-release.sh -p $DEV -N $REL -T $VER


## Step 5: Update and upload WORKFLOW WDLs (with WDL imports) to Terra


## Step 6: Create production workspace for release


## Step 7: Copy all methods/data to new workspace and configure
