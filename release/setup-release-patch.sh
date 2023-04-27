#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

start_dir=`pwd`    # must be invoked from the PANOPLY/release directory
panoply=$start_dir/..
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}Error:${reg}"
not="${grn}====>>${reg}" ## notification
modules_all=( $( ls -d $panoply/hydrant/tasks/panoply_* | xargs -n 1 basename ) )


### This script is for PATCH FIXING a release
### It executes setup-release.sh with standard dev parameters and the -P patch-fix flag
### NOTE: no analyses are performed WITHIN this script; it JUST a wrapper for setup-release.sh

display_usage() {
  echo "usage: ./setup-release-patch.sh -v <VER> [-h]"
  echo "       -v <VER>      specify version string (required)"
  echo "       <modules>     specify modules to be rebuilt; no flag required"
  echo "       -h  print usage"
  exit
}

while getopts ":v:M:h" opt; do
    case $opt in
        v) VER="$OPTARG";;
        h) display_usage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done
shift "$(( OPTIND - 1 ))" # discard arguments parsed by getopts
modules=( $@ ) # parse remaining arguments as modules

# Validate VER Parameter
if [[ -z $VER ]]; then
  echo -e "$err Release version not specified. Exiting.."
  exit 1
elif ! [[ "$VER" =~ ^[a-zA-Z0-9_\-]+$ ]]; then
  echo -e "$err Invalid version $VER. Must be alphanumeric with dash/unscore"
  exit 1
elif ! [[ -d version-$VER ]]; then
  echo -e "$err Invalid version $VER. Version must already exist."
  exit 1
fi

# Validate modules Parameter
if [[ -z $modules ]]; then
  echo -e "$err No modules specified. Exiting.."
  exit 1
fi
for mod in "${modules[@]}"; do
  if ! [[ ${modules_all[@]} =~ $mod ]]; then
    echo -e "$err Invalid module $mod. Module must match a task from ../hydrant/tasks"
    exit 1
  fi
done


### Run ./setup-release.sh with Patch-Flag, for these parameters

# default development and release namespaces (copied from ./release.sh)
DEV=broadcptacdev
REL=broadcptac
proj=broad-firecloud-cptac
prod="PANOPLY_Production_Pipelines_v$VER"
prod_all="PANOPLY_Production_Modules_v$VER"

./setup-release.sh -p $DEV -N $REL -T $VER -r $proj -w $prod_all -y $prod -P ${modules[@]}

