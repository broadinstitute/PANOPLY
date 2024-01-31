#!/bin/bash
#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

## Release a new version of the PTM-SEA Notebook
##  (must be invoked from the ptm-sea-notebook directory)
## Original Release Documentation can be found at <https://docs.google.com/document/d/1Y7ik_iYhx6Df7Jy9aqOZua-ZGNRnUU4iTnPgo_a8BDQ/>

red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}Error:${reg}"
warn="${grn}Warning:${reg}" ## notification

# ## users (for method permissions)
# proteomics_comp=(manidr@broadinstitute.org nclark@broadinstitute.org wcorinne@broadinstitute.org)
# proteomics_cptac=(GROUP_Broad_CPTAC@firecloud.org)


display_usage() {
  echo "usage: release.sh -v <VER> [-f] [-h]"
  echo "       -v <VER>  specify docker version"
  echo "       -C  clear cache and build docker from scratch"
  echo "       -h  print usage"
  exit
}

while getopts "v:Ch" opt; do
    case $opt in
        v) VER="$OPTARG";;
        C) C=true;; # clear cache
        h) display_usage;;
        \?) echo "$warn Invalid Option -$OPTARG" >&2;;
    esac
done

if [[ -z $VER ]]; then
  echo -e "$err Release version not specified. Exiting.."
  exit 1
fi

# ensure that VER has only numeric and dot
if ! [[ "$VER" =~ ^[0-9.]+$ ]]; then
  echo -e "$err Invalid version. Must be numeric, separated by dots (e.g. 0.5.1)"
  exit 1
fi


######

# should never change
project=broad-firecloud-cptac
ws="PTM-SEA_Notebook"


### DOCKER HUB
echo -e "$not Logging in to Docker Hub"
docker login

# build, tag, and push docker
mkdir src # make temporary src/ directory
cp ../src/panoply_ssgsea/panoply_ptmsea_functions.R src/panoply_ptmsea_functions.R # copy files into /src 
cp ../src/panoply_ssgsea_report/rmd-ssgsea-functions.R src/rmd-ssgsea-functions.R # copy files into /src  
cp ../src/panoply_ssgsea_report/rmd-ssgsea.r src/rmd-ssgsea.r # copy files into /src 
docker build . ${C:+--no-cache} -t "broadcptac/ptm-sea:$VER"
rm -r src # remove temporary src/ directory
docker tag broadcptac/ptm-sea:$VER gcr.io/broadcptac/ptm-sea:$VER # tag as gcr.io/
docker push gcr.io/broadcptac/ptm-sea:$VER # push docker


### WORKSPACE

# ## Update Workspace Permissions (likely not necessary, but doesn't hurt I suppose)
# fissfc space_set_acl -w $ws -p $project -r OWNER --users ${proteomics_comp[@]}
# # setting workspace to PUBLIC is currently not functional -- need to contact Terra support
# # fissfc space_set_acl -w $ws -p $project -r READER --users public
# fissfc space_set_acl -w $ws -p $project -r READER --users ${proteomics_cptac[@]}


## Update Workspace Description
# old_descr=`fissfc space_info -w $ws -p $project | jq -r '.workspace.attributes.description'`
sed "s/\([sS][eE][aA][:_]\)[0-9.]\{1,\}/\1$VER/g" README.md > temp.txt && mv temp.txt README.md # update description in .txt file
# use the python version of the firecloud API to upload a markdown file description
python << EOF
import firecloud.api as fapi
with open('README.md') as f:
  fapi.update_workspace_attributes(
      namespace="$project",
      workspace="$ws",
      attrs=[fapi._attr_set(attr='description', value=f.read())]).json()
EOF


### NOTEBOOK

## Update Notebook with new Docker-Verion
sed "s/ptm-sea:\([0-9.]\{1,\}\)/ptm-sea:$VER/g" PTM_SEA.ipynb > temp.txt && mv temp.txt PTM_SEA.ipynb

## Upload new 
bucket=`fissfc space_info -w $ws -p $project | jq -r '.workspace.bucketName'` # get bucket ID
gsutil cp PTM_SEA.ipynb gs://$bucket/notebooks/PTM_SEA_$VER.ipynb # copy notebook into bucket

