#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

start_dir=`pwd`    # must be invoked from the PANOPLY/release directory
panoply=$start_dir/..

modules_ws=PANOPLY_Modules_DEV
pipelines_ws=PANOPLY_Pipelines_DEV
project=broad-firecloud-cptac
yaml=master-parameters.yaml

./setup-release.sh \
  -p broadcptacdev \
  -T DEV \
  -N broadcptacdev \
  -w $modules_ws \
  -y $pipelines_ws \
  -r $project

for ws in $modules_ws $pipelines_ws
do
  bucket=`fissfc space_info -w $ws -p $project | jq -r '.workspace.bucketName'`
  gsutil cp $panoply/src/panoply_common/$yaml gs://$bucket/$yaml
done
