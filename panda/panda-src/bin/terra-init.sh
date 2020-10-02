#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

source config.sh
source $src/tedmint-lib.sh

echo -e "\n---------------------------"
echo -e "initializing Terra w-space "
echo -e "---------------------------\n"

fissfc space_new -p $project -w $wkspace
fissfc space_set_acl -p $project -w $wkspace -r OWNER --users $group
bucket=$( get_bucket $wkspace $project )
echo -e "bucket=\"$bucket\"" >> config.sh
