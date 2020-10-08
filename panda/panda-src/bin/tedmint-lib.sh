#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

set_globals()
{
  get_bucket $wkspace $project
  inp=$dst/pipeline-input
  blu='\033[1;34m'
  grn='\033[1;32m'
  yel='\033[1;33m'
  red='\033[0;31m'
  reg='\033[0m'   # No Color
  yon="${grn}[ enter yes | no ]"
  css="${yel}[ case-sensitive ]"
  inf="${blu}[ notification-> ]"
  err="${red}error.${reg}"
  not="${grn}----->${reg}" ## notification
}

trim()
{
  echo -e "$( echo -e "$1" | awk '{$1=$1};1' )"
}

ftrim()
{
  echo -e "$( echo -e "$1" | awk '{$1=$1};1' | awk '{print tolower($0)}' )"
}

get_bucket()
{
  bucket=$(fissfc space_info -w $wkspace -p $project |
          perl -nle 'print $& if m{"bucketName":.*?[^\\]"}' |
          awk -F: '{print $2}' |
          awk -F\" '{print $2}')
  echo "$bucket"
  # bucket global
}

