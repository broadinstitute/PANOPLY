#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
source config.sh
source $src/tedmint-lib.sh

get_update_list() 
{
  local flag arr os switch ss
  erase_list=
  add_list=
  osets=
 
  ## If no new Sample sets have been set; exit
  ##  TODO: Figure out exit codes
  echo -e "$sets"
 
  if [[ -z $sets ]]; then
    printf "$inf ${reg}You have not created any new sample sets in this run.\n"
    return
  fi
  
  ## Get the current set-lists from Terra for the wkspace
  osets=`fissfc sset_list -w $wkspace -p $project`
  IFS=$'\n' read -rd '' -a osets <<<"$osets"
  
  printf "$inf ${reg}Sets on Terra:"
  for os in "${osets[@]}"
  do
    printf " $os |"
  done
  
  printf "\n$css ${reg}Refresh sets from terra? [some/all/default:'no'] = "
  read switch; switch=$( trim "$switch" )
  
  ## Add all sample sets in config file to add_list
  if [[ $switch == 'all' ]]; then
    add_list+=( "${osets[@]}" )

  ## Add only some sample sets to add_list
  elif [[ $switch == 'some' ]]; then
    printf "$css ${reg}List terra sample-sets to refresh (;) = "
    read arr; arr=$( trim "$arr")
    IFS=';' read -ra arr <<<"$arr"
    add_list+=( "${arr[@]}" )
  fi
  
  if [[ ${#add_list[@]} -eq 0 ]]; then
    return
  fi
  
  ## Refine the manually entered sample set names in add_list
  add_list=( $( echo "${add_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )

  ## Get sample sets that have to deleted from Terra first in order to be "updated"
  ## Erase them one by one
  erase_list=( $( comm -12 <( printf '%s\n' "${osets[@]}" | LC_ALL=C sort) \
                           <( printf '%s\n' "${add_list[@]}" | LC_ALL=C sort) ))
  for ss in "${erase_list[@]}"
  do
      echo -e "yes\n" | fissfc sset_delete -w $wkspace -p $project -e $ss
      echo ""
  done
  
  ## Send back the sample sets that need to be uploaded in the add_list
  printf "$inf ${reg}Updating:"
  for os in "${add_list[@]}"
  do
    printf " $os |"
  done
  add_list=$( IFS=';'; echo "${add_list[*]}" )
}
