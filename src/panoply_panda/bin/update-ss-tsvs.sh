#!/bin/bash

echo -e "\n---------------------------"
echo -e "updating sample-set models --"
echo -e "---------------------------\n"

source config.sh
source $src/tedmint-lib.sh

online_sets=

read -r -d '' MSG_CONFIG_SS_404 <<- EOM
No sample sets in config.yaml..
Nothing to update..
Exiting..
EOM

read -r -d '' MSG_TERRA_SS <<- EOM
Sample sets found on Terra Workspace ($wkspace) ::
EOM

read -r -d '' MSG_OVERLAP_SS << EOM
The sample sets you in your config file exist on Terra.
Do you want to continue to delete all sample sets on Terra?
You can type "no" and rereun this with new sample set names.
(yes/default:'no') ::
EOM

read -r -d '' MSG_DELETE_SS << EOM
Delete all sets from terra? (yes/default:'no') ::
EOM

## get the names of the sets to be updated
set_names=
set_names_str="set_names=\""
IFS='%' read -ra ssets <<< "$sets"
for sset in "${ssets[@]}"
do
  IFS=',' read -ra temp <<< "$sset"
  set_names_str=$set_names_str"${temp[0]}"
  set_names_str=$set_names_str";"
  set_names+=( "${temp[0]}" )
done
set_names_str=$set_names"\""
echo -e $set_names >> config.sh

if [[ -z $sets ]]; then
  printf "$MSG_SS_404"
  exit
fi

online_sets=`fissfc sset_list -w $wkspace -p $project`
if [[ -z $online_sets ]]; then
  printf "$MSG_TERRA_SS (NULL)\n"
fi

if [[ ! -z $online_sets ]]; then
  IFS=$'\n' read -rd '' -a online_sets <<<  "$online_sets"
  printf "$MSG_TERRA_SS\n"
  for os in "${online_sets[@]}"
  do
    printf "\t* $os\n"
  done

  overlap=( $( comm -12 \
    <( printf '%s\n' "${online_sets[@]}" | LC_ALL=C sort) \
    <( printf '%s\n' "${set_names[@]}" | LC_ALL=C sort) ) )

  if [[ -z $overlap ]]; then
    printf "$MSG_OVERLAP"
    read switch; switch=$( trim "$switch" )
    if [[ $switch != "yes" ]]; then
      exit
    fi
  fi
  printf "$MSG_DELETE_SS\n"
  read switch; switch=$( trim "$switch" )
  if [[ $switch == 'yes' ]]; then
    erase_list+=( "${online_sets[@]}" )
  fi

  for ss in "${erase_list[@]}"
  do
    echo -e "yes\n" | fissfc sset_delete -w $wkspace -p $project -e $ss
    echo ""
  done
fi

Rscript --verbose $src/r-source/populate_sample_set.r \
  -g $gct_types \
  -c $csv_types \
  -s $dst/sample-set-members

mkdir -p sample-set-members-store
mv sample-set-members/* sample-set-members-store/.
fissfc entity_import -w $wkspace -p $project -f sample_set_membership.tsv
fissfc entity_import -w $wkspace -p $project -f sample_set_entity.tsv
