#!/bin/bash

# %prog arg1 arg2 arg3 arg4
# Steps before running this script: 
#  = should have data directory ne directory above this file's directory
#  = should have packages and r-util directories in this file's directory
#  = should have a directory names dependencies in this file's directory
#  = should have a file with the name (arg2).txt w/o newine separated entries
#    for each source file required to run this work-flow

# OPTIONS:  
# arg1 - subdirectory name [has to exist] in which to create workflow dir
# arg2 - workflow name [hydrant will create directory inside subdir]
# arg3 - docker namespace
# arg4 - firecloud namespace


wspace="/Users/rkothadi/Documents/Code/GithubClones/PGDAC"
while getopts ":s:n:i:w:d:f:m:" opt; do
    case $opt in 
        s) sub_dir="$OPTARG"
        ;;
        n) wf_name="$OPTARG"
        ;;
        i) docker_file="$OPTARG"
        ;;
        w) wdl="$OPTARG"
        ;; 
        d) docker_ns="$OPTARG"
        ;;
        f) fire_ns="$OPTARG"
        ;;
        m) modify_docker_from="$OPTARG"
        ;;
        \?) echo "Invalid Option -$OPTARG" >&2
        ;;
    esac
done

if [ "$sub_dir" = "" ] || [ "$wf_name" = "" ]; then
    echo "Usage: %prog -s dir_name -w workflow_name";
    echo "Program exiting..."
   exit
fi
cd $sub_dir
echo -e "[Task $wf_name]\nSrc=$wspace/src/$wf_name" > $wf_name"-config.txt"
hydrant init $wf_name -c $wf_name"-config.txt" -n 0
cd $wf_name/$wf_name
if [[ ! -z "$docker_file" || ! -z $modify_docker_from ]]; then
  cp $wspace/hydrant/$docker_file Dockerfile
  if [[ ! -z $modify_docker_from ]]; then
      sed -i '' "s|broadcptac/pgdac_common:1|broadcptac/$modify_docker_from:1|g" Dockerfile
  fi
  hydrant build -n $docker_ns
  hydrant push -n $docker_ns
fi

cd ..
cd ..

# edit inputs.json
# hydrant test manually
# hydrant install manually
cd - > /dev/null
