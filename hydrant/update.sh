#!/bin/bash

cd ..
pgdac=`pwd`
cd hydrant

display_usage() {
  echo "usage: ./update.sh -t [task_name] -n [docker_namespace] [-h]"
  echo "-t | string | task name"
  echo "-n | string | Docker namespace for building and pushing"
  echo "-g | string | Docker tag"
  echo "-h | flag   | Print Usage"
  exit
}

while getopts ":t:n:g:h" opt; do
    case $opt in 
        t) task="$OPTARG";;
        n) docker_ns="$OPTARG";;
        h) display_usage;;
        g) docker_tag="$OPTARG";;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

if [[ -z "$task" ]] || [[ -z "$docker_ns" ]]; then
  display_usage
  exit
fi

R CMD BATCH --vanilla "--args -p $pgdac -t $task" map_dependency.r
ftarget=$pgdac/hydrant/tasks/targets/$task-targets.txt
targets=`head -n 1 $ftarget`
IFS=';' read -ra tasks <<< "$targets"
for task in "${tasks[@]}"
do
  if [[ $task = r_util ]]; then
    #echo "here $docker_tag"
    ./setup.sh -t $task -n $docker_ns -g $docker_tag -b
  elif [[ $task = pgdac_common ]]; then
    ./setup.sh -t $task -p -n $docker_ns -y -b
  else
    ./setup.sh -t $task -n $docker_ns -y -b
  fi
done
