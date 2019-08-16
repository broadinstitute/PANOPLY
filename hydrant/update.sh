#!/bin/bash

cd ..
pgdac=`pwd`
cd hydrant

display_usage() {
  echo "usage: ./update.sh -t [task_name] -n [docker_namespace] [-h]"
  echo "-t | string | task name"
  echo "-n | string | Docker namespace for building and pushing"
  echo "-h | flag   | Print Usage"
  exit
}

while getopts ":t:n:h" opt; do
    case $opt in 
        t) task="$OPTARG";;
        n) docker_ns="$OPTARG";;
        h) display_usage;;
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
  ./setup.sh -t $task -n $docker_ns -y -b
done
