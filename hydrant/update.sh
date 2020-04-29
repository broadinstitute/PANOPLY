#!/bin/bash

cd ..
panoply=`pwd`
cd hydrant
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}error.${reg}"
not="${grn}----->${reg}" ## notification

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

if [[ -z "$task" ]]; then
  echo -e "$err Task not entered. Exiting."
  exit
fi
if [[ -z "$docker_ns" ]]; then
  echo -e "$err Docker namespace not entered. Exiting."
  exit
fi
if [[ -z "$docker_tag" ]]; then
  docker_tag=`git log -1 --pretty=%h`
fi


## create maps 
R CMD BATCH --vanilla "--args -p $panoply -t $task" map_dependency.r

<<"comment"
## prune local docker images
echo -e "$not Pruning local docker images to ensure new build..."
yes | docker system prune --all;
comment

## <<"comment"
## read targets for current task
ftarget=$panoply/hydrant/tasks/targets/$task-targets.txt
targets=`head -n 1 $ftarget`
IFS=';' read -ra targets <<< "$targets"
## comment

<<"comment"
## Build and push current parent task
./setup.sh -t $task -n $docker_ns -y -b -g $docker_tag -x -u
sleep 60
comment

## <<"comment"
## build and push targets
base_url="https://registry.hub.docker.com/v2/repositories/"
for target in "${targets[@]}"
do
  dockerfile="tasks/$target/$target/dockerfile"
  str=( $( grep "FROM" $dockerfile | cut -d' ' -f2 ) )
  dns=( $( echo -e $str | cut -d'/' -f1 ) )
  par=( $( echo -e $str | cut -d'/' -f2 | cut -d':' -f1 ) )
  tag=( $( echo -e $str | cut -d'/' -f2 | cut -d':' -f2 ) )
  url=$base_url$dns/$par/tags
  latest_tag=( $( curl -s -S "$url" | jq '."results"[]["name"]' | \
                    sed -n 1p | cut -d'"' -f2 ) )
  echo -e "dns:$dns,par:$par,tag:$tag,lat:$latest_tag,task:$target"
  #sed -i '' "s|$tag|$latest_tag|g" $dockerfile;
  ## ./setup.sh -t $target -n $docker_ns -y -b -g $docker_tag -x -u
done
## comment
