#!/bin/bash

cd ..
pgdac=`pwd`
cd hydrant

primary=$pgdac/hydrant/primary-dockerfile
secondary=$pgdac/hydrant/secondary-dockerfile

freshtask() {
 echo "Creating fresh task...";
 ( cd $pgdac/hydrant/tasks/;
   rm -r $task;
   mkdir -p configs;
   echo -e "[Task $task]\nSrc=$pgdac/src/$task" > configs/$task"-config.txt";
   hydrant init $task -c configs/$task"-config.txt" -n 0 )
}

copysrc() {
  echo "Copying source files to docker dir...";
  cp -R $pgdac/src/$task/* $pgdac/hydrant/tasks/$task/$task/src/.;
}

dockertemplate() {
  dockerfile=$1;
  cp $dockerfile $pgdac/hydrant/tasks/$task/$task/Dockerfile;
}

copywdl() {
  wdl=$1
  cp $wdl $pgdac/hydrant/tasks/$task/$task/pgdac_$task.wdl
}

buildpushdocker() {
  cd $pgdac/hydrant/tasks/$task/; 
  if [[ -z "$docker_tag" ]]; then
    hydrant build -n $docker_ns
    hydrant push -n $docker_ns
  else
    docker build -t $docker_ns/$wf_name:$docker_tag .
    docker push $docker_ns/$wf_name:$docker_tag
  fi
}

editdockerfile() {
  ( cd $pgdac/hydrant/tasks/$task/$task;
    sed -i '' "s|broadcptac/r-util:2|broadcptac/$base_task:1|g" Dockerfile;
    sed -i '' "s|broadcptac/pgdac_common:1|broadcptac/$base_task:1|g" Dockerfile; )
}

copycommonutils(){
  echo "Copying common data, packages, r-utilities to docker dir..."
  ( cd $pgdac/hydrant/tasks/$task/$task/;
    mkdir -p data;
    mkdir -p packages;
    mkdir -p R-utilities;
    cp $pgdac/data/* data/.;
    cp $pgdac/hydrant/packages/* packages/.;
    cp $pgdac/hydrant/r-util/* R-utilities/.;
    echo $'!data\n!packages\n!R-utilities' >> .dockerignore )
}

display_usage() {
  echo "usage: ./setup.sh -t [task_name] [-f] [-p] [-s] [-c [custom_docker_file]] "
  echo "                      [-w [wdl]] [-m [base_task]] [-n [docker_namespace]] "
  echo "                      [-g [docker_tag_num]] [-y] [-b] [-h]"
  echo -e "\nNOTE: If using -f, use it after -t and before all other options.\n\n"
  echo "-t | string | task name"
  echo "-f | flag   | Erase task directory if exists; Initialize a new task directory space"
  echo "-p | flag   | Use the primary-dockerfile template which builds on top of r-util docker; "
  echo "            | Copy common utilities to the docker as well"
  echo "-s | flag   | Use the secondary-dockerfile which builds on top of pgdac_common docker"
  echo "-c | string | Use the custom dockerfile with its full path specified in the argument"
  echo "-w | string | Copy the wdl with its full path specified in the argument"
  echo "-m | string | Edit the dockerfile to replace pgdac_common with the task specified in the argument"
  echo "-n | string | Docker namespace for building and pushing"
  echo "-g | string | Docker tag number for building and pushing"
  echo "-y | flag   | Copy source files from the PGDAC src directory to the hydrant's docker src directory"
  echo "-b | flag   | Build and push the docker to dockerhub with the specified docker namespace"
  echo "-h | flag   | Print Usage"
}

while getopts ":t:c:w:n:m:g:psfybh" opt; do
    case $opt in 
        t) task="$OPTARG";;
        p) dockertemplate $primary; copycommonutils;;
        s) dockertemplate $secondary;;
        c) docker_custom="$OPTARG"; dockertemplate $docker_custom;;
        w) wdl="$OPTARG"; copywdl $wdl;; 
        n) docker_ns="$OPTARG";;
        m) base_task="$OPTARG"; editdockerfile $base_task;;
        g) docker_tag="$OPTARG";;
        f) freshtask;;
        y) copysrc;;
        b) buildpushdocker;;
        h) display_usage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done
