#!/bin/bash

cd ..
pgdac=`pwd`

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
  mkdir -p $pgdac/hydrant/tasks/$task/$task/src/;
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
  cd $pgdac/hydrant/tasks/$task/$task/;
  if [[ -z "$docker_tag" ]]; then
    docker_tag=1
  fi
  echo "Building $task locally...";
  echo "!data\n!packages\n!R-utilities" > .dockerignore;
  docker build --rm --no-cache -t $docker_ns/$task:$docker_tag . > buildlog.txt 2>&1;
  echo "Pushing $task to dockerhub...";
  docker push $docker_ns/$task:$docker_tag >> buildlog.txt 2>&1;
  rm -rf R-utilities;
  rm -rf data;
  rm -rf packages;
  rm -rf src;
  rm .dockerignore;
  rm buildlog.txt;
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
    mkdir -p src
    cp $pgdac/data/* data/.;
    echo $'!data' >> .dockerignore )
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

print_log() {
  echo "Cleaning up..."
  mkdir -p $pgdac/hydrant/logs
  stamp=$(date "+%S%M%H%d%m%Y")
  logfile="$pgdac/hydrant/logs/$task-log-$stamp.txt"
  echo "Log-File: Parameter List" > $logfile;
  echo "task=$task" >> $logfile;
  echo "docker_custom=$docker_custom" >> $logfile;
  echo "wdl=$wdl" >> $logfile;
  echo "docker_ns=$docker_ns" >> $logfile;
  echo "base_task=$base_task" >> $logfile;
  echo "docker_tag=$docker_tag" >> $logfile;
}
print_log
