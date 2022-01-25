#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#

cd ..
panoply=`pwd`
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}error.${reg}"
not="${grn}----->${reg}" ## notification

primary=$panoply/hydrant/primary-dockerfile
secondary=$panoply/hydrant/secondary-dockerfile

# copies source code from ../src/ location to the 
# docker location
copySrc() {
  echo -e "$not Copying source files to docker dir...";
  mkdir -p $panoply/"third-party-modules"/$task/$task/src/;
  cp -R $panoply/src/$task/* $panoply/"third-party-modules"/$task/$task/src/.;
}

# copies appropriate WDL to the task
copyWdl() {
  wdl=$1
  cp $wdl $panoply/"third-party-modules"/$task/$task/panoply_$task.wdl
}

# builds a new docker with the default docker tag set to 
# the latest commit hash from github
buildDocker() {
  cd $panoply/"third-party-modules"/$task/$task/;
  if [[ $x_flag != "true" ]]; then
    echo -e "$not Pruning docker images on this system to ensure new build..."
    yes | docker system prune --all;
  fi
  echo -e "$not Building $task locally...";
  echo -e "!data\n!packages\n!R-utilities\n!src" > .dockerignore;

  dname_1="$docker_ns/$task:$docker_tag"
  dname_2="$docker_ns/$task:latest"
  docker build --rm --no-cache -t $dname_1 -t $dname_2 . ;
  docker images | grep "$task"
}

displayUsage() {
  echo ""
  echo "usage: ./setup.sh -t [task_name] " 
  echo "                     [-f] [-p] [-s] "
  echo "                     [-c [custom_docker_file]] "
  echo "                     [-w [wdl]] "
  echo "                     [-m [base_task]:[tag]] "
  echo "                     [-n [docker_namespace]] "
  echo "                     [-g [docker_tag_num]] "
  echo "                     [-y] [-b] [-h] [-x] [-z]"
  echo ""
  echo "==============================================="
  echo "| -t | string | Task name"
  echo "| -f | flag   | Erase task directory if exists; Initialize a new task directory space"
  echo "| -p | flag   | Use the primary-dockerfile template"
  echo "| -s | flag   | Use the secondary-dockerfile template"
  echo "| -c | string | Use the custom dockerfile with its full path specified in the argument"
  echo "| -w | string | Copy the wdl with its full path specified in the argument"
  echo "| -m | string | Replace panoply_common:<ver> with argument"
  echo "|    |        | Note: 'broadcptac' is the default namespace. ( Uncustomizable for now )"
  echo "| -n | string | Docker namespace"
  echo "| -g | string | Manually overrides docker tag number which is set by default to the "
  echo "|    |        | latest commit hash in the repository."
  echo "| -y | flag   | Add PGDAC/src/task/ files to docker or update them"
  echo "| -b | flag   | Build docker"
  echo "| -r | flag   | Replace Docker Namespace and Tag in WDL"
  echo "| -a | flag   | Update WDL on terra"
  echo "| -d | string | WDL Snapshot Comment which is set by default to the docker tag ID"
  echo "| -u | flag   | Push docker to dockerhub with the specified docker namespace"
  echo "| -x | flag   | Do not prune dockers from local system before building"
  echo "| -z | flag   | Cleanup task directory"
  echo "| -e | falg   | Cleanup hydrant and 'tasks' directory"
  echo "| -h | flag   | Print Usage"
  echo "==============================================="
  exit
}

while getopts ":t:c:w:n:m:g:d:praesfybuxzh" opt; do
    case $opt in
        t) task="$OPTARG"; wf_name="$task";;
        p) p_flag="true";;
        s) s_flag="true";;
        c) docker_custom="$OPTARG"; dockerTemplate $docker_custom;;
        w) wdl="$OPTARG"; copyWdl $wdl;;
        n) docker_ns="$OPTARG";;
        m) base_task_with_tag="$OPTARG";;
        g) docker_tag="$OPTARG";;
        d) wdl_snapshot_comment="$OPTARG";;
        f) f_flag="true";;
        y) y_flag="true";;
        b) b_flag="true";;
        u) u_flag="true";;
        e) e_flag="true";;
        r) r_flag="true";;
        a) a_flag="true";;
        x) x_flag="true";; ##if calling from update.sh
        z) z_flag="true";;
        h) displayUsage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

main()
{

  if [[ $z_flag == "true" ]]; then
    cleanup
    exit
  fi

  if [[ $e_flag == "true" ]]; then
    cleanup_hydrant
    exit
  fi

  if [[ -z "$task" ]]; then
    echo -e "$err Task name not entered. Exiting."
    exit
  fi

  ## fresh task
  if [[ $f_flag == "true" ]]; then
    freshTask;
  fi 

  ## dockerfile copy
  if [[ $p_flag == "true" ]]; then
    dockerTemplate $primary;
  fi

  if [[ $s_flag == "true" ]]; then
    dockerTemplate $secondary;
    rmData;  
  fi

  if [[ ! -z "$docker_custom" ]]; then
    dockerTemplate $docker_custom;
    rmData;
  fi

  if [[ ! -z "$wdl" ]]; then
    copyWdl $wdl;
  fi

  ## edit dockerfile
  if [[ ! -z "$base_task_with_tag" ]]; then
    editDockerfile $base_task_with_tag;
  fi

  ## copy source
  if [[ $y_flag == "true" ]]; then
    copySrc;
  fi 

  ## docker tag
  if [[ -z $docker_tag ]]; then
    docker_tag=`git log -1 --pretty=%h`
  fi

  ## build docker
  if [[ $b_flag == "true" ]]; then
    if [[ ! -z "$docker_ns" ]]; then
      buildDocker;
    elif [[ -z "$docker_ns" ]]; then
      echo -e "$err Docker namespace not entered. Exiting."
      exit
    fi
  fi 

  ## push docker
  if [[ $u_flag == "true" ]]; then
    if [[ ! -z "$docker_ns" ]]; then
      pushDocker;
    elif [[ -z "$docker_ns" ]]; then
      echo -e "$err Docker namespace not entered. Exiting."
      exit
    fi
  fi 

  ## replace docker reference in the WDL
  if [[ $r_flag == "true" ]]; then
    replaceDockerInWdl;
  fi

  ## update the WDL on Terra
  if [[ $a_flag == "true" ]]; then
    updateWdlOnTerra;
  fi

}

main