#!/bin/bash

cd ..
pgdac=`pwd`
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}error.${reg}"
not="${grn}----->${reg}" ## notification

primary=$pgdac/hydrant/primary-dockerfile
secondary=$pgdac/hydrant/secondary-dockerfile

# deletes existing folder if present, 
# creates a config file for hydrant, 
# and creates a new task with hydrant
freshtask() {
 echo -e "$not Creating fresh task...";
 ( cd $pgdac/hydrant/tasks/;
   if [ -d "$task" ]; then rm -r $task; fi
   mkdir -p configs;
   echo -e "[Task $task]" > configs/$task"-config.txt";
   hydrant init $task -c configs/$task"-config.txt" -n 0
   rm -rf $task/$task/src; )
}

# copies source code from ../src/ location to the 
# docker location
copysrc() {
  echo -e "$not Copying source files to docker dir...";
  mkdir -p $pgdac/hydrant/tasks/$task/$task/src/;
  cp -R $pgdac/src/$task/* $pgdac/hydrant/tasks/$task/$task/src/.;
}


# copies appropriate dockerfile to the task
dockertemplate() {
  dockerfile=$1;
  cp $dockerfile $pgdac/hydrant/tasks/$task/$task/Dockerfile;
}

# copies appropriate WDL to the task
copywdl() {
  wdl=$1
  cp $wdl $pgdac/hydrant/tasks/$task/$task/pgdac_$task.wdl
}

# builds a new docker with the default docker tag set to 
# the latest commit hash from github
builddocker() {
  cd $pgdac/hydrant/tasks/$task/$task/;
  if [[ $x_flag != "true" ]]; then
    echo -e "$not Pruning docker images on this system to ensure new build..."
    yes | docker system prune --all;
  fi
  echo -e "$not Building $task locally...";
  echo "!data\n!packages\n!R-utilities" > .dockerignore;
  if [[ $task == "pgdac_utils" ]]; then
    git clone https://github.com/broadinstitute/proteomics-Rutil.git
    mv proteomics-Rutil R-utilities
  fi
  if [[ $task == "pgdac_common" ]]; then
    mkdir -p data
    cp -r $pgdac/data/* data/.
  fi
  docker build --rm --no-cache -t $docker_ns/$task:$docker_tag . ;
  docker images | grep "$task"
}

# login into dockerhub and push the built docker
pushdocker()
{
  docker login
  cd $pgdac/hydrant/tasks/$task/$task/;
  echo -e "$not Pushing $task:$docker_tag to dockerhub...";
  docker push $docker_ns/$task:$docker_tag;
  open https://hub.docker.com/repository/docker/$docker_ns/$task
}

# replace the docker namespace and docker tag in the existing WDL
# to the current docker namespace and docker tag
replaceinwdl(){
  ( cd $pgdac/hydrant/tasks/$task;
    wdl_dns=`grep "/$task:" $task.wdl | cut -d'"' -f2 | cut -d'/' -f 1`
    wdl_tag=`grep "/$task:" $task.wdl | cut -d'"' -f2 | cut -d':' -f 2`
    sed -i '' "s|$wdl_dns/$task:$wdl_tag|$docker_ns/$task:$docker_tag|g" $task.wdl; )
}

editdockerfile() {
  ( cd $pgdac/hydrant/tasks/$task/$task;
    sed -i '' "s|broadcptac/pgdac_common:.*|broadcptac/$base_task_with_tag|g" Dockerfile; )
}

display_usage() {
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
  echo "| -m | string | Replace pgdac_common:<ver> with argument"
  echo "|    |        | Note: 'broadcptac' is the default namespace. ( Uncustomizable for now )"
  echo "| -n | string | Docker namespace"
  echo "| -g | string | Manually overrides docker tag number which is set by default to the "
  echo "|    |        | latest commit hash in the repository."
  echo "| -y | flag   | Add PGDAC/src/task/ files to docker or update them"
  echo "| -b | flag   | Build docker"
  echo "| -u | flag   | Push docker to dockerhub with the specified docker namespace"
  echo "| -x | flag   | Do not prune dockers from local system before building"
  echo "| -z | flag   | Cleanup children directories before Github commit"
  echo "| -h | flag   | Print Usage"
  echo "==============================================="
  exit
}

while getopts ":t:c:w:n:m:g:prsfybuxzh" opt; do
    case $opt in
        t) task="$OPTARG"; wf_name="$task";;
        p) p_flag="true";;
        s) s_flag="true";;
        c) docker_custom="$OPTARG"; dockertemplate $docker_custom;;
        w) wdl="$OPTARG"; copywdl $wdl;;
        n) docker_ns="$OPTARG";;
        m) base_task_with_tag="$OPTARG";;
        g) docker_tag="$OPTARG";;
        f) f_flag="true";;
        y) y_flag="true";;
        b) b_flag="true";;
        u) u_flag="true";;
        r) r_flag="true";;
        x) x_flag="true";; ##if calling from update.sh
        z) z_flag="true";;
        h) display_usage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

rmdata()
{
  cd $pgdac/hydrant/tasks/$task/$task/;
  if [ -d "data" ]; then
    rm -rf data;
  fi
}

cleanup()
{
  cd $pgdac/hydrant;
  find . -name "tests" -type d -exec rm -rf {} \;
  find . -name "src" -type d -exec rm -rf {} \;
  find . -name "hydrant.cfg" -exec rm {} \;
  find . -name "R-utilities" -type d -exec rm {} \;
  rm *.Rout
}


main()
{

  if [[ $z_flag == "true" ]]; then
    cleanup
    exit
  fi

  if [[ -z "$task" ]]; then
    echo -e "$err Task name not entered. Exiting."
    exit
  fi

  ## fresh task
  if [[ $f_flag == "true" ]]; then
    freshtask;
  fi 

  ## dockerfile copy
  if [[ $p_flag == "true" ]]; then
    dockertemplate $primary;
  fi

  if [[ $s_flag == "true" ]]; then
    dockertemplate $secondary;
    rmdata;  
  fi

  if [[ ! -z "$docker_custom" ]]; then
    dockertemplate $docker_custom;
    rmdata;
  fi

  if [[ ! -z "$wdl" ]]; then
    copywdl $wdl;
  fi

  ## edit dockerfile
  if [[ ! -z "$base_task_with_tag" ]]; then
    editdockerfile $base_task_with_tag;
  fi

  ## copy source
  if [[ $y_flag == "true" ]]; then
    copysrc;
  fi 

  ## docker tag
  if [[ -z $docker_tag ]]; then
    docker_tag=`git log -1 --pretty=%h`
  fi

  ## build docker
  if [[ $b_flag == "true" ]]; then
    if [[ ! -z "$docker_ns" ]]; then
      builddocker;
    elif [[ -z "$docker_ns" ]]; then
      echo -e "$err Docker namespace not entered. Exiting."
      exit
    fi
  fi 

  ## push docker
  if [[ $u_flag == "true" ]]; then
    if [[ ! -z "$docker_ns" ]]; then
      pushdocker;
    elif [[ -z "$docker_ns" ]]; then
      echo -e "$err Docker namespace not entered. Exiting."
      exit
    fi
  fi 

  # replace docker reference in the WDL
  if [[ $r_flag == "true" ]]; then
    replaceinwdl;
  fi

}

main
