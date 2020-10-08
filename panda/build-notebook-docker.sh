#! /bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


default_ns="broadcptacdev"

cd ..
panoply=`pwd`


displayUsage() {
  echo ""
  echo "usage: ./build-notebook-docker.sh "
  echo "                     [-n [docker_namespace]] "
  echo "                     [-g [docker_tag_num]] "
  echo "                     [-a] [-u] [-h] "
  echo ""
  echo "==============================================="
  echo "| -n | string | Docker namespace (defaults to broadcptacdev)"
  echo "| -g | string | Overrides docker tag (defaults to github commit hash)"
  echo "|    |        | :latest tag always included"
  echo "| -a | flag   | Build both panda_config_libs and panda (default: only panda)"
  echo "| -u | flag   | Push docker to dockerhub with the specified docker namespace"
  echo "| -h | flag   | Print Usage"
  echo "==============================================="
  exit
}


while getopts ":n:gauh" opt; do
  case $opt in
    n) docker_ns="$OPTARG";;
    g) docker_tag="$OPTARG";;
    a) a_flag="true";;
    u) u_flag="true";;
    h) displayUsage;;
    \?) echo "Invalid Option -$OPTARG" >&2;;
  esac
done

## docker namespace
if [[ -z $docker_ns ]]; then
  docker_ns=$default_ns
fi

## docker tag
if [[ -z $docker_tag ]]; then
  docker_tag=`git log -1 --pretty=%h`
fi


## build PANDA base docker and final Terra startup notebook docker
# base docker
if [[ $a_flag == "true" ]]; then
  cd $panoply/panda
  base_docker1="$docker_ns/panda_config_libs:$docker_tag"
  base_docker2="$docker_ns/panda_config_libs:latest"
  docker build --rm --no-cache -t $base_docker1 -t $base_docker2 .
fi

# final docker
cp $panoply/src/panoply_common/master-parameters.yaml $panoply/panda/panda-src/defaults/.
cd $panoply/panda/panda-src
final_docker1="$docker_ns/panda:$docker_tag"
final_docker2="$docker_ns/panda:latest"
docker build --rm --no-cache -t $final_docker1 -t $final_docker2 .


## push dockers if requested
if [[ $u_flag == "true" ]]; then
  docker login
  echo -e "Pushing images to dockerhub...";
  docker push $docker_ns/panda_config_libs:$docker_tag
  docker push $docker_ns/panda_config_libs:latest
  docker push $docker_ns/panda:$docker_tag
  docker push $docker_ns/panda:latest
fi
