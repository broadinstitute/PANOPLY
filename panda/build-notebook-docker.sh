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
  echo "| -u | flag   | Push docker to dockerhub/gcr.io with the specified namespace"
  echo "| -h | flag   | Print Usage"
  echo "==============================================="
  exit
}


while getopts "n:g:auh" opt; do
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
  # copy panoply-Rutil repository
  git clone https://github.com/broadinstitute/proteomics-Rutil.git
  mv proteomics-Rutil R-utilities
  base_docker1="broadcptacdev/panda_config_libs:$docker_tag" # panda Dockerfile references broadcptacdev namespace explicitly
  base_docker2="broadcptacdev/panda_config_libs:latest"  # panda Dockerfile references broadcptacdev namespace explicitly
  docker build --rm --no-cache -t $base_docker1 -t $base_docker2 .
  rm -rf R-utilities # cleanup
fi

# final docker
cp $panoply/src/panoply_common/master-parameters.yaml $panoply/panda/panda-src/defaults/. # copy in code
cd $panoply/panda/panda-src
final_docker1="$docker_ns/panda:$docker_tag"
final_docker2="$docker_ns/panda:latest"
docker build --rm --no-cache -t $final_docker1 -t $final_docker2 .


## push dockers if requested
if [[ $u_flag == "true" ]]; then
  docker login
  echo -e "Pushing images to dockerhub...";
  if [[ $a_flag == "true" ]]; then
    docker push broadcptacdev/panda_config_libs:$docker_tag # panda Dockerfile references broadcptacdev namespace explicitly
    docker push broadcptacdev/panda_config_libs:latest # panda Dockerfile references broadcptacdev namespace explicitly
  fi
  docker push $docker_ns/panda:$docker_tag
  docker push $docker_ns/panda:latest
  
  # also push to gcr.io
  docker tag $docker_ns/panda:$docker_tag gcr.io/broadcptac/panda:$docker_tag # push to broadcptac, there is no broadcptacdev on gcr.io
  docker push gcr.io/broadcptac/panda:$docker_tag # push to broadcptac, there is no broadcptacdev on gcr.io
fi
