#! /bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#


default_ns="broadcptacdev"

cd ..
panoply=`pwd`


# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color


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
  echo "| -d | flag   | Pull current version of Hallmark geneset and PTM-signature databases."
  echo "| -h | flag   | Print Usage"
  echo "==============================================="
  exit
}


while getopts "n:g:auhd" opt; do
  case $opt in
    n) docker_ns="$OPTARG";;
    g) docker_tag="$OPTARG";;
    a) a_flag="true";;
    u) u_flag="true";;
    d) d_flag="true";;
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


## update ssGSEA and PTM-SEA databases
if [[ $d_flag == "true" ]]; then
  cd $panoply/panda

  ## clone ssGSEA repository
  git clone https://github.com/broadinstitute/ssGSEA2.0.git

  ## update Hallmarks Pathway
  hallmark=`ls $panoply/panda/ssGSEA2.0/db/msigdb/h.all.v*`
  if [[ -n $hallmark ]]; then
    hallmark_old=`ls $panoply/panda/panda-src/defaults/h.all.v*`
    if [[ -n $hallmark_old ]]; then
      rm $hallmark_old
    fi
    echo $GREEN Updating Hallmarks Geneset DB to \'`basename $hallmark`\' $NC
    cp $hallmark $panoply/panda/panda-src/defaults/
  else
    echo $RED Could not find new hallmark pathway database in ssGSEA2.0 repository $NC
  fi

  ## update PTM-Signature Database
  ptmsig_ver=`ls -v $panoply/panda/ssGSEA2.0/db/ptmsigdb/ | tail -n 1`
  if [[ -n $ptmsig_ver ]]; then
    ptmsig_fl=`ls $panoply/panda/ssGSEA2.0/db/ptmsigdb/$ptmsig_ver/all/ptm.sig.db.all.flanking.human.$ptmsig_ver.gmt`
    ptmsig_uni=`ls $panoply/panda/ssGSEA2.0/db/ptmsigdb/$ptmsig_ver/all/ptm.sig.db.all.uniprot.human.$ptmsig_ver.gmt`
    if [[ -n $ptmsig_fl && -n $ptmsig_uni ]]; then
      ptmsig_old=`ls $panoply/panda/panda-src/defaults/ptm.sig.db.all*`
      if [[ -n $ptmsig_old ]]; then
        rm $ptmsig_old
      fi
      echo $GREEN Updating PTM-Signature DBs to \'`basename $ptmsig_fl`\' and \'`basename $ptmsig_uni`\' $NC
      cp $ptmsig_fl $panoply/panda/panda-src/defaults/
      cp $ptmsig_uni $panoply/panda/panda-src/defaults/
    else
      echo $RED Could not find new PTM-Signature database\(s\) in ssGSEA2.0 repository $NC
    fi
  else
    echo $RED Could not find PTM-Signature version directory in ssGSEA2.0 repository $NC
  fi

  ## cleanup
  rm -rf $panoply/panda/ssGSEA2.0
fi



# final docker
cp $panoply/src/panoply_common/master-parameters.yaml $panoply/panda/panda-src/defaults/. # copy in parameters
cp $panoply/src/panoply_metaboanalyst/pathway_db/compound_db.qs $panoply/panda/panda-src/defaults/. # copy in metabolite ID mapping
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
