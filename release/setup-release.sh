#!/bin/bash
start_dir=`pwd`    # must be invoked from the PANOPLY/release directory
panoply=$release_dir/..
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}error.${reg}"
not="${grn}----->${reg}" ## notification

base_url="https://registry.hub.docker.com/v2/repositories/"

display_usage() {
  echo "usage: ./release.sh -p pull_dns -T release_tag "
  echo "                    [-N [release_dns]] [-h]"
  echo "-N | string | Docker namespace for building and pushing release"
  echo "-T | string | Docker release tag"
  echo "-p | string | Docker namespace to pull latest dev images from"
  echo "-h | flag   | Print Usage"
  exit
}

# replace the docker namespace and docker tag in the existing WDL
# to the specified docker namespace and docker tag
replaceDockerInWdl() {
  task=$1
  new_ns=$2
  new_tag=$3
  wdl_dns=`grep "/$task:" $task.wdl | cut -d'"' -f2 | cut -d'/' -f 1`
  wdl_tag=`grep "/$task:" $task.wdl | cut -d'"' -f2 | cut -d':' -f 2`
  sed -i '' "s|$wdl_dns/$task:$wdl_tag|$new_ns/$task:$new_tag|g" $task.wdl
}


while getopts ":p:T:N:h" opt; do
    case $opt in
        p) pull_dns="$OPTARG";;
        T) release_tag="$OPTARG";;
        N) release_dns="$OPTARG";;
        h) display_usage;;
        \?) echo "Invalid Option -$OPTARG" >&2;;
    esac
done

if [[ -z $pull_dns ]]; then
  echo -e "$err Development Docker Namespace required. Exiting.."
  exit 1
fi

if [[ -z $release_tag ]]; then
  echo -e "$err Release Tag required. Exiting.."
  exit 1
fi

if [[ -z $release_dns ]]; then
  release_dns=$pull_dns
fi

modules=( $( ls -d $panoply/tasks/panoply_* | xargs -n 1 basename ) )
release_dir=version-$release_tag

mkdir -p $release_dir;
yes | docker system prune --all

for mod in "${modules[@]}"
do
  mkdir -p $release_dir/$mod
  url=$base_url$pull_dns/$mod/tags
  lat=( $( curl -s -S "$url" | \
             jq '."results"[]["name"]' | \
             sed -n 1p | cut -d'"' -f2 ) )
  if [[ -z $lat ]]; then
    continue;
  fi
  
  mkdir -p $release_dir/$mod;
  cd $release_dir/$mod;
  
  # build release docker
  echo -e "FROM $pull_dns/$mod:$lat" > Dockerfile
  docker build --rm --no-cache -t $release_dns/$mod:$release_tag . ;
  docker images | grep "$mod"
  docker login
  docker push $release_dns/$mod:$release_tag
  
  # copy and update task WDL
  wdl_dir=$panoply/tasks/panoply_$mod/panoply_$mod
  wdl=panoply_$mod.wdl
  if [[ -f $wdl_dir/$wdl ]]; then
    cp $wdl_dir/$wdl $wdl
    replaceDockerInWdl panoply_$mod $release_dns $release_tag
    # install method with synopsis/documentation
    # get method snapshot id using fissfc meth_list and save release snapshot ids
  fi
  cd $start_dir;
done
