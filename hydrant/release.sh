#!/bin/bash
panoply=`pwd`
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}error.${reg}"
not="${grn}----->${reg}" ## notification

display_usage() {
  echo "usage: ./release.sh -p pull_dns -T release_tag "
  echo "                    [-N [release_dns]] [-h]"
  echo "-N | string | Docker namespace for building and pushing release"
  echo "-T | string | Docker release tag"
  echo "-p | string | Docker namespace to pull latest dev images from"
  echo "-h | flag   | Print Usage"
  exit
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

base_url="https://registry.hub.docker.com/v2/repositories/"
modules=( $( ls -d $panoply/tasks/panoply_* | xargs -n 1 basename ) )

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

mkdir -p releases;
yes | docker system prune --all

for mod in "${modules[@]}"
do
  mkdir -p releases/$mod
  url=$base_url$pull_dns/$mod/tags
  lat=( $( curl -s -S "$url" | \
             jq '."results"[]["name"]' | \
             sed -n 1p | cut -d'"' -f2 ) )
  if [[ -z $lat ]]; then
    continue;
  fi
  
  mkdir -p releases/$mod;
  cd releases/$mod;
  echo -e "FROM $pull_dns/$mod:$lat" > Dockerfile
  docker build --rm --no-cache -t $release_dns/$mod:$release_tag . ;
  docker images | grep "$mod"
  docker login
  docker push $release_dns/$mod:$release_tag
  cd $panoply;
done
