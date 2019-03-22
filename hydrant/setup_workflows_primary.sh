#!/bin/bash

# Just build Docker


display_usage() {
    echo "usage: %prog [-s <sub_dir>] [-n <wf_name>] [[optional]]"
    echo "-s | <sub_dir>            | string | Directory inside PGDAC/hydrant/"
    echo "-n | <wf_name>            | string | Workflow name that you wish to keep"
    echo "-i | <docker_file_temp>   | string | Name of docker template inside PGDAC/hydrant/"
    echo "-o | <docker_file_notemp> | path   | Full path of the dockerfile if one doesn't want to use the template"
    echo "-w | <wdl>                | path   | If you wish to overwrite the WDL template inside task directory"
    echo "-d | <docker_ns>          | string | Docker namespace; for example rkothadi / broadcptac / etc"
    echo "-f | <fire_ns>            | string | Firecloud Namespace"
    echo "-m | <modify_from_docker> | string | Build a new docker from this base docker when there are shared files"
    echo "-t | <docker_tag>         | number | Docker Tag"
    echo "-r | <recreate_dir>       | flag   | Destroy original and recreate directory structure"
    echo "-b | <build_docker>       | flag   | Build and Push docker"
    echo "-u | <common_util>        | flag   | Add utilities, packages, and data to docker directory"
    echo "-h | <help_flag>          | flag   | Print Usage"
}

# Points to PGDAC/
cd ..
wspace=`pwd`
cd hydrant

while getopts ":s:n:i:w:d:f:m:t:o:hrbu" opt; do
    case $opt in 
        s) sub_dir="$OPTARG"
           ;;
        n) wf_name="$OPTARG"
           ;;
        i) docker_file_temp="$OPTARG"
           ;;
        o) docker_file_notemp="$OPTARG"
           ;;
        w) wdl="$OPTARG"
           ;; 
        d) docker_ns="$OPTARG"
           ;;
        f) fire_ns="$OPTARG"
           ;;
        m) modify_docker_from="$OPTARG"
           ;;
        t) docker_tag="$OPTARG"
           ;;
        r) recreate_dir="yes"
           ;;
        b) build_docker="yes"
           ;;
        h) help_flag="yes"
           ;;
        u) common_util="yes"
           ;;
        \?) echo "Invalid Option -$OPTARG" >&2
           ;;
    esac
done

# Check if <sub_dir> and <wf_name> exist. Exit otherwise.
if [ "$sub_dir" = "" ] || [ "$wf_name" = "" ] || [[ ! -z "$help_flag" ]]; then
    display_usage
    echo "Program exiting..."
   exit
fi

cd $sub_dir

# If recreate_flag exists destroy and recreate hydrant structure
# else copy only source files to docker dir
if [[ ! -z "$recreate_dir" ]]; then
    if [ -d "$wf_name" ]; then
        read -p "Directory with the same workflow name exists. Erase and create new? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && rm -r $wf_name || exit 1
    fi

    echo -e "[Task $wf_name]\nSrc=$wspace/src/$wf_name" > $wf_name"-config.txt"
    hydrant init $wf_name -c $wf_name"-config.txt" -n 0
else
    echo "Copying source files to docker dir..."
    cp $wspace/src/$wf_name/* $wf_name/$wf_name/src/.
fi

# Enter <wf_name>/<wf_name>/ 
# 0. set <docker_file> from the correct variable
# 1. Copy dockerfile to <wf_name>/<wf_name> IF <docker_file> exists
# 2. Create data, packages, r-util subdirectories and copy files within IF <modify_from_docker> is unset
# 3. IF <modify_from_docker> is SET, modify existing Dockerfile's base docker
# 4. Build and push docker to dockerhub based on if <docker_file> and / or <docker_tag> is set

if [[ ! -z "$docker_file_temp" ]]; then
    docker_file=$wspace/hydrant/$docker_file_temp
fi

if [[ ! -z "$docker_file_notemp" ]]; then
    docker_file=$docker_file_notemp 
fi

cd $wf_name/$wf_name

if [[ ! -z "$docker_file" ]]; then
   cp $docker_file Dockerfile
fi

# If common_util flag exists copy the utilties to docker dir
if [[ ! -z "$common_util" ]]; then
    echo "Copying Data, Packages, R-utilities to docker dir..."
    mkdir -p data
    mkdir -p packages
    mkdir -p R-utilities
    cp $wspace/data/* data/.
    cp $wspace/hydrant/packages/* packages/.
    cp $wspace/hydrant/r-util/* R-utilities/.
    echo $'!data\n!packages\n!R-utilities' >> .dockerignore
fi

# Edit docker base in dockerfile
if [[ ! -z "$modify_docker_from" ]]; then
    sed -i '' "s|broadcptac/r-util:2|broadcptac/$modify_docker_from:1|g" Dockerfile
    sed -i '' "s|broadcptac/pgdac_common:1|broadcptac/$modify_docker_from:1|g" Dockerfile
fi

# Build & Push Docker
if [[ ! -z "$build_docker" ]]; then
    if [[ -z "$docker_tag" ]]; then
        hydrant build -n $docker_ns
        hydrant push -n $docker_ns
    else
        docker build -t $docker_ns/$wf_name:$docker_tag .
        docker push $docker_ns/$wf_name:$docker_tag
    fi
fi

cd ..


# Overwrite WDL template with your own
if [ ! -z "$wdl" ]; then
    cp $wdl $wf_name.wdl
    hydrant validate
fi

cd - > /dev/null
