#!/bin/bash
#
# Copyright (c) 2020 The Broad Institute, Inc. All rights reserved.
#
start_dir=`pwd`    # must be invoked from the PANOPLY/release directory
panoply=$start_dir/..
red='\033[0;31m'
grn='\033[0;32m'
reg='\033[0m'
err="${red}Error:${reg}"
not="${grn}====>>${reg}" ## notification

## users (for method permissions)
proteomics_comp=(manidr@broadinstitute.org nclark@broadinstitute.org wcorinne@broadinstitute.org)
proteomics_cptac=(GROUP_Broad_CPTAC@firecloud.org)

## documentation location (assumes wiki repo is available in path)
doc_dir="$panoply/../PANOPLY.wiki"
generic_docs="PANOPLY-Tutorial.md Navigating-Results.md PANOPLY-without-Terra.md"


display_usage() {
  echo "usage: ./setup-release.sh -p pull_dns -T release_tag -N release_dns "
  echo "                          -w workspace_all -y workspace_pipelines -r project [-h]"
  echo "-N | string | Docker namespace for building and pushing release"
  echo "-T | string | Docker release tag"
  echo "-p | string | Docker namespace to pull latest dev images from"
  echo "-r | string | Workspace project id"
  echo "-w | string | Workspace to populate with all methods"
  echo "-y | string | Workspace to populate pipeline/workflow methods"
  echo "-P | string | Patch flag. When toggled, only specified modules will be rebuilt."
  echo "-R | string | Patch flag. When toggled, dockers will be rebuilt for specified modules."
  echo "        [...] Trailing arguments will be interpretted as modules to patch-fix. Should only be used alongside -P flag."
  echo "-h | flag   | Print Usage"
  exit
}

# replace the docker namespace and docker tag in the existing WDL
# to the specified docker namespace and docker tag
replaceDockerInWdl() {
  task=$1
  new_ns=$2
  new_tag=$3
  old_ns=$4

  # wdl_dns=`grep "/$task:" $task.wdl | cut -d'"' -f2 | cut -d'/' -f 1`
  # wdl_tag=`grep "/$task:" $task.wdl | cut -d'"' -f2 | cut -d':' -f 2`
  # sed -i '' "s|$wdl_dns/$task:$wdl_tag|$new_ns/$task:$new_tag|g" $task.wdl
  
  # need to handle WDLs where the docker name(s) do not match $task
  sed -i'' -e "s|\(docker.*\)\"$old_ns/\(panoply_.*\):.*\"|\1\"$new_ns/\2:$new_tag\"|g" $task.wdl
  rm $task.wdl-e
}


replaceWdlImports() {
  wdl_f=$1
  repl=$2
  modules=$( grep import $wdl_f | sed -n 's|.*:\(panoply_.*\)/versions.*|\1|p' )
  for mod in $modules
  do
    new_snap=$( grep ":$mod/" $repl )
    sed -i'' -e "s|https.*:$mod/.*descriptor|$new_snap|" $wdl_f
  done
  rm $wdl_f-e
}


put_method_config() {
  w=$1
  p=$2
  m=$3

  exists=`fissfc config_list -w $w -p $p | cut -f 2 | grep "^$m$" | wc -l`
  if [ $exists -eq 1 ]; then
    fissfc config_delete -c $m -w $w -p $p -n $release_dns
  fi

  fissfc config_put -w $w -p $p -c $m-template.json 
}


configure_primary_workflow() {
  ws=$1  # workspace
  wp=$2  # project
  wf=$3  # workflow
  
  bucket=`fissfc space_info -w $ws -p $wp | jq -r '.workspace.bucketName'`
  
  if [ "$wf" == "panoply_main" ]; then
    # configure panoply_main for one-click execution using data model metadata
    cat $wf-template.json |  \
      jq '.inputs."panoply_main.yaml" = $val' --arg val "this.parameters" |  \
      jq '.inputs."panoply_main.input_cna" = $val' --arg val "this.cna_ss" |  \
      jq '.inputs."panoply_main.input_rna" = $val' --arg val "this.rna_ss" |  \
      jq '.inputs."panoply_main.sample_annotation" = $val' --arg val "this.annotation_ss" |  \
      jq '.inputs."panoply_main.groups_file" = $val' --arg val "this.groups_ss" |  \
      jq '.inputs."panoply_main.run_cmap" = $val' --arg val "\"false\"" |  \
      jq '.inputs."panoply_main.annotation_pathway_db" = $val' --arg val "this.gseaDB" |  \
      jq '.inputs."panoply_main.geneset_db" = $val' --arg val "this.gseaDB" |  \
      jq '.inputs."panoply_main.ptm_db" = $val' --arg val "this.ptmseaDB" > new-template.json
    mv new-template.json $wf-template.json  
    put_method_config $ws $wp $wf
    
    # in addition to a generic workflow include versions for proteome and phosphoproteome
    for t in proteome phosphoproteome 
    do
      cat $wf-template.json |  \
        jq '.name = $val' --arg val "${wf}_${t}" |  \
        jq '.inputs."panoply_main.ome_type" = $val' --arg val "\"${t}\"" |  \
        jq '.inputs."panoply_main.input_pome" = $val' --arg val "this.${t}_ss" > ${wf}_${t}-template.json
      put_method_config $ws $wp ${wf}_${t} 
    done
  fi
  
  if [ "$wf" == "panoply_unified_workflow" ]; then
    # configure panoply_unified_workflow for one-click execution using data model metadata
    cat $wf-template.json |  \
      jq '.inputs."panoply_unified_workflow.yaml" = $val' --arg val "this.parameters" |  \
      jq '.inputs."panoply_unified_workflow.cna_data" = $val' --arg val "this.cna_ss" |  \
      jq '.inputs."panoply_unified_workflow.rna_data" = $val' --arg val "this.rna_ss" |  \
      # jq '.inputs."panoply_unified_workflow.sample_annotation" = $val' --arg val "this.annotation_ss" |  \
      jq '.inputs."panoply_unified_workflow.groups_file" = $val' --arg val "this.groups_ss" |  \
      jq '.inputs."panoply_unified_workflow.run_cmap" = $val' --arg val "\"false\"" |  \
      jq '.inputs."panoply_unified_workflow.run_nmf" = $val' --arg val "\"true\"" |  \
      jq '.inputs."panoply_unified_workflow.run_ptmsea" = $val' --arg val "\"false\"" |  \
      jq '.inputs."panoply_unified_workflow.nmf.gene_set_database" = $val' --arg val "this.gseaDB" |  \
      jq '.inputs."panoply_unified_workflow.pome.geneset_db" = $val' --arg val "this.gseaDB" |  \
      jq '.inputs."panoply_unified_workflow.pome.ptm_db" = $val' --arg val "this.ptmseaDB" |  \
      jq '.inputs."panoply_unified_workflow.acetyl_ome" = $val' --arg val "this.acetylome_ss" |  \
      jq '.inputs."panoply_unified_workflow.prote_ome" = $val' --arg val "this.proteome_ss" |  \
      jq '.inputs."panoply_unified_workflow.phospho_ome" = $val' --arg val "this.phosphoproteome_ss" |  \
      jq '.inputs."panoply_unified_workflow.ubiquityl_ome" = $val' --arg val "this.ubiquitylome_ss" |  \
      jq '.inputs."panoply_unified_workflow.pome.annotation_pathway_db" = $val' --arg val "this.gseaDB" > new-template.json
    mv new-template.json $wf-template.json  
    put_method_config $ws $wp $wf
  fi
  
  if [ "$wf" == "panoply_nmf_workflow" ]; then
    # configure panoply_nmf_workflow with template for paired-array input
    cat $wf-template.json |  \
      jq '.inputs."panoply_nmf_workflow.ome_pairs" = $val' --arg val '[ { "left":"prot", "right":this.proteome_ss }, { "left":"pSTY", "right":this.phosphoproteome_ss }, { "left":"acK", "right":this.acetylome_ss }, { "left":"ubK", "right":this.ubiquitylome_ss }, { "left":"RNA", "right":this.rna_ss }, { "left":"CNA", "right":this.cna_ss } ]' \
      > new-template.json
    mv new-template.json $wf-template.json  
    put_method_config $ws $wp $wf
  fi
} 


installMethod() {
  meth=$1
  meth_wdl=$2
  type=$3   # tasks or workflows
  
  echo -e "$not ... Installing method $meth"

  # install method on FireCloud/Terra with synopsis/documentation
  # markdown documentation does not render in Terra -- use a URL to point to
  #  the release-specific documentation on GitHub (added first line)
  syn=""
  doc=""
  if [[ -f $panoply/hydrant/$type/$meth/$meth-synopsis.txt ]]; then
    syn=`cat $panoply/hydrant/$type/$meth/$meth-synopsis.txt`
  fi
  orig_doc=`ls $doc_dir/*$meth.md`
  if [[ -f $orig_doc ]]; then
    cp $orig_doc $meth.md    
    # documentation too large error for some modules -- hence just use URL for method
    echo -e "Documentation at https://github.com/broadinstitute/PANOPLY/blob/$release_ver/release/$release_dir/$meth/$meth.md\n" > $meth-doc-URL.txt
    fissfc meth_new -m $meth -n $release_dns -d $meth_wdl -c "Snapshot for Release v$release_tag" -s "$syn" --doc $meth-doc-URL.txt
  else
    fissfc meth_new -m $meth -n $release_dns -d $meth_wdl -c "Snapshot for Release v$release_tag" -s "$syn"
  fi
    
  # get method snapshot id and save release snapshot ids
  snap=$(fissfc meth_list -m $meth -n $release_dns | sort -n -k3 | tail -1 | cut -f3)

  # set method permissions
  fissfc meth_set_acl -m $meth -n $release_dns -i $snap -r OWNER --users ${proteomics_comp[@]}
  fissfc meth_set_acl -m $meth -n $release_dns -i $snap -r READER --users public
  
  # add to snapshot list
  echo "https://api.firecloud.org/ga4gh/v1/tools/$release_dns:$meth/versions/$snap/plain-WDL/descriptor" >> $snapshots
  
  # copy method to workspace(s)
  #  there is no (?) direct way to copy a method to a workspace
  #  instead, create an empty config template and install it in the workspace
  #  (delete method config if it exists)
  fissfc config_template -m $meth -n $release_dns -i $snap -t sample_set |  \
    sed 's/\"EDITME.*\"/""/' | jq '.name = $val' --arg val $meth > $meth-template.json
  
  if [[ ("$type" = "workflows") && ("$meth" = "panoply_main" || "$meth" = "panoply_unified_workflow") ]]; then
    # only panoply_main and panoply_unified_workflow in $wkspace_pipelines
    configure_primary_workflow $wkspace_pipelines $project $meth
  elif [[ ("$meth" = "panoply_nmf_workflow") ]]; then
    # also configure panoply_nmf_workflows to have paired-array template
    configure_primary_workflow $wkspace_all $project $meth
  else
    put_method_config $wkspace_all $project $meth
  fi
}


while getopts ":p:T:N:r:w:y:PRh" opt; do
    case $opt in
        p) pull_dns="$OPTARG";;
        T) release_tag="$OPTARG";;
        N) release_dns="$OPTARG";;
        r) project="$OPTARG";;
        w) wkspace_all="$OPTARG";;
        y) wkspace_pipelines="$OPTARG";;
        P) patch_flag=TRUE;;
        R) rebuild_docker_flag=TRUE;;
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

if [[ -z $project ]]; then
  echo -e "$err Workspace project ID required (usually broad-firecloud-cptac). Exiting.."
  exit 1
fi

if [[ -z $wkspace_all ]]; then
  echo -e "$err Production workspace name for all methods required. Exiting.."
  exit 1
fi

if [[ -z $wkspace_pipelines ]]; then
  release_dns=$pull_dns
fi

if [[ -n $patch_flag ]]; then
  # parse modules argument
  shift "$(( OPTIND - 1 ))" # discard arguments parsed by getopts
  modules=($@) # parse remaining arguments as modules
  
  #validate modules argument
  if [[ -z $modules ]]; then
    echo -e "$err No modules specified. Exiting.."
    exit 1
  fi
  modules_all=( $( ls -d $panoply/hydrant/tasks/panoply_* | xargs -n 1 basename ) )
  for mod in "${modules[@]}"; do # check that all input modules are valid tasks
    if ! [[ ${modules_all[@]} =~ $mod ]]; then
      echo -e "$err Invalid module $mod. Module must match a task from ../hydrant/tasks"
      exit 1
    fi
  done

  echo -e "$not Patching module(s) ${modules[@]} for release-$release_tag."
  echo -e "$not All workflows will be rebuilt."
  if [[ $rebuild_docker_flag ]]; then # if we are rebuilding dockers
    echo -e "$not Dockers will be rebuilt."
  fi
  echo -e "$not Documentation will be rereleased."
fi




## DOCKER HUB
echo -e "$not Logging in to Docker Hub"
base_url="https://registry.hub.docker.com/v2/repositories/"
docker login


## WORKSPACES 
createWkSpace() {
  ws=$1
  
  if fissfc space_exists -w $ws -p $project -q; then
    echo -e "$not Creating workspace $ws"
    fissfc space_new -w $ws -p $project
  else
    echo -e "$not Workspace $ws exisits. Updating permissions"
  fi
  # set permissions
  fissfc space_set_acl -w $ws -p $project -r OWNER --users ${proteomics_comp[@]}
  # setting workspace to PUBLIC is currently not functional -- need to contact Terra support
  # fissfc space_set_acl -w $ws -p $project -r READER --users public
  fissfc space_set_acl -w $ws -p $project -r READER --users ${proteomics_cptac[@]}
  
  # set workspace public using FireCloud API -- eg
  # curl -X POST "https://api.firecloud.org/api/methods/broadcptac/panoply_main_copy/1/permissions" \
  #      -H "accept: application/json" -H "Authorization: Bearer \
  #         $(gcloud auth --account=manidr@broadinstitute.org print-access-token)" \
  #      -H "Content-Type: application/json" -d "[{\"user\":\"public\",\"role\":\"READER\"}]"
  
}
# workspace for pipelines + all modules
createWkSpace $wkspace_all 
# workspace for pipelines only
createWkSpace $wkspace_pipelines 


## DIRECTORY CREATION / CLEANUP
release_dir=version-$release_tag
release_ver=release-$release_tag
snapshots="$panoply/release/$release_dir/snapshot-ids.txt"

if [[ -z $patch_flag ]]; then # during full release
  mkdir -p $release_dir # create release directory
  rm -f $snapshots # delete all snapshots
  yes | docker system prune --all
  modules=( $( ls -d $panoply/hydrant/tasks/panoply_* | xargs -n 1 basename ) ) # rebuild all modules
fi # otherwise, modules are pulled from trailing arguments ( see display_usage() )

for mod in "${modules[@]}"
do
  echo -e "$not Processing task $mod"

  if [[ -n $patch_flag ]]; then # if we are patch-fixing
    sed -i'' -e "/:$mod\/versions/d" $release_dir/snapshot-ids.txt # delete relevant snapshots from snapshot-ids.txt
    rm -r ./$release_dir/$mod # delete relevant version folders
  fi
  
  if [[ $rebuild_docker_flag ]]; then # if we are rebuilding dockers
    wd=`pwd` # record og dir
    cd $panoply/hydrant # change dir to location of setup.sh
    ./setup.sh -t $mod -n $pull_dns -y -b -u -x # rebuild
    ./setup.sh -t $mod -z    # cleanup
    cd $wd # change back to og dir
  fi

  url=$base_url$pull_dns/$mod/tags
  lat=( $( curl -s -S "$url" | \
             jq '."results"[]["name"]' | \
             sed -n 1p | cut -d'"' -f2 ) )
  if [[ -z $lat ]]; then
    lat="NO_DOCKER"
  else
    # this module present -- set tag to "latest"
    lat="latest"
  fi
  
  mkdir -p $release_dir/$mod
  cd $release_dir/$mod
  
  if [ "$lat" != "NO_DOCKER" ]; then
    # build release docker
    echo -e "FROM $pull_dns/$mod:$lat" > Dockerfile
    docker build --rm --no-cache -t $release_dns/$mod:$release_tag .
    docker images | grep "$mod"
    docker push $release_dns/$mod:$release_tag
  fi
  
  # copy and update task WDL, install method and save snapshot id
  wdl_dir=$panoply/hydrant/tasks/$mod/
  wdl=$mod.wdl
  if [[ -f $wdl_dir/$wdl ]]; then
    cp $wdl_dir/$wdl $wdl
    replaceDockerInWdl $mod $release_dns $release_tag $pull_dns
    installMethod $mod $wdl "tasks"
  fi
  cd $start_dir
done


## WORKFLOWS
workflows=( $( ls -d $panoply/hydrant/workflows/panoply_* | xargs -n 1 basename |  grep -vE 'panoply_nmf_workflow|panoply_main|panoply_unified_workflow' ) ) # remove panoply_nmf_workflow, panoply_main, and panoply_unified_workflow
workflows+=( panoply_nmf_workflow panoply_main panoply_unified_workflow ) # add to end of array, to ensure that these are built last
for wk in "${workflows[@]}"
do
  echo -e "$not Processing workflow $wk"

  if [[ -n $patch_flag ]]; then # if we are patch-fixing
    sed -i'' -e "/:$wk\/versions/d" ./$release_dir/snapshot-ids.txt # delete relevant snapshots from snapshot-ids.txt
    rm -r ./$release_dir/$wk # delete relevant version folders
  fi

  mkdir -p $release_dir/$wk
  cd $release_dir/$wk

  # copy and update WDL imports 
  wdl_dir=$panoply/hydrant/workflows/$wk/
  wdl=$wk.wdl
  if [[ -f $wdl_dir/$wdl ]]; then
    cp $wdl_dir/$wdl $wdl
    replaceWdlImports $wdl $snapshots
    installMethod $wk $wdl "workflows"
  fi
  cd $start_dir
done


## DOCUMENTATION
if [[ -z $patch_flag ]]; then # only make directory on initial-release
  mkdir -p $release_dir/docs
fi
for f in $generic_docs
do
  cp $doc_dir/$f $release_dir/docs/$f
done
