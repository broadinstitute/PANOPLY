#!/bin/bash

src=`pwd`

createCopy(){
  cd $src
  rm -rf tasks-panoply
  mkdir -p tasks-panoply
  cp -r tasks/* tasks-panoply/.
}

renameFilesFolders(){
  cd $src/tasks-panoply/
  for dir in pgdac_*/;
  do
    old_dir=$dir
    new_dir="$(echo "$dir" | sed s/pgdac/panoply/)"

    ## rename the sub-task
    if [[ -d $old_dir/$old_dir ]]; then
      mv $old_dir/$old_dir $old_dir/$new_dir
    else
      echo "$old_dir/$old_dir does not exist."
    fi

    bn_old_dir="$(basename $old_dir)"
    bn_new_dir="$(basename $new_dir)"

    ## rename the WDLs
    if [[ -f $old_dir/$bn_old_dir.wdl ]]; then
      mv $old_dir/$bn_old_dir.wdl $old_dir/$bn_new_dir.wdl
    else
      echo "$old_dir/$bn_old_dir.wdl does not exist."
    fi

    ## rename the task folders
    mv $old_dir $new_dir
  done

  ## rename the exception
  mv ../tasks-panoply/panoply_mo_nmf/pgdac_mo_nmf_task_1 \
    ../tasks-panoply/panoply_mo_nmf/panoply_mo_nmf_task_1
}

mergeToTasks(){
  cd $src
  rm -rf tasks
  mv tasks-panoply tasks
}

createCopy
renameFilesFolders
mergeToTasks
