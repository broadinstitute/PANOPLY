#! /usr/bin/env bash

USAGE="Usage: "`basename $0`" [-x <file/pattern to exclude> [-x <file/pattern to exclude> ...]] <package name>"

EXCLUDE_STR='"fc-[a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9]/*" lost+found/\* "tmp.[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]/*" exec.sh'

while getopts 'x:h' opt; do
  case $opt in
    x) EXCLUDE_STR="$EXCLUDE_STR '$OPTARG'";;
    h) echo $USAGE; exit 0;;
   \?) echo $USAGE 1>&2; exit 1;;
    esac
done

shift $((OPTIND-1))

if [ $# -ne 1 ]; then
  echo "Expected one positional argument, found $#" 1>&2
  echo $USAGE 1>&2
  exit 1
fi

PACKAGE=$1

eval "zip -r $PACKAGE . -x $EXCLUDE_STR"
