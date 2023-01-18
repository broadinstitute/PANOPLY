#!/bin/bash
ssid=$1
ssid=${ssid##*/}
ssid=${ssid%.zip}
src=`pwd`
mkdir -p psms
psms="$src/psms"

( cd $ssid;
  find "." -name "psm.tsv" -not -path '*/\.*' > $src/allpsms.txt )

cat $src/allpsms.txt

while IFS= read -r psm; do
  IFS="/" read -ra psm_split <<< "$psm"
  echo "${psm[@]}" >> $src/print.log
  newpsm=${psm_split[1]}"_"${psm_split[2]}
  cp $ssid/$psm $src/psms/$newpsm
done < $src/allpsms.txt

zip psms.zip psms/*.tsv
