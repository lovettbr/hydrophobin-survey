#!/bin/bash

files=$(ls . | grep "faa")

mkdir renamed

for sequences in $files
do

name=$(echo "$sequences" | cut -f 1 -d '.')

sed "s|>|>${name}_|g" $sequences > renamed/"$name".faa

done
