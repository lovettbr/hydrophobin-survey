#!/bin/bash

folder=$1
out=$2

files=$(ls $folder)

mkdir "$out"

for sequences in $files
do

name=$(echo "$sequences" | cut -f 1 -d '.')

hmmsearch --tblout "$out"/"$name"_pfam.out ../Profiles/Hydrophobins.hmm "$folder"/$sequences

awk '/PF01185.19/{ print $1 }' "$out"/"$name"_pfam.out  > "$out"/"$name"_hyd1.txt
seqtk subseq "$folder"/$sequences "$out"/"$name"_hyd1.txt > "$out"/"$name"_hyd1.faa

awk '/PF06766.12/{ print $1 }' "$out"/"$name"_pfam.out > "$out"/"$name"_hyd2.txt
seqtk subseq "$folder"/$sequences "$out"/"$name"_hyd2.txt > "$out"/"$name"_hyd2.faa

done
