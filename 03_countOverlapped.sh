#!/bin/bash

sp=$SPECIES
input_dir="/data/sp11/05_mapping/isoformsID/$sp"
output_dir="/data/sp11/05_mapping/summary/$sp"
mkdir -p "$output_dir"

result20="${output_dir}/${sp}.depth20.overlap.txt"
result30="${output_dir}/${sp}.depth30.overlap.txt"
result40="${output_dir}/${sp}.depth40.overlap.txt"

python 03_countOverlapped.sh $input_dir 20 $result20
python 03_countOverlapped.sh $input_dir 30 $result30
python 03_countOverlapped.sh $input_dir 40 $result40
