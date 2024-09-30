#!/bin/bash

echo> output3.txt
for file in GCA*.fna
do
   echo "Processing $file"
   echo $file >> output3.txt
   python trans_dna.py $file >> output3.txt
done
