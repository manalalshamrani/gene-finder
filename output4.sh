echo "" > output4.txt
for file in GCA*.fna
do
   echo "Processing $file"
   echo $file >> output4.txt
   python trans_dna_filtered.py -f $file -l 100  >> output4.txt
done
