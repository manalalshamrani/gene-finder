echo "" > output5.txt
for file in GCA*.fna
do
   echo "Processing $file"
   echo $file >> output5.txt
   python trans_dna_filtered2.py -f $file -l 100  >> output5.txt
done
