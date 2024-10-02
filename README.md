# gene-finder

## 1. Code
The gene finder tool is implemented in [gene_finder.py](https://github.com/manalalshamrani/gene-finder/blob/master/gene_finder.py).
The results file is q1.txt.

## 2. Code for 1,2 combined
The extended code to include reverse complments is implmented in [gene_finder_2.py](https://github.com/manalalshamrani/gene-finder/blob/master/gene_finder_2.py).
The results file is output.txt.

## 3. Code for Rosalind
The Rosalind code is implmented in [trans_dna.py](https://github.com/manalalshamrani/gene-finder/blob/master/trans_dna.py).
The results file is output2.txt. 

## 4. 14: screenshot of command line
```bash
#!/bin/bash                                                                                                                                                                                                 

echo> output3.txt
for file in GCA*.fna
do
   echo "Processing $file"
   echo $file >> output3.txt
   python trans_dna.py $file >> output3.txt
done
```
The results file is output3.txt. 

## 5. Code for 5 updated 
                                                                                                                                                                                               
```bash
#!/bin/bash  
echo "" > output4.txt
for file in GCA*.fna
do
   echo "Processing $file"
   echo $file >> output4.txt
   python trans_dna_filtered.py -f $file -l 100  >> output4.txt
done
```

The  code to include length is implmented in [trans_dna_filtered.py](https://github.com/manalalshamrani/gene-finder/blob/master/trans_dna_filtered.py).
The results file is output4.txt. 


## 6. Cofe for 6 updated
The  code to include Shine-Dalgarno sequence is implmented in [trans_dna_filtered2.py](https://github.com/manalalshamrani/gene-finder/blob/master/trans_dna_filtered2.py).

Analysis done on: GCA_000007125.1_ASM712v1_genomic.fna
```bash
[alshammm@login509-02-l data]$ wc -l /home/alshammm/output.q6.txt
115 /home/alshammm/output.q6.txt
```



#
This HW was done with the help of Dalia, Layan, Yazeed and Haoling.
BlackBox AI was used to debugg code.
