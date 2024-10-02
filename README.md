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
Used to run 14 files:
```bash
#!/bin/bash                                                                                                                                                                              
       
for file in GCA*.fna
do
   echo "Processing $file"
   echo $file
   python gene_finder4.py -f $file
done
```
Analysis done on 14 files:
```bash
[alshammm@login509-02-l ~]$ wc -l output.q4.txt
765195 output.q4.txt
```

## 5. Code for 5 updated 
                                                                                                                                                                                         The  code to include length is implmented in [trans_dna_filtered.py](https://github.com/manalalshamrani/gene-finder/blob/master/trans_dna_filtered.py).
Analysis done on: GCA_000007125.1_ASM712v1_genomic.fna
  
```bash
[alshammm@login509-02-l data]$ wc -l /home/alshammm/output.q5.txt
29026 /home/alshammm/output.q5.txt
```


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
