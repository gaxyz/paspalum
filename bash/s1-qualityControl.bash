#!/usr/bin/env bash
set -e
set -u
set -o pipefail
trimmomaticPath=~/Software/Trimmomatic-0.38
adapterFile=$trimmomaticPath/adapters/TruSeq3-PE.fa


originalDir=$PWD



#LISTA = le meto de entrada una lista con las muestras

LISTA=$(readlink -f $1) # lista de nombres de muestras
SAMPLEDIR=$2
nTHREADS=10
## 1 ## CORRO FASTQC EN RAW FASTQs



cd $SAMPLEDIR

#mkdir 0_raw_fastq
mkdir -p 1_trimmomatic
mkdir -p 2_trimmed_fastqc

cd 0_raw_fastq

while read -r line 
do

	I="$line"
       
        echo $I
        fastqc --threads $nTHREADS ${I}_R1.fastq.gz ${I}_R2.fastq.gz

done < "$LISTA" 

cd ..

#exit 1
#so far so good




## 2 ## CORRO MULTIQC PARA READS 1 y 2 

multiqc --title "MultiQC: FastQC on raw reads -READ#1" --ignore *2_fastqc.* --filename 0_multiqc_raw_fastqc_R1.html  0_raw_fastq

multiqc --title "MultiQC: FastQC on raw reads -READ#2" --ignore *1_fastqc.* --filename 0_multiqc_raw_fastqc_R2.html  0_raw_fastq


#exit 1
#so far so good

############aca

## 3 ## CORRO TRIMMOMATIC 0.36 EN RAW FASTQs ### Tengo que ver cuales son los adaptadores



cd 1_trimmomatic



while read -r line
do
	I="$line"
        echo "$I trimmomatic"
	java -jar $trimmomaticPath/trimmomatic-0.38.jar PE -threads $nTHREADS -phred33 ../0_raw_fastq/${I}_R1.fastq.gz ../0_raw_fastq/${I}_R2.fastq.gz ${I}_R1.fq.gz ${I}_1U.fq.gz ${I}_2.fq.gz ${I}_2U.fq.gz ILLUMINACLIP:$adapterFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

	cat ${I}_1U.fq.gz ${I}_2U.fq.gz > ${I}_U.fq.gz
	rm ${I}_1U.fq.gz ${I}_2U.fq.gz



done < "$LISTA"

cd ..

#exit 1
#so far so good


## 4 ## CORRO FASTQC EN TRIMMED FASTQs



while read -r line
do
	I="$line"
        echo $I
        fastqc --threads $nTHREADS --outdir 2_trimmed_fastqc 1_trimmomatic/${I}_R1.fq.gz 1_trimmomatic/${I}_R2.fq.gz
done < "$LISTA"

#exit 1
# so far so good

# 5 ## CORRO MULTIQC PARA READS 1 y 2 

multiqc --title "MultiQC: FastQC on Trimmomatic reads -READ#1" --ignore *2_fastqc.* --filename 2_multiqc_trimmed_fastqc_R1.html  2_trimmed_fastqc

multiqc --title "MultiQC: FastQC on Trimmomatic reads -READ#2" --ignore *1_fastqc.* --filename 2_multiqc_trimmed_fastqc_R2.html  2_trimmed_fastqc



