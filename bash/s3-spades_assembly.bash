#!/bin/usr/env bash

LISTA=$1
READSDIR=$2
OUTDIR=$3
THREADS=12




while read -r linea
do

	SAMPLE=$linea

	ILEAVED=${READSDIR}/${SAMPLE}_PairedInterleaved.fastq
	ORPHAN=${READSDIR}/${SAMPLE}_Orphan.fastq

	mkdir -p ${OUTDIR}/${SAMPLE}_ASSEMBLY
	

	spades.py  --12 $ILEAVED  -s $ORPHAN -t $THREADS -m 16  --careful  -o ${OUTDIR}/${SAMPLE}_ASSEMBLY


	
done <- "$LISTA"


