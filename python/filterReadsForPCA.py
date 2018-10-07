#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import subprocess
from Bio import SeqIO




def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


helpMessage= "Ha you thought there was help here, how naive. The world is a cruel place."


if sys.argv[1] == "-h":
    eprint(helpMessage)
    exit()


#-----------------------Setting things up--------------------------------------------------

cAssemblyRef = sys.argv[1]               # directorio de base de datos de blast
outputDirectory = sys.argv[2]        # directorio donde contener todos los archivos de salida
sampleName = sys.argv[3]                
sampleDir = sys.argv[4]
dbCodename = sys.argv[5]


cAssemblyRef = os.path.abspath( cAssemblyRef )
outputDirectory = os.path.abspath( outputDirectory )
sampleDir = os.path.abspath( sampleDir ) 


identity_cutoff = 90    #cutoff de identidad para blast
overlap_cutoff = 0.98   #cutoff de overlap para blast 
blastThreads = "4"        # numero de nucleos a usar en blast
mainDir = os.getcwd()   # grabar directorio donde fue llamado el script
os.chdir( sampleDir )
outFasta = "{0}_hits_{1}.fna".format(sampleName, dbCodename)

#------------------------------------------------------------------------------------------



r1FQ = sampleName + "_1" + ".fq"
r2FQ = sampleName + "_2" + ".fq"
upFQ = sampleName + "_U" + ".fq"

eprint( "# Reading sample {0} FQFiles: {1}, {2}, {3}".format( sampleName, r1FQ, r2FQ, upFQ ) )

readtypes = [ "_1", "_2", "_U" ]
eprint( "## Converting to fasta format" )
for suffix in readtypes:
    
    fastaTransform = "seqtk seq -a " + sampleName + suffix + ".fq" + ">"  + sampleName + suffix + ".fna" 
    subprocess.call( fastaTransform , shell=True )







#Identificar reads a filtrar--------------------------------------------------------------------------------------------------------



eprint( "## Looking for reads in BLAST DB..." )

for readtype in readtypes:  # para cada r1, r2, unpaired hace un blast y genera una lista de reads seleccionados


    outReadsSet = set()   # inicializar set de reads de interes (homologia con mito o cloro)
    suffix = readtype 
    fileName = sampleName + suffix + ".fna"

    blastOutName = "{0}_hits{1}.blastn.out.tmp".format( sampleName, suffix )               # defino nombre de output de blast 
    outReadsListFilename = "{0}_hits{1}.tsv.tmp".format( sampleName, suffix ) 
    

    blastFormat = '"6 qseqid sseqid qstart qend sstart send  qlen slen length pident evalue"'  # formato de output de blast personalziado
    #################    0     1      2     3     4     5     6    7     8      9     10     # corresponde a c/u



    blastCommand = [ "blastn" , "-query" , fileName , "-db" , cAssemblyRef ,  "-out" ,  blastOutName ,
                    "-evalue", "1e-10" , "-dust" , "no" , "-num_threads" , blastThreads , "-max_target_seqs", "1",
                    "-best_hit_score_edge" , "0.05" , "-best_hit_overhang" , "0.25" ,
                    "-perc_identity" , "50" ,"-outfmt" , blastFormat ]

    

    string = ""
    for i in blastCommand:
        string = string + " " + i


    
    eprint("### Blasting...{0}".format(fileName) )
           
           
           
           
    subprocess.call( string, shell = True ) # empiezo blast

    eprint("#### Reading BLAST output...")
    with open( blastOutName, 'r' ) as blastOutfile:             # leyendo el archivo de output de blast
        
        for line in blastOutfile:          

            qseqid, sseqid, qstart, qend, sstart, send, qlen, slen, length, pident, evalue = line.split()
            overlap = int(qlen) / int(length)
            
            if overlap >= overlap_cutoff and float(pident) >= identity_cutoff :   
                outReadsSet.add( qseqid )
                
    eprint("#### Printing read ids into {0}".format( outReadsListFilename ) )

    with open( outReadsListFilename, 'w' ) as output:

        for item in outReadsSet:
            output.write( item + "\n" )   # imprimo set de interes

    #cargar fasta entero para extraer 

    id_dict = {}
    for seq in SeqIO.parse( fileName, "fasta"):
        id_dict[seq.id] = seq.seq

    with open( outReadsListFilename , 'r') as hitList:
        with open( outFasta, 'a' ) as output:

            for line in hitList:
                identifier = line.split()[0]
                output.write( identifier + suffix + "\n" )
                output.write( str(id_dict[identifier]) + "\n" )



