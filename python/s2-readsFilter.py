#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import subprocess
from Bio import SeqIO


def interleave_pairs(R1, R2, outPairs, outOrphans):

    input_forward_filename = R1
    input_reverse_filename = R2
    output_pairs_filename = outPairs
    output_orphans_filename = outOrphans

    #asumo que no tienen suffix

    print("Indexing forward file...")
    forward_dict =  SeqIO.index( input_forward_filename, "fastq" ) 
    print("Indexing reverse file...")
    reverse_dict = SeqIO.index( input_reverse_filename, "fastq" )
    print("Outputing pairs and forward reads...")
    with open(output_pairs_filename, "wb" ) as pair_handle:
        with open(output_orphans_filename, "wb"  ) as orphan_handle:

            for key in forward_dict:
                if key in reverse_dict:
                    pair_handle.write( forward_dict.get_raw( key ) )
                    pair_handle.write( reverse_dict.get_raw( key ) )
                else:
                    orphan_handle.write(forward_dict.get_raw( key ) )
    
    print( "Outputting reverse orphans..." )
    with open( output_orphans_filename, 'ab' ) as orphan_handle:
        for key in reverse_dict:
            if key not in forward_dict:
                orphan_handle.write( reverse_dict.get_raw( key ) )


    print("Done!")
             


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

helpMessage = '''
############################################Mixed chloroplast assembly###############################################

This script filters reads with homology to specified chloroplast and mithochondria databases using BLAST.



Assumes installed software(added to path):


-blast
-seqtk



Usage: <this>.py <filterAssemblyblastDBDir>  <outDir> <sample_fastq_list>


fastq files can't be compressed
Absolute paths are necessary

######################################################################################################################
'''



if sys.argv[1] == "-h":
    eprint(helpMessage)
    exit()


#-----------------------Setting things up--------------------------------------------------

cAssemblyRef = sys.argv[1]               # directorio de base de datos de blast de mitocondrias
outputDirectory = sys.argv[2]        # directorio donde contener todos los archivos de salida
sampleList = sys.argv[3]
sampleDir = sys.argv[4]


#---------------
#cAssemblyRef = "DB/mixedDB/cAssemblyRef"              
#outputDirectory = "test/chloroReads/"        
#sampleList = "samples.txt"
#sampleDir = "test/1_trimmomatic/"
#sample = "HKYCTBBXX_Paspalum1_S1_L001"
   


cAssemblyRef = os.path.abspath( cAssemblyRef )
outputDirectory = os.path.abspath( outputDirectory )
sampleList = os.path.abspath( sampleList )
sampleDir = os.path.abspath( sampleDir ) 

#sample list example
'''
sample_name R1-name R2-name unpaired-name
.
.
.



'''

identity_cutoff = 90    #cutoff de identidad para blast
overlap_cutoff = 0.98   #cutoff de overlap para blast 
blastThreads = "12"        # numero de nucleos a usar en blast
mainDir = os.getcwd()   # grabar directorio donde fue llamado el script
os.chdir( sampleDir )



with open(sampleList, 'r') as samples:

    for sample in samples:
        

        sampleName = sample.split()[0]
                     

        r1FQ = sampleName + "_R1" + ".fq.gz"
        r2FQ = sampleName + "_R2" + ".fq.gz"
        upFQ = sampleName + "_U" + ".fq.gz"

        eprint( "# Reading sample {0} FQFiles: {1}, {2}, {3}".format( sampleName, r1FQ, r2FQ, upFQ ) )

        readtypes = [ "_R1", "_R2", "_U" ]
        eprint( "## Decompressing and converting to fasta format" )
        for suffix in readtypes:
            
            decompress = ["pigz", "-d","-k","-p","4",sampleName + suffix + ".fq.gz" ]
            subprocess.call( decompress )


            fastaTransform = "seqtk seq -a " + sampleName + suffix + ".fq" + ">"  + sampleName + suffix + ".fna" 
            subprocess.call( fastaTransform , shell=True )





        #Identificar reads a filtrar--------------------------------------------------------------------------------------------------------

        

        eprint( "## Looking for chloroplastic and mithochondrial reads..." )

        for readtype in readtypes:  # para cada r1, r2, unpaired hace un blast y genera una lista de reads seleccionados


            outReadsSet = set()   # inicializar set de reads de interes (homologia con mito o cloro)
            suffix = readtype 
            fileName = sampleName + suffix + ".fna"

            blastOutName = "{0}_cAssemblyReads{1}.blastn.out.tmp".format( sampleName, suffix )               # defino nombre de output de blast 

            outReadsListFilename = "{0}_cAssemblyReadsList{1}.tsv.tmp".format( sampleName, suffix ) # defino nombre de archivo de salida

            blastFormat = '"6 qseqid sseqid qstart qend sstart send  qlen slen length pident evalue"'  # formato de output de blast personalziado
            #################    0     1      2     3     4     5     6    7     8      9     10     # corresponde a c/u



            blastCommand = [ "blastn" , "-query" , fileName , "-db" , cAssemblyRef ,  "-out" ,  blastOutName ,
                            "-evalue", "1e-10" , "-dust" , "no" , "-num_threads" , blastThreads , "-max_target_seqs", "1",
                            "-best_hit_score_edge" , "0.05" , "-best_hit_overhang" , "0.25" ,
                            "-perc_identity" , "50" ,"-outfmt" , blastFormat ]
            
            

            string = ""
            for i in blastCommand:
                string = string + " " + i


            # defino el comando de blast a usar
            eprint("### Blasting...{0}".format(fileName) )
                   
                   
                   
                   
            subprocess.call( string, shell = True ) # empiezo blast

            eprint("#### Reading BLAST output...")
            with open( blastOutName, 'r' ) as blastOutfile:             # leyendo el archivo de output de blast
                
                for line in blastOutfile:          

                    qseqid, sseqid, qstart, qend, sstart, send, qlen, slen, length, pident, evalue = line.split()
                    overlap = int(qlen) / int(length)
                    
                    if overlap >= overlap_cutoff and float(pident) >= identity_cutoff :   
                        outReadsSet.add( qseqid )
                        
            eprint("#### Printing reads into {0}".format( outReadsListFilename ) )

            with open( outReadsListFilename, 'w' ) as output:

                for item in outReadsSet:
                    output.write( item + "\n" )   # imprimo set de interes





        #-----------------------------------------------------------------------------------------------------------------------------------------





		#Filtrado de reads de interes----------------------------------------------------------------------------------------------------
		######De cada lista r1, r2, unpaired, a partir de sus respectivos fq, 
		######selecciona los reads y los ordena por r1,r2 / huerfano  

        eprint("## Filtering and sorting reads from original fastq files...")

		#aca agarro los fastq y los parseo segun su id
		# si su id es uno de los reads que tienen homologia con cloroplastos


        for readtype in readtypes: # para cada tipo


            suffix = readtype

            fileName = sampleName + suffix + ".fq"

            outReadsListFilename = "{0}_cAssemblyReadsList{1}.tsv.tmp".format( sampleName, suffix )

            eprint("### Filtering...")
            outReadsFilename = "{0}_FilteredReads{1}.fastq".format(sampleName, suffix)
            filterCommand =  "seqtk subseq {0} {1} > {2}".format(fileName,outReadsListFilename,outReadsFilename) # defino comando de seqtk

            subprocess.call( filterCommand , shell = True ) # esto no lo tendria que hacer... pero lo permito por velocidad, supongo


        for readtype in readtypes:


            rm_fna = [ "rm" , sampleName + suffix + ".fna" ]
            subprocess.call(rm_fq)

		# tengo que hacer interleaving e identifiacion de reads huerfanos
		#interleave_pairs tiene que estar incluido en este script ?

        
        outReadsFilenameR1 = "{0}_FilteredReads_R1.fastq".format( sampleName )
        outReadsFilenameR2 =  "{0}_FilteredReads_R2.fastq".format( sampleName )
        

        pairedInterlFileName = "{0}_PairedInterleaved.fastq".format(sampleName)
        orphanFileName = "{0}_Orphan.fastq" 

        eprint("### Interleaving and getting orphans...")

        interleave_pairs( outReadsFilenameR1, outReadsFilenameR2, pairedInterlFileName, orphanFileName )

        # me falta concatenar los unpaired con los orphans











   

