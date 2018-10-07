#!/usr/bin/env python3
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
               


