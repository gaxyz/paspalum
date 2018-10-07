#!/usr/bin/env bash

mithoGenome = $1
chloroGenome = $2


mkdir chloroplastAssemblbyDB

cat $mithoGenome > chloroplastAssemblbyDB/cAssemblyRef.fna
cat $chloroGenome >> chloroplastAssemblbyDB/cAssemblyRef.fna


cd chloroplastAssemblbyDB 

makeblastdb -in cAssemblyRef.fna -parse_seqids -dbtype nucl -title cAssemblyRef -out cAssemblyRef

