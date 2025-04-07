#!/bin/bash

### This code extract all the assemblies that compose a pan genome graph (.gfa) in distict fasta files ###
# see PG-SCUnK webpage for details
# Developer : Tristan CUMER - t.cumer.sci[at]gmail.com

# Version : 1.0 - 05/04/2025

# Function to display help message
usage() {
    echo "Usage: $0 -p <panGenome> -t <tempDir> -o <outDir> -@ <threads>"
    exit 1
}

# Parse input arguments
while getopts "p:t:o:@:" opt; do
    case ${opt} in
        p) PG=${OPTARG} ;;       # panGenome file
        t) TEMPDIR=${OPTARG} ;;   # temporary directory
        o) OUTDIR=${OPTARG} ;;
        @) THREADS=${OPTARG} ;;   # number of threads
        *) usage ;;               # show usage if wrong input is given
    esac
done

# Check if all required parameters are provided
if [ -z "${PG}" ] || [ -z "${TEMPDIR}" ] || [ -z "${THREADS}" ]; then
    usage
fi

# Generate a random name
RandName=$(echo $RANDOM | md5sum | head -c 6)
echo "Random ID :" ${RandName}

# Create the temporary directory
# mkdir -p ${TEMPDIR}/${RandName}
mkdir -p ${OUTDIR}

#### Process the Pan Genome ####
echo "Extracting the multifasta from the graph"

# Extract the fasta with all the assemblies & index it
odgi paths -t ${THREADS} -f -i ${PG} > ${TEMPDIR}/${RandName}.full.fa

# Index the panGenome fasta using samtools
# Assuming samtools is available, loading necessary modules
samtools faidx ${TEMPDIR}/${RandName}.full.fa

#### Process the multifasta ####
echo "Splitting of the mulifasta per haplotype"

# Identify Haplotype names assuming PanSN formatting
cut -f 1 ${TEMPDIR}/${RandName}.full.fa.fai | awk -F"#" '{print $1"#"$2}' | sort | uniq > ${TEMPDIR}/${RandName}.ListHaplotypes.txt

# Loop over the haplotypes to extract independent haplotype files
(cat ${TEMPDIR}/${RandName}.ListHaplotypes.txt) | while read Haplo
do 
    echo "Processing $Haplo"
    samtools faidx ${TEMPDIR}/${RandName}.full.fa $(cut -f 1 ${TEMPDIR}/${RandName}.full.fa.fai | grep $Haplo | tr '\n' ' ') > ${OUTDIR}/${Haplo}.${RandName}.fasta
done 

# Cleanup
rm ${TEMPDIR}/${RandName}.full.fa
rm ${TEMPDIR}/${RandName}.full.fa.fai
rm ${TEMPDIR}/${RandName}.ListHaplotypes.txt
rmdir ${TEMPDIR}

echo "Processing completed. Haplotypes extracted to ${OUTDIR}"
