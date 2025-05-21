#!/bin/bash

### This code extract all the assemblies that compose a pan genome graph (.gfa) in distict fasta files ###
# see PG-SCUnK webpage for details : https://github.com/cumtr/PG-SCUnK
# Developer : Tristan CUMER - t.cumer.sci[at]gmail.com

# Version : 21/05/2025

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
echo "[INFO] Random ID :" ${RandName}

# Create the temporary directory
mkdir -p ${OUTDIR}
# mkdir -p ${TEMPDIR}/${RandName}

#### Process the Pan Genome ####

echo "[INFO] Testing if the graph is in a valid format"

vg validate ${PG} >${TEMPDIR}/${RandName}.validate.out.txt 2>${TEMPDIR}/${RandName}.validate.err.txt

if ! grep -q "graph: valid" "${TEMPDIR}/${RandName}.validate.err.txt"; then
  echo "       Graph is not valid, please make sure your graph follows gfa formating specs."
  rm ${TEMPDIR}/${RandName}.validate.out.txt 
  rm ${TEMPDIR}/${RandName}.validate.err.txt
  exit 1
fi

if grep -q "graph: valid" "${TEMPDIR}/${RandName}.validate.err.txt"; then
  echo "       Graph is valid."
  rm ${TEMPDIR}/${RandName}.validate.out.txt 
  rm ${TEMPDIR}/${RandName}.validate.err.txt
fi

echo "[INFO] Extracting the multifasta from the graph"

# Extract the fasta with all the assemblies & index it
vg paths --extract-fasta -x ${PG} > ${TEMPDIR}/${RandName}.full.fa

# Index the panGenome fasta using samtools
samtools faidx ${TEMPDIR}/${RandName}.full.fa

#### Process the multifasta ####
echo "[INFO] Splitting of the mulifasta per haplotype"

# Identify Haplotype names assuming PanSN formatting
cut -f 1 ${TEMPDIR}/${RandName}.full.fa.fai | awk -F"#" '{print $1"#"$2}' | sort | uniq > ${TEMPDIR}/${RandName}.ListHaplotypes.txt

# Loop over the haplotypes to extract independent haplotype files
(cat ${TEMPDIR}/${RandName}.ListHaplotypes.txt) | while read Haplo
do 
    echo "       Processing $Haplo"
    samtools faidx ${TEMPDIR}/${RandName}.full.fa $(cut -f 1 ${TEMPDIR}/${RandName}.full.fa.fai | grep $Haplo | tr '\n' ' ') > ${OUTDIR}/${Haplo}.${RandName}.fasta
done 

# Cleanup
rm ${TEMPDIR}/${RandName}.full.fa
rm ${TEMPDIR}/${RandName}.full.fa.fai
rm ${TEMPDIR}/${RandName}.ListHaplotypes.txt

echo "[INFO] Processing completed. Haplotypes extracted to: ${OUTDIR}/"
