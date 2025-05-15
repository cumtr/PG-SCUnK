#!/bin/bash

### This code extract all the assemblies that compose a pan genome graph (.gfa) in distict fasta files ###
# see PG-SCUnK webpage for details
# Developer : Tristan CUMER - t.cumer.sci[at]gmail.com

# Version : 1.1 - 07/05/2025

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
# # test if the graph is in the good format
# test_gfa1_format() {
#     awk '
#         BEGIN { has_header = 0; is_gfa1 = 0; has_p = 0 }
#         $1 == "H" {
#             has_header = 1;
#             if ($0 ~ /VN:Z:1\.0/) {
#                 is_gfa1 = 1
#             }
#         }
#         $1 == "P" { has_p = 1 }
#         END {
#             if (!has_header) {
#                 print "[FAIL] No header found → cannot determine GFA version. Please check your .gfa file";
#                 exit 3;
#             }
#             if (!is_gfa1) {
#                 print "[FAIL] Header indicates NOT GFA 1.0. consider converting your file into gfa1 using : vg convert -gfW <panGenome>";
#                 exit 1;
#             }
#             if (!has_p) {
#                 print "[FAIL] No P lines found → invalid GFA 1.0 for paths";
#                 exit 2;
#             }
#             print "[OK] GFA 1.0 format detected with P lines";
#             exit 0;
#         }
#     ' $1
# }
# test_gfa1_format ${PG} || exit 1

vg validate test.convert.gfa >${TEMPDIR}/${RandName}.validate.out.txt 2>${TEMPDIR}/${RandName}.validate.err.txt

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
