#!/bin/bash


# set -e # exit on error
start=`date +%s`

###########################################
####           Initialisation          ####
###########################################

#### Initialize variables ####
input_assemblies_directory=""
temp_directory=""
output_base=""
kmer_size=100
# threads=1
verbose=0

# Function to display help message
usage() {
    echo "Usage: $0 -a <path/to/assemblies/> -o <outputDir/outputBasename> -t <WorkDir> (-k <kmer_size> -v)"
    exit 1
}

# Use getopts to parse the flags
while getopts "a:o:t:k:v" flag; do
  case "${flag}" in
    a) input_assemblies_directory=${OPTARG} ;;
    o) output_base=${OPTARG} ;;   
    t) temp_directory=${OPTARG} ;;     
    k) kmer_size=${OPTARG} ;;
    v) verbose=1 ;;           # -v flag for verbose mode (no argument)
    *) usage ;; 
  esac
done

if [[ $verbose -eq 1 ]]; then
    set -x
    echo "Verbose mode enabled: 'set -x' is active, all commands will be printend to outputDir/outputBasename.log file)"
fi


#### Display the parsed options ####
echo ""
echo "Output Directory: $output_base"
echo "kmer Size: $kmer_size"
if [[ $verbose -eq 1 ]]; then
  echo "Verbose mode is ON"
fi

RandName=$(xxd -u -l 3 -p /dev/urandom)
echo ""
echo "Internal random ID :" ${RandName}
echo ""

#### Create necessary directories ####
mkdir -p ${temp_directory}/${RandName}
output_dir="$(dirname "$output_base")"
mkdir -p "$output_dir" "$temp_directory"

#### Remove output file if already existing ####
[ -e $output_base.Sharedkmers.log ] && rm $output_base.Sharedkmers.log
[ -e $output_base.Sharedkmers.txt ] && rm $output_base.Sharedkmers.txt
[ -e $output_base.Sharedkmers.Hist.txt ] && rm $output_base.Sharedkmers.Hist.txt

#### Define usefull variables ####
if [ "$verbose" == 0 ]; then verboseKMC="-hp"; fi
if [ "$verbose" == 0 ]; then verboseKMCtools="-hp -v"; fi

################################################################
#####          Identify uniq kmers in each assembly         ####
#####              and merge across assemblies              ####
################################################################

start=$(date +%s)
echo "[$((`date +%s`-start)) sec] Starting k-mer extraction and intersection ..."

# Initialize variables
first_file=true
intersect_file=""

# Process each assembly file iteratively
for FILE in ${input_assemblies_directory}/*.fasta; do
  NAME=$(basename "${FILE%.*}")
  echo "[INFO] Processing file: $NAME"

  # Extract kmers for the current file
  temp_kmer_db="${temp_directory}/${RandName}/$NAME"
  kmc -t1 -cm -k"$kmer_size" -ci0 -cx1 -fm ${verboseKMC:+"$verboseKMC"} "$FILE" "$temp_kmer_db" "${temp_directory}/${RandName}/" > $output_base.Sharedkmers.log out 2>&1
  # set kmer counter to 1 (mandatory for later counting)
  temp_kmer_db_1="${temp_directory}/${RandName}/$NAME_1"
  kmc_tools transform "$temp_kmer_db" set_counts 1 "$temp_kmer_db_1"

  if $first_file; then
    # First file: copy kmer database to initialize the union file
    mv "$temp_kmer_db_1.kmc_suf" "${temp_directory}/${RandName}/union_initial.kmc_suf"
    mv "$temp_kmer_db_1.kmc_pre" "${temp_directory}/${RandName}/union_initial.kmc_pre"
    union_file="${temp_directory}/${RandName}/union_initial"
    first_file=false
  else
    # Perform union with the current kmer set
    next_union_file="${temp_directory}/${RandName}/union_temp_$NAME"
    kmc_tools -t1 ${verboseKMCtools:+"$verboseKMCtools"} simple "$union_file" "$temp_kmer_db_1" union "$next_union_file" -ocsum

    # Clean up previous union files after generating the new one
    rm "$temp_kmer_db.kmc_suf" "$temp_kmer_db.kmc_pre"

    # Update union_file to the new union
    previous_union_file="$union_file"
    union_file="$next_union_file"
  fi

  # Clean up kmc files for the current assembly
  [ -e "$previous_union_file.kmc_suf" ] && rm "$previous_union_file.kmc_suf" "$previous_union_file.kmc_pre"
done

# Finalize the results
final_kmer_file="${temp_directory}/${RandName}/${RandName}_allUniqKmers"
mv "$union_file.kmc_suf" "$final_kmer_file.kmc_suf"
mv "$union_file.kmc_pre" "$final_kmer_file.kmc_pre"

kmc_tools -t1 ${verboseKMCtools:+"$verboseKMCtools"} transform "$final_kmer_file" dump "$output_base.Sharedkmers.txt"
kmc_tools -t1 ${verboseKMCtools:+"$verboseKMCtools"} transform "$final_kmer_file" histogram "$output_base.Sharedkmers.hist.txt"

# Clean up final temporary intersection file
# rm "$final_kmer_file.kmc_suf" "$final_kmer_file.kmc_pre"

echo "[$((`date +%s`-start)) sec] k-mer extraction and intersection: DONE"

# Cleanup: remove unnecessary temporary files
rm -rf "${temp_directory}/${RandName}"


