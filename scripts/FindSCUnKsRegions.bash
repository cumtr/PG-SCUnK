#!/bin/bash

## Parse command-line options ##
usage() {
    echo "Usage: $0 -b <basename> -r <reference> -t <tempDir> -o <outDir> -@ <threads>"
    exit 1
}

while getopts ":b:r:t:o:@:" opt; do
  case $opt in
    b) BASENAME="$OPTARG";;
    r) REFGENOME="$OPTARG";;
    t) TEMPDIR="$OPTARG";;
    o) OUTDIR="$OPTARG";;
    @) THREADS="$OPTARG";;
    *) usage ;;
  esac
done

## Check required params ##
if [ -z "$BASENAME" ] || [ -z "$REFGENOME" ] || [ -z "$TEMPDIR" ] || [ -z "$OUTDIR" ] || [ -z "$THREADS" ]; then
    usage
fi

mkdir -p "$TEMPDIR"
mkdir -p "$OUTDIR"

BASE="$(basename "$BASENAME")"
BASE_REFGENOME="$(basename "$REFGENOME")"

echo ""
echo "Input File: $BASENAME"
echo "Output Directory: $OUTDIR"
# if [[ $verbose -eq 1 ]]; then
#   echo "Verbose mode is ON"
# fi

RandName=$(xxd -u -l 3 -p /dev/urandom)
echo ""
echo "Internal random ID :" ${RandName}
echo ""

## Check if reference is indexed ##
if [ ! -f "${REFGENOME}.bwt" ]; then
    echo "[INFO] Reference index not found. Indexing with bwa..."
    bwa index "$REFGENOME"
else
    echo "[INFO] Reference index found."
fi

## List of k-mer types ##
TYPES=("unique" "duplicated" "split")

## Process each k-mer file ##
for TYPE in "${TYPES[@]}"; do
    KMER_FILE="${BASENAME}.${TYPE}.txt"
    
    if [ ! -f "$KMER_FILE" ]; then
        echo "Warning: $KMER_FILE does not exist. Skipping."
        continue
    fi
    echo "       "
    echo "[INFO] Processing $KMER_FILE"

    # Filenames for intermediate & output
    FASTQ_FILE="${TEMPDIR}/${BASE}.${TYPE}.fastq"
    SAM_FILE="${TEMPDIR}/${BASE}.${TYPE}.sam"
    BAM_FILE="${TEMPDIR}/${BASE}.${TYPE}.sorted.bam"
    BED_FILE="${TEMPDIR}/${BASE}.${TYPE}.bed"
    MERGED_BED_FILE="${OUTDIR}/${BASE}.SCUnKs.${TYPE}.bed"

    #  Convert kmers to FASTQ 
    echo "       Converting k-mers to FASTQ for $TYPE SCUnKs"
    awk '{qual=""; for(i=1;i<=length($0);i++) qual=qual "I"; printf("@seq%d\n%s\n+\n%s\n", NR, $0, qual)}' "$KMER_FILE" > "$FASTQ_FILE"

    #  Map with bwa 
    echo "       Mapping $TYPE SCUnKs"
    bwa mem -t "$THREADS" "$REFGENOME" "$FASTQ_FILE" > "$SAM_FILE" 2>/dev/null

    #  Convert SAM to sorted BAM 
    echo "       Converting and sorting BAM for $TYPE SCUnKs"
    samtools view -@ "$THREADS" -bS "$SAM_FILE" | samtools sort -@ "$THREADS" -o "$BAM_FILE"
    samtools index "$BAM_FILE"

    #  Extract mapped regions 
    echo "       Extracting mapped regions for $TYPE SCUnKs"
    bedtools bamtobed -i "$BAM_FILE" > "$BED_FILE"
    bedtools merge -i "$BED_FILE" > "$MERGED_BED_FILE"

    rm "$FASTQ_FILE" "$SAM_FILE" "$BAM_FILE"

    echo "       Finished processing $TYPE."
    echo "       "
done

echo "       "
echo "[INFO] Generating output files"

LEN=$(grep -v ">" "$REFGENOME" | tr -d '\n' | wc -c)

R --vanilla <<EOF > ${OUTDIR}/${BASE}.SCUnKs.log

# loading the bed files
TabUniq = read.table(paste0("${OUTDIR}","/","${BASE}",".SCUnKs.unique.bed"), comment.char = "@", sep = "\t")
TabDup = read.table(paste0("${OUTDIR}","/","${BASE}",".SCUnKs.duplicated.bed"), comment.char = "@", sep = "\t")
TabSplit = read.table(paste0("${OUTDIR}","/","${BASE}",".SCUnKs.split.bed"), comment.char = "@", sep = "\t")

# computing proportions of the gneome in the different categories
PartUniq = (sum(TabUniq[,3]-TabUniq[,2])/${LEN})*100
PartDup = (sum(TabDup[,3]-TabDup[,2])/${LEN})*100
PartSplit = (sum(TabSplit[,3]-TabSplit[,2])/${LEN})*100

PartSCUnKs = sum(PartUniq, PartDup, PartSplit)

# Plotting the repartition of the SCUnKs along the genome
png(paste0("${OUTDIR}/${BASE}",".SCUnKs.position.png") , width = 12, height = 3, units = "in", res = 300)

par(mar=c(5,5,3,2))

plot(0,0, xlim = c(1, ${LEN}), ylim = c(0,0.7), type ="n", 
     axes = F, xlab = "Position along the genome", ylab ="",
     main = "")
mtext( paste0("${BASE_REFGENOME}", " - ",round(PartSCUnKs,digits = 2)," % SCUnKs [",round(PartUniq,digits = 2),"U / ",round(PartDup,digits = 2),"D / ",round(PartSplit,digits = 2),"S]"),
       line = 0.5, cex = 1.1)
axis(1, col = "grey70")
axis(2, at = c(0.1, 0.35, 0.6), las = 2, labels = c("Split", "Duplicated", "Unique"), lwd = 0, line = -1.1)

rect(1, 0, ${LEN}, 0.2, lwd = 1.5, border = NA, col = "grey95")
rect(1, 0.25, ${LEN}, 0.45, lwd = 1.5, border = NA, col = "grey95")
rect(1, 0.50, ${LEN}, 0.7, lwd = 1.5, border = NA, col = "grey95")

for(i in 1:nrow(TabUniq)){
  rect(TabUniq[i,2],0.5, TabUniq[i,3], 0.7, border = "grey40", col = "grey40", lwd=0)
}
for(i in 1:nrow(TabDup)){
  rect(TabDup[i,2],0.25, TabDup[i,3], 0.45, border = "darkorange", col = "darkorange", lwd=0)
}
for(i in 1:nrow(TabSplit)){
  rect(TabSplit[i,2],0, TabSplit[i,3], 0.2, border = "darkblue", col = "darkblue", lwd=0)
}
rect(1, 0, ${LEN}, 0.2, lwd = 1.5, border = "grey70", col = NA)
rect(1, 0.25, ${LEN}, 0.45, lwd = 1.5, border = "grey70", col = NA)
rect(1, 0.50, ${LEN}, 0.7, lwd = 1.5, border = "grey70", col = NA)

dev.off()

# Writing the SCUnKs in proportion of the reference genome
SCUnKs = t(c(round(PartSCUnKs,4),
           round(PartUniq,4),
           round(PartDup,4),
           round(PartSplit,4)))
colnames(SCUnKs) = c("#%_GenomeInSCUnKs", "%_UniqueSCUnKs", "%_DuplicatedSCUnKs", "%_SplitSCUnKs")
write.table(SCUnKs, paste0("${OUTDIR}/${BASE}",".SCUnKs.proportion.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

EOF

echo "       "
echo "[INFO] Finished."
