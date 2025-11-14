#!/usr/bin/env python3
import argparse

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def load_kmers(kmer_file):
    kmers = []
    with open(kmer_file) as f:
        for line in f:
            kmer = line.strip().split()[0]  # first column
            if kmer:
                kmers.append(kmer)
    return kmers

def parse_fasta(fasta_file):
    sequences = {}
    with open(fasta_file) as f:
        seq_id = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = ''.join(seq_chunks)
                seq_id = line[1:]  # remove '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_id:
            sequences[seq_id] = ''.join(seq_chunks)
    return sequences

def find_kmers_in_fasta(kmer_file, fasta_file, output_file):
    kmers = load_kmers(kmer_file)
    sequences = parse_fasta(fasta_file)

    # Prepare dictionary: kmer -> list of sequence IDs
    kmer_hits = {k: [] for k in kmers}

    for kmer in kmers:
        rev_comp = reverse_complement(kmer)
        for seq_id, seq in sequences.items():
            if kmer in seq or rev_comp in seq:
                kmer_hits[kmer].append(seq_id)

    # Write output
    with open(output_file, 'w') as out:
        for kmer, ids in kmer_hits.items():
            if ids:
                out.write(f"{kmer}\t{','.join(ids)}\n")

def main():
    parser = argparse.ArgumentParser(description="Find k-mers and their reverse complements in FASTA sequences.")
    parser.add_argument("-d", "--kmers", required=True, help="File with k-mers (one per line, first column).")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file with sequences.")
    parser.add_argument("-o", "--output", default="matches.txt", help="Output file (default: matches.txt).")
    args = parser.parse_args()

    find_kmers_in_fasta(args.kmers, args.fasta, args.output)

if __name__ == "__main__":
    main()
