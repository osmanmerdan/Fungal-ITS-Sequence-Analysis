# Mucorales Species Identification
Bioinformatic scripts to analyze ITS (Internal transcribed spacer) sequences to identify Mucorales order fungi. 
Analysis steps:
1. Concatenating fasta files
2. Gathering basic statistics
3. Building BLAST database
4. Parsing BLAST output
5. Creating multiple sequence alignment and phylogenetic ML tree

# ---------------------------- #
# !!!WARNING!!!:
# This pipeline assumes that seqs file names and seq headers are same
# Example; sample file name = ERR789 sample fasta header = >ERR789
# All individual seq files are stored in ./seqs directory
# Required programs:
# seqkit , blast, MAFFT , and iqtree2 in the iqtree directory
# ---------------------------- # 
