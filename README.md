# Mucorales Species Identification
Bioinformatic scripts to analyze ITS (Internal transcribed spacer) sequences to identify Mucorales order fungi. 

**Files:**
- mucorales-turkey.sh : main bash script to analyze FASTA sequences
- blast_result_parser.py: Python 3 script to filter out BLAST hits

**Analysis steps:**
1. Concatenating fasta files
2. Gathering basic statistics
3. Building BLAST database
4. Parsing BLAST output
5. Creating multiple sequence alignment and phylogenetic ML tree

_!!!WARNING!!!_

This pipeline assumes that seqs file names and seq headers are same.

For example; sample file name = ERR789.fasta sample fasta header = >ERR789

All individual seq files are stored in the "./seqs" directory

_Required programs_:

seqkit , blast, MAFFT , and iqtree2 in the iqtree directory
