![Mushrooms](https://github.com/osmanmerdan/Fungal-ITS-Sequence-Analysis/assets/91588759/cc61b58d-7b5b-4ff0-ada0-5694c73331c8)

# Fungal Species Identification with ITS Sequencing
Bioinformatic scripts to analyze ITS (Internal transcribed spacer) sequences to identify Mucorales order fungi. 

**Files:**
- _its_pipeline.sh_: main bash script to analyze FASTA sequences
- _blast_result_parser.py_: Python 3 script to filter out BLAST hits

**Analysis steps:**

![Untitled-2](https://github.com/osmanmerdan/Fungal-ITS-Sequence-Analysis/assets/91588759/c6114cdf-8b3a-4132-9a28-b19b5950a791)


**Notes:**

- This pipeline assumes that seqs file names and seq headers are the same.
- For example, Sample file name: ERR789.fasta, Sample fasta header: >ERR789
- All individual seq files are stored in the "./seqs" directory

**Used Command Line Programs:**
- Seqkit, Blast, MAFFT, and Iqtree2
