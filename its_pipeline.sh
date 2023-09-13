#!/bin/bash
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Author: Osman Merdan, MD
# Uludag University Medical School Department of Medical Microbiology 
# Last Updated : 6 July 2023
# Bioinformatic analysis pipeline for the project:
# Retrospective determination of species distribution and antifungal susceptibility of Mucorales order fungi isolated from clinical specimens.

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

# Stop on any error 
# -e: Exit immediately if a command exits with a non-zero status.
# -u  Treat unset variables as an error when substituting.
# -x  Print commands and their arguments as they are executed.
set -ue

# 1- Concetanete All Individual Fasta Files
# -------------------------------------------------------------------------------------------------------------------------
cat ./seqs/*.fasta >> ./seqs.fasta 

# 2- Remove Any Blank Lines or Spaces From Multi Fasta File
# -------------------------------------------------------------------------------------------------------------------------
seqkit seq --threads 6 --line-width 50 seqs.fasta > tmpfile && mv tmpfile seqs.fasta


# 3- Get Summary Stats Of The Concatenated Fasta File
# -------------------------------------------------------------------------------------------------------------------------
 seqkit stats seqs.fasta  --tabular --all > stats.tsv

# 4- Get Stats About Each Fasta Sequence In Multi-Fasta File 
# -------------------------------------------------------------------------------------------------------------------------
 seqkit fx2tab\
 seqs.fasta\
 --base-content A\
 --base-content T\
 --base-content G\
 --base-content C\
 --base-content N\
 --base-count A\
 --base-count T\
 --base-count G\
 --base-count C\
 --base-count N\
 --gc\
 --gc-skew\
 --length\
 --header-line\
 --name\
 --threads 6\
 > seqstats.tsv

# Quality check step
echo "**** Quality Control ****"
cat stats.tsv
echo "========================="


# 5- Creating Blast Database for Fungal ITS Sequences
# -------------------------------------------------------------------------------------------------------------------------
read -rp "Do you want to build or update the blast database (y/n): " choice
if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
    echo "Building/Updating blast database"

# Download latest refseq Fungal ITS Sequences
    wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz
    gunzip fungi.ITS.fna.gz
    mkdir -p blastdb
    makeblastdb\
     -input_type fasta\
     -dbtype nucl\
     -parse_seqids\
     -in fungi.ITS.fna\
     -out blastdb/ITS_Database_2023
# To get info about database
    echo "**** Database Info ****" 
    blastdbcmd -db ./blastdb/ITS_Database_2023 -info
    echo "========================="
else
    echo "Continuing..."
fi
read -rp "Do you want to continue blast analysis(y/n): " choice
if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
    echo "Continuing..."
else
    echo "Exiting."
    exit 0
fi
# ~~~Alternative Method~~~~
#mkdir -p blastdb
# The blast package includes a script that can be used to download ready made databases.
# List all available blast databases: 
#  update_blastdb.pl --showall 
# Download ITS_RefSeq_Fungi database 
#cd blastdb || exit
#update_blastdb.pl ITS_RefSeq_Fungi --decompress
# Download the taxonomy database.
#update_blastdb.pl taxdb --decompress
# Set BLASTDB variable to enable blast to understand database names
#cd .. 
#export BLASTDB=$BLASTDB:./blastdb
# With this method it is possible to get scientific names of blast hits with output format 6. 

# 6- Running BLAST
# -------------------------------------------------------------------------------------------------------------------------
# Running BLAST with blast archive output
echo "***** BLAST search started at $(date) *****"
mkdir -p blast-results
blastn\
 -task megablast\
 -db blastdb/ITS_Database_2023\
 -query seqs.fasta\
 -out blast-results/blast_result_archive.asn\
 -outfmt 11\
 -max_target_seqs 10\
 -num_threads 6
# Reformatting as standart alignment blast output
blast_formatter\
 -archive blast-results/blast_result_archive.asn\
 -outfmt 0\
 -out blast-results/blast_result_standart.txt
# Reformatting as flat query anchored output with seq identifiers
blast_formatter\
 -archive blast-results/blast_result_archive.asn\
 -outfmt 3\
 -out blast-results/blast_result_flat_qanchored.txt
# Tabular output for downstream analysis
blast_formatter\
 -archive blast-results/blast_result_archive.asn\
 -outfmt "6 qacc saccver stitle qlen slen length qstart qend sstart send qcovs pident nident mismatch gaps bitscore evalue"\
 -max_target_seqs 5\
 -out blast-results/blast_result_tabular.tsv

# 7- Control of Blast Output 
# -------------------------------------------------------------------------------------------------------------------------
echo "******* Number of hits per sample *******"
< blast-results/blast_result_tabular.tsv cut -f 1 | sort | uniq  -c 
echo "******** Total number of samples with at least one hit *********"
< blast-results/blast_result_tabular.tsv cut -f 1 | sort | uniq  -c | wc -l

# 8- Filtering Blast Output Based on Predefined Filtering Criteria
# -------------------------------------------------------------------------------------------------------------------------
# Using blast result parser python script to filter out blast results. 
# This scrpit outputs downstream friendly tabuler format (21 columns)
echo "***** BLAST tabular output filtering started at $(date) *****"
python blast_result_parser.py blast-results/blast_result_tabular.tsv
echo "***** BLAST tabular output filtering finished at $(date) *****"

read -rp "Do you want to create MSA for genus groups and create phylogenetic tree(y/n): " choice
if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
    echo "Continuing..."
    
else
    echo "Exiting."
    exit 0
fi


# 9- Extracting Ref Sequence of Top Unique Hits to Use in Phylogenetic Tree
# -------------------------------------------------------------------------------------------------------------------------
 < blast-results/parsed_blast_result.tsv tail -n +2 | cut -f 2 | sort | uniq > blast-results/unique_top_hits_acc.txt
echo "***** Extracting refrence sequences of top unique hits at $(date) *****"
mkdir -p refs
while read -r i;
do
 blastdbcmd -db blastdb/ITS_Database_2023 -entry "$i" -outfmt "%f" > ./refs/"$i"
done < blast-results/unique_top_hits_acc.txt
echo "***** Reference sequences of top unique hits were extacted at $(date) *****"
rm -f blast-results/unique_top_hits_acc.txt


# 10- Creating Individual MultiFasta File For Each Genera Name
# -------------------------------------------------------------------------------------------------------------------------
mkdir -p genus-groups

# Extracting the unique genera names from tabular filtered blast result file
< ./blast-results/parsed_blast_result.tsv tail -n +2 | cut -f 20 | sort | uniq\
 > blast-results/unique_genus_names.txt
# Removing any unidentified hit 
< ./blast-results/unique_genus_names.txt grep -v "unidentified" \
 > tmpfile && mv tmpfile "./blast-results/unique_genus_names.txt"
# First awk :
#  Finding sample names and sample paths matching each unique genera name
# Second awk:
#  Finding ref names and ref paths matching each unique genera name and
# Concatenating sample and ref sequences 
mkdir -p genus-groups
while read -r i;
do 
 awk -F '\t'\
 -v x="$i"\
 '$20==x {print "./seqs/"$1".fasta"}'\
 ./blast-results/parsed_blast_result.tsv | sort | uniq > temp2

 #'{split($3, words, " "); if (words[1] == x) print "./seqs/"$1".fasta"}'\

 awk -F '\t'\
 -v x="$i"\
 '$20==x {print "./refs/"$2}'\
 ./blast-results/parsed_blast_result.tsv | sort | uniq >> temp2

#'{split($3, words, " "); if (words[1] == x) print "./refs/"$2}'\

 xargs cat < temp2 > ./genus-groups/"$i"
done < blast-results/unique_genus_names.txt

rm -f temp2
rm -f blast-results/unique_genus_names.txt 


# 11- Creating MSA With MAFFT for Genera Groups to Use in Phylogenetic Tree
# -------------------------------------------------------------------------------------------------------------------------
# Ref file headers are long and causes problems during visualizations.
# To get accession numbers and also species names in headers use :
#  seqkit replace --pattern "^(\S+).+?(\w+ \w+).*$" --replacement "\$1-\$2" input > output
# To get only NCBI accessions and sample names in headers use:
#  seqkit seq --threads 6 --only-id input > output
# mafft:
#  Returns MSA format fasta file 
#  --auto: selects best algorthm for the data and computer resources
#  --clustalout: for clustal format output
mkdir -p alignments
echo "****** MSA creation with MAFFT started at $(date) ******"
sleep 5
for i in ./genus-groups/*
do
seqkit seq --threads 6 --only-id "$i" > tmpfile && mv tmpfile "$i"
mafft --auto --reorder --adjustdirectionaccurately --thread 6 "$i"> ./alignments/"$(basename "$i")"
done
echo "****** MSA creation with MAFFT finished at $(date) ******"
sleep 5

rm -rf genus-groups
rm -rf refs


# 12- Creating Genera Based Phylogenetic Trees with IQTree-2
# -------------------------------------------------------------------------------------------------------------------------
# -m MPF: automatic model selection by ModelFinder Plus
# -B 1000: assessing branch support with ultrafast bootstrap approximation
# --name to print full name
echo "****** Phylogenetic tree construction with IQTree2 started at $(date) ******"
sleep 5
mkdir -p ./trees
for i in ./alignments/*
    do
    echo "***** Processing $i Alignment File *****"
    if iqtree/iqtree2 -s "$i" -T 6 -B 1000 -m MFP --prefix ./trees/"$(basename "$i")"; then
        echo "***** $i Tree File Created *****"
        sleep 10
    else
        echo "***** An error occured while processing $i file *****"
        sleep 10
        echo "**** Skipping to next file ****"
        sleep 5
    fi
done
echo "****** Phylogenetic tree construction with IQTree2 finished at $(date) ******"
sleep 5
echo "~~~~~~~~~~ Analysis Complete ~~~~~~~~~~"


