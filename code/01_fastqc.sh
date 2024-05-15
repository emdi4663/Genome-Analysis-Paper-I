#!/bin/bash-l

#Modules needed
module load bioinfo-tools
module load FastQC/0.11.9

#Creating links to Illumina files
ln -s /proj/uppmax2024-2-7/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina/* ~/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw

cd ~/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw

#running FASTQC on Illumina raw data
forward_read=E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
reverse_read=E745-1.L500_SZAXPI015146-56_2_clean.fq.gz
fastqc $forward_read $reverse_read -o /home/emdi4663/Genome-Analysis-Paper-I/analyses/01_preprocessing/fastqc_raw

