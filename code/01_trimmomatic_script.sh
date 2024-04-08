#!/bin/bash-l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 00:00:00
#SBATCH-J trimmomatic_preprocessing 
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

# Load modules
module load bioinfo-tools
module load trimmomatic/0.39

# Your commands
java -jar /sw/bioinfo/trimmomatic/0.39/snowy/trimmomatic-0.39.jar PE
-threads 2 
-phred33
/home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
/home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz
/home/emdi4663/Genome-Analysis-Paper-I/analyses/01_preprocessing/fastqc_trimmed/trimmed_P_1
/home/emdi4663/Genome-Analysis-Paper-I/analyses/01_preprocessing/fastqc_trimmed/trimmed_U_1
/home/emdi4663/Genome-Analysis-Paper-I/analyses/01_preprocessing/fastqc_trimmed/trimmed_P_2
/home/emdi4663/Genome-Analysis-Paper-I/analyses/01_preprocessing/fastqc_trimmed/trimmed_U_2
ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.39/snowy/adapters/TruSeq3-PE.fa:2:30:10 
LEADING:3 
TRAILING:3 
SLIDINGWINDOW:4:15

