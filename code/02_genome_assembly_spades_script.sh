#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 03:00:00
#SBATCH-J genome_assembly_spades
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

# Load modules
module load bioinfo-tools spades/3.15.5

# Your command
zcat /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/PacBio_raw/* > /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/PacBio_total.fastq
gzip /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/PacBio_total.fastq

spades.py \
-1 /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
-2 /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
--pacbio /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/PacBio_total.fastq.gz \
-o /home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/spades

rm /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/PacBio_total.fastq.gz
