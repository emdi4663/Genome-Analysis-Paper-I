#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 05:00:00
#SBATCH-J genome_assembly_flye
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

# Load modules
module load bioinfo-tools Flye/2.9.1

# Your commands
flye --pacbio-raw /home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/PacBio_raw/* \
--out-dir /home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/flye --threads 2
