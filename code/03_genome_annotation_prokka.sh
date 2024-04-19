#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 01:00:00
#SBATCH-J spades_prokka_annotation
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

#Load the modules
module load bioinfo-tools
module load prokka/1.45-5b58020

#The variables 
spades_assembled_genome=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/spades/scaffolds.fasta
output_path=/home/emdi4663/Genome-Analysis-Paper-I/analyses/03_annotation

#The command
prokka \
$spades_assembled_genome \
--outdir $output_path

#Make Prokka .dff output to a readable format
#less -S <XX.dff> output 
