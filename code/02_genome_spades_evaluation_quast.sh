#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 01:00:00
#SBATCH-J spades_evaluation_quast
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

#Load the modules
module load bioinfo-tools
module load  quast/5.0.2

#The variables
spades_assembly=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/spades/scaffolds.fasta
genome_ref=/home/emdi4663/Genome-Analysis-Paper-I/data/genome_reference/*.fna.gz
output_path=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly


#The commands
python /sw/bioinfo/quast/5.0.2/snowy/bin/quast.py \
$spades_assembly \
-r $genome_ref \
-o $output_path
