#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 00:30:00
#SBATCH-J spades_blast_visualisation
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

#Load the modules
module load bioinfo-tools
module load blast/2.15.0+ 
module load artemis/16.0.0


#The variables
spades_query=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/spades/scaffolds.fasta
genome_ref=/home/emdi4663/Genome-Analysis-Paper-I/data/genome_reference/*.fna.gz

output_path=/home/emdi4663/Genome-Analysis-Paper-I/code

#The alignment command
blastn -query $spades_query -subject  <(gunzip -c $genome_ref) -outfmt 6 -out $output_path/spades_blast.out

##Alignment with first 3 contigs (edited assmebly)
edited_spades=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/spades/spades_filtered_conc.fasta 
blastn -query $edited_spades -subject  <(gunzip -c $genome_ref) -outfmt 6 -out $output_path/spades_conc_blast.out