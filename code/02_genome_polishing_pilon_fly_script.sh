#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 01:00:00
#SBATCH-J flye_polishing_pilon
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

#Load the modules
module load bioinfo-tools
module load Pilon/1.24
module load bwa/0.7.17
module load samtools/1.19

#Defining the directiories as variables
flye_assembly=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/flye/assembly.fasta
illumina_reads=/home/emdi4663/Genome-Analysis-Paper-I/data/raw_data/Illumina_raw/

work_dir=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly

cd $work_dir 

#Aligning the short paired-end reads (Illumina) to the assembled genome

#Indexing the assembled genome that the reads are to be aligend to
bwa index $flye_assembly

#Performing the alignment between the assembled genome and PacBio reads, and saving as .sam file 
bwa mem $flye_assembly -$illumina_reads/* > aligned_pacbio.sam

#Converting .sam to .bam files 
samtools view aligned_pacbio.sam -b -o aligned_pacbio.bam
samtools sort aligned_pacbio.bam -o aligned_pacbio_sorted.bam

#Indexing the .bam file (needed for Pilon)
samtools index aligned_pacbio_sorted.bam

#Creating Pilon working directory if it already doesn't exist
[ ! -d $work_dir/pilon ] && mkdir $work_dir/pilon

#Running Pilon with the assembled genome as the first argument, followed by the aligned reads to that genome. 
java -Xmx16G -jar /sw/bioinfo/Pilon/1.24/snowy/pilon.jar \
--genome $flye_assembly \
--frags aligned_pacbio_sorted.bam \
--output $work_dir/pilon

