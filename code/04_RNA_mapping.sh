#!/bin/bash -l

#SBATCH-A uppmax2024-2-7
#SBATCH-M snowy
#SBATCH-p core
#SBATCH-n 2
#SBATCH-t 04:00:00
#SBATCH-J RNA_mapping
#SBATCH--mail-type=ALL
#SBATCH--mail-user emdi4663@student.uu.se
#SBATCH--output=%x.%j.out

#Load the modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.19

#Defining the directiories as variables
spades_assembly=/home/emdi4663/Genome-Analysis-Paper-I/analyses/02_genome_assembly/spades/scaffolds.fasta
RNA_seq_reads=/home/emdi4663/Genome-Analysis-Paper-I/data/RNA_seq

work_dir=/home/emdi4663/Genome-Analysis-Paper-I/analyses/04_diff_expression_analysis

cd $work_dir 

#Indexing the assembled genome that the reads are to be aligend to
bwa index $spades_assembly

#Aligning the RNA reads to the assembled genome. 2 conditions, 3 pairs of reads in each.
for sample in "69" "70" "71" ; do
        forward_read="$RNA_seq_reads/trim_paired_ERR17979${sample}_pass_1.fastq.gz"
        reverse_read="$RNA_seq_reads/trim_paired_ERR17979${sample}_pass_2.fastq.gz"
        bwa mem $spades_assembly $forward_read $reverse_read > serum_${sample}_mapped.sam
done

for sample in "72" "73" "74" ; do
        forward_read="$RNA_seq_reads/trim_paired_ERR17979${sample}_pass_1.fastq.gz"
        reverse_read="$RNA_seq_reads/trim_paired_ERR17979${sample}_pass_2.fastq.gz"
        bwa mem $spades_assembly $forward_read $reverse_read > bh_${sample}_mapped.sam
done

#Converting the .sam files to .bam, followed by sorting of the .bam files 
for samfile in $work_dir/*_mapped.sam; do
	samtools view ${samfile} -b -o $(basename -s ".sam" ${samfile}).bam       
	bamfile="$(basename -s ".sam" ${samfile}).bam"
        rm ${samfile}
	samtools sort ${bamfile} -o $(basename -s ".bam" ${bamfile})_sorted.bam
	rm ${bamfile}	
done

#Have to save ti to /proj/uppmax2024-2-7/private