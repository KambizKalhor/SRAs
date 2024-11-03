#!/bin/bash
#SBATCH --account=asteen_1130
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nodes=2
#SBATCH --time=01:00:00
#SBATCH --array=1-2

#############################
# this is the file contains all the SRR names, srr names should be diveded by a newline(enter)
srr_input=$1
#############################
# this is the directory to save all the results
output_directory=$2
#############################
echo "job started"
echo "this is job ${SLURM_ARRAY_TASK_ID}"
#############################
# now make a variable to feed as input to main script
line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" $srr_input)
#############################







# PART-one: download the SRA files
## load needed modules
module purge
module load gcc/8.3.0
module load sratoolkit/2.11.0

## make a directory for SRA raw file
mkdir -p $output_directory/00_raw_read_from_SRA

## download srr files
prefetch -p ${line} -O $output_directory/00_raw_read_from_SRA
#############################










# PART-two: fastq-dump
## make a directory for fastq_dump
mkdir -p $output_directory/01_fastq_dump_result
fastq-dump --split-3 -O $output_directory/01_fastq_dump_result/${line}/ $output_directory/00_raw_read_from_SRA/${line}/*.sra
#############################









# PART-three: Quality Control
## load needed modules
module load gcc/11.3.0
module load fastqc/0.11.9

## make a directory to save Quality Control results
mkdir -p $output_directory/02_first_quality_control_results/${line}/

## execute quality control
fastqc -t $SLURM_CPUS_PER_TASK -o  $output_directory/02_first_quality_control_results/${line}/  $output_directory/01_fastq_dump_result/${line}/*.fastq
#############################








# PART-Four
## running trimmomatic
### first load the required modules
module load gcc/11.3.0
module load trimmomatic/0.39

### to see trimmomatic help use
trimmomatic -h

### make a directory for trim results
mkdir -p $output_directory/03_trimmomatic_results/${line}/

### run trimmomatic
#### here in this code the most important variables to change is LEADING:x TRAILING:x SLIDINGWINDOW:x:xx MINLEN:xx
#### you should spend sometime to adjust variables above depending to your input
trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 -trimlog $output_directory/03_trimmomatic_results/${line}/trimlog.log $output_directory/01_fastq_dump_result/${line}/*.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_UNPAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_UNPAIRED_R2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:13:20 MINLEN:40
#############################






# PART-Five: Second Quality Control
## load needed modules
module load gcc/11.3.0
module load fastqc/0.11.9

## make a directory to save Quality Control results
mkdir -p $output_directory/04_second_quality_control_results/${line}/

## execute quality control
fastqc -t $SLURM_CPUS_PER_TASK -o  $output_directory/04_second_quality_control_results/${line}/  $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq
#############################





# PART-Six: assembely using Spades
## module load
module load gcc/8.3.0
module load spades/3.13.0

## make a directory to save the assemply results
mkdir -p $output_directory/05_Spades_results/${line}

## run Spades
#### here in this code the most important variables to change is --cov-cutoff auto
#### you should spend sometime to adjust variables above depending on input
spades.py --meta --threads $SLURM_CPUS_PER_TASK --memory $SLURM_MEM_PER_NODE --pe1-1 $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq --pe1-2 $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq -o $output_directory/05_Spades_results/${line}/



#############################
# PART-Seven: metaQuast
# make a directory for results
mkdir -p $output_directory/06_metaQuast_results/${line}
#run metaQuast
metaQuast/quast-5.2.0/metaquast.py $output_directory/05_Spades_results/${line}/scaffolds.fasta -o $output_directory/06_metaQuast_results/${line}/


#############################
# PART-Eight: bwa
# load modules
module purge
module load gcc/11.3.0
module load bwa/0.7.17

# make a directory for results and copy the scaffolds.fasta as input to that directory
mkdir -p $output_directory/07_bwa/${line}
cp $output_directory/05_Spades_results/${line}/scaffolds.fasta $output_directory/07_bwa/${line}/scaffolds.fasta

# paths to scaffold and R1 and R2
scaffold_fasta_path=$output_directory/07_bwa/${line}/scaffolds.fasta
R1_path=$output_directory/01_fastq_dump_result/${line}/*_1.fastq
R2_path=$output_directory/01_fastq_dump_result/${line}/*_2.fastq

# run index (make the genome accessible for alignment)
bwa index $scaffold_fasta_path

# alignment steps (run this one on original fastq files before trimming-> 01_fastq_dump_result)
bwa aln -n 0 -t 16 $scaffold_fasta_path $R1_path > $output_directory/07_bwa/${line}/R1.sai
bwa aln -n 0 -t 16 $scaffold_fasta_path $R2_path > $output_directory/07_bwa/${line}/R2.sai


# make path for new outputs
R1_sai_path=$output_directory/07_bwa/${line}/R1.sai
R2_sai_path=$output_directory/07_bwa/${line}/R2.sai
# i should read about this outputs, it uses R1 and database to align them


# combine them
bwa sampe $scaffold_fasta_path $R1_sai_path $R2_sai_path $R1_path $R2_path > $output_directory/07_bwa/${line}/combined_R1_R2.sam


#############################
# PART-Nine: samtools

# load modules
module purge
module load gcc/9.2.0
module load samtools/18.0.4


# run sotneareiodhn -> this command will produce a .fai file
samtools faidx $scaffold_fasta_path



# path to output of the last command
fai_output=$output_directory/07_bwa/${line}/*.fai
sam_output=$output_directory/07_bwa/${line}/combined_R1_R2.sam


samtools import $fai_output $sam_output $output_directory/07_bwa/${line}/R1_R2.bam
samtools sort $output_directory/07_bwa/${line}/R1_R2.bam -o $output_directory/07_bwa/${line}/R1_R2.sorted.bam
samtools index $output_directory/07_bwa/${line}/R1_R2.sorted.bam

#############################
# PART-Ten: pilon


# load modules
module purge
module load gcc/8.3.0
module load pilon/1.22

# make directory for pilon results
mkdir -p $output_directory/08_pilon_results/${line}

# make a path for the input
sorted_sam_output=$output_directory/07_bwa/${line}/R1_R2.sorted.bam

# run pilon
pilon --genome $scaffold_fasta_path --frags $sorted_sam_output --fix all,local,breaks --changes --outdir $output_directory/08_pilon_results/${line}  > $output_directory/08_pilon_results/${line}/pilon.log




#############################
# PART-Eleven: CheckM2

# first	activate conda environment
# if you haven't install the CheckM2, please find installation on my github
# also you need	to download a database before running this commands
module purge
eval "$(conda shell.bash hook)"
conda activate checkm2

# make a directory for CheckM2 results
mkdir -p $output_directory/09_CheckM2_results/${line}

# run CheckM2
checkm2 predict -t $SLURM_CPUS_PER_TASK -x .fasta -i $output_directory/08_pilon_results/${line}/pilon.fasta -o $output_directory/09_CheckM2_results/${line}/


