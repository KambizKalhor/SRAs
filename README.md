# SRAs
here is my personalized bioinformatic workflow starting from short reads from SRA database on CARC (HPC of USC)

# A Practical Guide to my Automated bioinformatic workflow in CARC
![GitHub forks](https://img.shields.io/github/forks/KambizKalhor/my_personalized_bioinformatic_workflow)
![GitHub User's stars](https://img.shields.io/github/stars/KambizKalhor)
![X (formerly Twitter) Follow](https://img.shields.io/twitter/follow/KambizKalhor)

![Website](https://img.shields.io/website?url=https%3A%2F%2Fadsteen.github.io%2F&up_message=Steen%20Lab&up_color=ccd5ae&style=for-the-badge&logo=ocean)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/KambizKalhor/my_personalized_bioinformatic_workflow?style=for-the-badge&color=e9edc9)
![GitHub commit activity](https://img.shields.io/github/commit-activity/t/KambizKalhor/my_personalized_bioinformatic_workflow?style=for-the-badge&color=fefae0)
![GitHub Created At](https://img.shields.io/github/created-at/KambizKalhor/my_personalized_bioinformatic_workflow?style=for-the-badge&color=faedcd)
![GitHub contributors](https://img.shields.io/github/contributors/KambizKalhor/my_personalized_bioinformatic_workflow?style=for-the-badge&color=d4a373)


![Static Badge](https://img.shields.io/badge/Bash-%23D989C3?style=for-the-badge&logo=linux&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/Python-%234BBF9E?style=for-the-badge&logo=python&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/R-%23D9C873?style=for-the-badge&logo=R&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/Git-%23F29D7E?style=for-the-badge&logo=GIT&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/docker-%23F27E7E?style=for-the-badge&logo=docker&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/Obsidian-%23D93D93?style=for-the-badge&logo=Obsidian&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/gitkraken-%2359C1D9?style=for-the-badge&logo=gitkraken&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/LaTeX-%2333A67C?style=for-the-badge&logo=LaTeX&logoColor=black&labelColor=%23FFE4C4)
![Static Badge](https://img.shields.io/badge/mysql-%23F27A5E?style=for-the-badge&logo=mysql&logoColor=black&labelColor=%23FFE4C4)

This repository provides a comprehensive guide and scripts to help myself to keep track of my bioinformatic workflow and organize my scripts.
in this guide, I'm going to explain my code line by line and try to maximize my clean coding skills.

---


# 📑 Table of Contents

- #️⃣ [Batch Job Script Header](#batch-job-script-header)
- 📥 [Inputs and Array Prerequisites](#inputs-and-array-prerequisites)
- 🌐 [PART-One: Download the SRA Files](#part-one-download-the-sra-files)
- 🔨 [PART-Two: Fastq-dump](#part-two-fastq-dump)
- 🔍 [PART-Three: Quality Control](#part-three-quality-control)
- ✂️ [PART-Four: Trimmomatic](#part-four-trimmomatic)
- 🔍 [PART-Five: Second Quality Control](#part-five-second-quality-control)
- 🧬 [PART-Six: Assembly Using Spades](#part-six-assembly-using-spades)
- 📝 [PART-Seven: metaQuast]()

## #️⃣ CARC Batch job script header
---
```
#!/bin/bash
#SBATCH --account=asteen_1130
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nodes=2
#SBATCH --time=01:00:00
#SBATCH --array=1-2
```

## 📥 Inputs and Array Prerequisites
### inputs
in this part you get two arguments, first is the input file(example provided in this repository `job_input_short.txt`) and the output argument is the directory you chosed to save all results
```
# this is the file contains all the SRR names, srr names should be diveded by a newline(enter)
srr_input=$1

# this is the directory to save all the results
output_directory=$2
```

### organize job as array
```
# let's start separate inputs as array
echo "job started"
echo "this is job ${SLURM_ARRAY_TASK_ID}"

# now make a variable to feed as input to main script
line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" $srr_input)
```

## 🌐 PART-One: Download the SRA files
### load needed modules
```
module purge
module load gcc/8.3.0
module load sratoolkit/2.11.0

```

### make a directory for SRA raw file
```
mkdir -p $output_directory/00_raw_read_from_SRA
```

### download srr files
```
prefetch -p ${line} -O $output_directory/00_raw_read_from_SRA
```

## 🔨 PART-Two: Fastq-dump
### make a directory for fastq_dump
```
mkdir -p $output_directory/01_fastq_dump_result
```
### run fastq-dump
```
fastq-dump --split-3 -O $output_directory/01_fastq_dump_result/${line}/ $output_directory/00_raw_read_from_SRA/${line}/*.sra
```

## 🔍 PART-Three: Quality Control
### load needed modules
```
module load gcc/11.3.0
module load fastqc/0.11.9
```

### make a directory to save Quality Control results
```
mkdir -p $output_directory/02_first_quality_control_results/${line}/
```

### execute quality control
```
fastqc -t $SLURM_CPUS_PER_TASK -o  $output_directory/02_first_quality_control_results/${line}/  $output_directory/01_fastq_dump_result/${line}/*.fastq
```

## ✂️ PART-Four: Trimmomatic
### first load the required modules
```
module load gcc/11.3.0
module load trimmomatic/0.39
```

### to see trimmomatic help use
```
trimmomatic -h
```

### make a directory for trim results
```
mkdir -p $output_directory/03_trimmomatic_results/${line}/
```

### run trimmomatic
here in this code the most important variables to change is `LEADING:x`, `TRAILING:x`, `SLIDINGWINDOW:x:xx` and `MINLEN:xx`
you should spend sometime to adjust variables above depending to your input
```
trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 -trimlog $output_directory/03_trimmomatic_results/${line}/trimlog.log $output_directory/01_fastq_dump_result/${line}/*.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_UNPAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_UNPAIRED_R2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:13:20 MINLEN:40
```

## 🔍 PART-Five: Second Quality Control
### load needed modules
```
module load gcc/11.3.0
module load fastqc/0.11.9
```

### make a directory to save Quality Control results
```
mkdir -p $output_directory/04_second_quality_control_results/${line}/
```
### execute quality control
```
fastqc -t $SLURM_CPUS_PER_TASK -o  $output_directory/04_second_quality_control_results/${line}/  $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq
```

## 🧬 PART-Six: Assembely using Spades
### module load
```
module load gcc/8.3.0
module load spades/3.13.0
```

### make a directory to save the assemply results
```
mkdir -p $output_directory/05_Spades_results/${line}
```

## run Spades
here in this code the most important variables to change is `--cov-cutoff auto`
you should spend sometime to adjust variables above depending on input
```
spades.py --meta --threads $SLURM_CPUS_PER_TASK --memory $SLURM_MEM_PER_NODE --pe1-1 $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq --pe1-2 $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq -o $output_directory/05_Spades_results/${line}/
```

## 📝 PART-Seven: metaQuast
### make a directory for results
```
mkdir -p $output_directory/06_metaQuast_results/${line}
```
### run metaQuast
```
quast-5.2.0/metaquast.py $output_directory/05_Spades_results/${line}/scaffolds.fasta -o $output_directory/06_metaQuast_results/${line}/
```



