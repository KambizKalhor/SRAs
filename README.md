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
Table of contents
=================

<!--ts-->
   * :hash: [Batch Job Script Header](#hash-batch-job-script-header)
   * 📥 [Inputs and Array Prerequisites](#inputs-and-array-prerequisites)
   * 🌐 [PART-One: Download the SRA Files](#part-one-download-the-sra-files)
   * 🔨 [PART-Two: Fastq-dump](#part-two-fastq-dump)
   * 🔍 [PART-Three: Quality Control](#part-three-quality-control)
   * :scissors: [PART-Four: Trimmomatic](#scissors-part-four-trimmomatic)
   * 🔍 [PART-Five: Second Quality Control](#part-five-second-quality-control)
   * 🧬 [PART-Six: Assembly Using metaSPAdes](#part-six-assembly-using-metaspades)
   * :mega: [PART-Seven: metaQuast](#mega-part-seven-metaquast)
   * :calling: [PART-Eight: bwa](#calling-part-eight-bwa)
   * :wrench: [PART-Nine: samtools](#wrench-part-nine-samtools)
   * :money_with_wings: [PART-Ten: pilon](#money_with_wings-part-ten-pilon)
   * :white_check_mark: [PART-Eleven: CheckM2](#white_check_mark-part-eleven-checkm2)
<!--te-->

Overview
====
![workflow](README_Figures_SRA/workflow_SRA.png)


:hash: [Batch Job Script Header](#batch-job-script-header)
====

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

📥[Inputs and Array Prerequisites](#inputs-and-array-prerequisites)
=====
### inputs
in this part you get two arguments, first is the input file(example provided in this repository [job_input_short.txt](job_input_short.txt)) and the output argument is the directory you chosed to save all results
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

🌐[PART-One: Download the SRA Files](#part-one-download-the-sra-files)
=====
in this step, we download the SRA files using the text file contains a list of SRR numbers.

> [!NOTE]  
> [What is SRA database?](https://en.wikipedia.org/wiki/Sequence_Read_Archive)
>
> [learn more on SRA search results](READMORE/SRA_database.md)



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

🔨[PART-Two: Fastq-dump](#part-two-fastq-dump)
====
> [!TIP]
> The fastq-dump command is a utility from the SRA Toolkit (Sequence Read Archive Toolkit) used to convert SRA files (Sequence Read Archive format) into FASTQ format. The FASTQ format is commonly used for storing sequencing reads with quality scores and is compatible with most downstream bioinformatics tools, unlike the SRA format, which is more specialized.
### make a directory for fastq_dump
```
mkdir -p $output_directory/01_fastq_dump_result
```
### run fastq-dump
```
fastq-dump --split-3 -O $output_directory/01_fastq_dump_result/${line}/ $output_directory/00_raw_read_from_SRA/${line}/*.sra
```

🔍[PART-Three: Quality Control](#part-three-quality-control)
====
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

:scissors: [PART-Four: Trimmomatic](#scissors-part-four-trimmomatic)
====

**Trimmomatic** is a versatile tool used in bioinformatics for **trimming and quality filtering of sequencing reads**. we primarily apply it to **FASTQ files** obtained from last step to improve data quality before downstream analysis. The main tasks include removing low-quality bases, adapters, and unwanted sequences, and filtering out poor-quality reads.

## Key Functions of Trimmomatic
1. **Adapter Removal**: Removes adapters or other unwanted sequences introduced during the library preparation process.
2. **Quality Trimming**: Trims bases from the ends of reads that fall below a certain quality threshold, ensuring that only high-quality sequence data is retained.
3. **Length Filtering**: Discards reads that fall below a specified length after trimming to retain only informative reads.
4. **Sliding Window Trimming**: Trims reads based on a sliding window approach, where the average quality in a defined window must meet a threshold.

## Commonly Used Flags in Trimmomatic

Here's a summary of key flags you can use in Trimmomatic:

| Flag                    | Description |
|-------------------------|-------------|
| **PE** or **SE**        | Specifies whether the input is **paired-end** (PE) or **single-end** (SE) data. Use `PE` for paired-end reads with both forward and reverse files, and `SE` for single-end reads. |
| **ILLUMINACLIP**        | Removes adapter sequences. Takes four parameters: the path to the adapter file, the maximum mismatch count, the palindrome clip threshold, and the simple clip threshold. <br>Example: `ILLUMINACLIP:adapters.fa:2:30:10` |
| **LEADING**             | Trims bases from the start (5' end) of a read if they are below a specified quality threshold. <br>Example: `LEADING:20` |
| **TRAILING**            | Trims bases from the end (3' end) of a read if they are below a specified quality threshold. <br>Example: `TRAILING:20` |
| **SLIDINGWINDOW**       | Performs a sliding window trimming, where a window of a specified size is applied to the read, and if the average quality in that window falls below the set threshold, the read is trimmed at that point. <br>Example: `SLIDINGWINDOW:4:20` (uses a 4-base window with a quality threshold of 20) |
| **MINLEN**              | Discards reads that fall below a specified minimum length after trimming. <br>Example: `MINLEN:36` |
| **CROP**                | Trims reads to a specified length from the start (5' end), useful when you want to retain only the first part of a read. <br>Example: `CROP:100` (keeps only the first 100 bases of each read) |
| **HEADCROP**            | Trims a specified number of bases from the start (5' end) of each read. <br>Example: `HEADCROP:15` (removes the first 15 bases of each read) |
| **AVGQUAL**             | Discards reads with an average quality below a specified threshold, ensuring overall read quality. <br>Example: `AVGQUAL:30` |


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
```
trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 -trimlog $output_directory/03_trimmomatic_results/${line}/trimlog.log $output_directory/01_fastq_dump_result/${line}/*.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_UNPAIRED_R1.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_UNPAIRED_R2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:13:20 MINLEN:40
```

🔍[PART-Five: Second Quality Control](#part-five-second-quality-control)
====
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

🧬[PART-Six: Assembly Using metaSPAdes](#part-six-assembly-using-metaspades)
====
**MetaSPAdes** is an assembly tool designed for metagenomic data, focusing on constructing genomes from mixed community DNA sequencing data.
After running MetaSPAdes, you get several output files that represent different aspects of the assembly. Here’s an overview of the primary output files and their purposes:

| Output File/Folder                | Description |
|-----------------------------------|-------------|
| **contigs.fasta**                 | Contains assembled contigs (continuous sequences) derived from overlapping reads. This file represents the main output and contains sequences that can be analyzed further to identify potential organisms and gene content. |
| **scaffolds.fasta**               | Similar to `contigs.fasta`, but with additional connections between contigs based on paired-end information. Scaffolds are longer sequences that can include gaps where MetaSPAdes inferred connections between contigs. |
| **before_rr.fasta**               | The contigs generated before repeat resolution. Repeat resolution is a step that helps in handling regions with repetitive sequences to improve contig uniqueness and accuracy. |
| **assembly_graph.fastg**          | A graph-based representation of the assembly. Each node represents a sequence (often a contig or part of a contig), and each edge represents connections between sequences based on overlap. This file is useful for visualizing complex assembly paths, especially for highly diverse metagenomic samples. |
| **assembly_graph_with_scaffolds.gfa** | A graph representation in **GFA** format that includes scaffolds. It provides a more comprehensive view of the assembly, with nodes and edges representing sequences and their connections, allowing for complex genome visualization and analysis. |
| **KmerCount.log**                 | Logs the k-mer (subsequence) counts used in the assembly. This file provides details on the different k-mer sizes used and the frequency of each k-mer, helping in quality control analysis by showing how well-covered certain regions are. |
| **params.txt**                    | Contains the parameters used in the MetaSPAdes run, providing documentation of the input settings, which is useful for reproducibility and tracking specific configuration details in the assembly. |
| **spades.log**                    | The primary log file for the MetaSPAdes process. It records the steps and any issues encountered during assembly, which is useful for troubleshooting and tracking the assembly process. |
| **warnings.log**                  | Lists any warnings encountered during the assembly, such as low-coverage areas or potential contaminants. This helps identify potential problems with the assembly quality. |

- **Primary Output**: The `contigs.fasta` file is the main outcome, containing assembled sequences.
- **Supporting Outputs**: Files like `assembly_graph.fastg` and `scaffolds.fasta` provide additional structural and linkage information.
- **Log Files**: Files like `spades.log`, `warnings.log`, and `KmerCount.log` are critical for assessing the quality of the assembly and troubleshooting issues.

These outputs provide a comprehensive set of data to help analyze, evaluate, and understand the structure and quality of the assembled metagenomic sequences. in Cameron's class we used scaffold.fasta as main output to do the downstream analysis.



### module load
```
module load gcc/8.3.0
module load spades/3.13.0
```

### make a directory to save the assemply results
```
mkdir -p $output_directory/05_Spades_results/${line}
```

### run metaSPAdes
here in this code the most srtdarsdrs
you should spend sometime to adjust variables above depending on input
```
spades.py --meta --threads $SLURM_CPUS_PER_TASK --memory $SLURM_MEM_PER_NODE --pe1-1 $output_directory/03_trimmomatic_results/${line}/TRIMMED_FORWARD_PAIRED_R1.fastq --pe1-2 $output_directory/03_trimmomatic_results/${line}/TRIMMED_REVERSE_PAIRED_R2.fastq -o $output_directory/05_Spades_results/${line}/
```

:mega: [PART-Seven: metaQuast](#part-seven-metaQuast)
====
### first you need to download metaquast from the following website [HERE](https://quast.sourceforge.net/metaquast.html)
### you need to unzip it using 
```
tar -zxvf filename.tar.gz
```

### make a directory for results
```
mkdir -p $output_directory/06_metaQuast_results/${line}
```
### run metaQuast
```
metaQuast/quast-5.2.0/metaquast.py $output_directory/05_Spades_results/${line}/scaffolds.fasta -o $output_directory/06_metaQuast_results/${line}/
```


:calling: [PART-Eight: bwa](#part-eight-bwa)
====

### load modules
```
module purge
module load gcc/11.3.0
module load bwa/0.7.17
```
### make a directory for results and copy the scaffolds.fasta as input to that directory
```
mkdir -p $output_directory/07_bwa/${line}
cp $output_directory/05_Spades_results/${line}/scaffolds.fasta $output_directory/07_bwa/${line}/scaffolds.fasta
```

### paths to scaffold and R1 and R2
```
scaffold_fasta_path=$output_directory/07_bwa/${line}/scaffolds.fasta
R1_path=$output_directory/01_fastq_dump_result/${line}/*_1.fastq
R2_path=$output_directory/01_fastq_dump_result/${line}/*_2.fastq
```

### run index (make the genome accessible for alignment)
```
bwa index $scaffold_fasta_path
```

### alignment steps (run this one on original fastq files before trimming-> 01_fastq_dump_result)
```
bwa aln -n 0 -t 16 $scaffold_fasta_path $R1_path > $output_directory/07_bwa/${line}/R1.sai
bwa aln -n 0 -t 16 $scaffold_fasta_path $R2_path > $output_directory/07_bwa/${line}/R2.sai
```

### make path for new outputs
```
R1_sai_path=$output_directory/07_bwa/${line}/R1.sai
R2_sai_path=$output_directory/07_bwa/${line}/R2.sai
# i should read about this outputs, it uses R1 and database to align them
```
### combine them
```
bwa sampe $scaffold_fasta_path $R1_sai_path $R2_sai_path $R1_path $R2_path > $output_directory/07_bwa/${line}/combined_R1_R2.sam
```

:wrench: [PART-Nine: samtools](#part-nine-samtools)
====

### load modules
```
module purge
module load gcc/9.2.0
module load samtools/18.0.4
```

### run sotneareiodhn -> this command will produce a .fai file
```
samtools faidx $scaffold_fasta_path
```
### path to output of the last command
```
fai_output=$output_directory/07_bwa/${line}/*.fai
sam_output=$output_directory/07_bwa/${line}/combined_R1_R2.sam

samtools import $fai_output $sam_output $output_directory/07_bwa/${line}/R1_R2.bam
samtools sort $output_directory/07_bwa/${line}/R1_R2.bam -o $output_directory/07_bwa/${line}/R1_R2.sorted.bam
samtools index $output_directory/07_bwa/${line}/R1_R2.sorted.bam
```

:money_with_wings: [PART-Ten: pilon](#part-ten-pilon)
====

### load modules
```
module purge
module load gcc/8.3.0
module load pilon/1.22
```
### make directory for pilon results
```
mkdir -p $output_directory/08_pilon_results/${line}
```
### make a path for the input
```
sorted_sam_output=$output_directory/07_bwa/${line}/R1_R2.sorted.bam
```
# warning -> bam or sam in command above?
### run pilon
```
pilon --genome $scaffold_fasta_path --frags $sorted_sam_output --fix all,local,breaks --changes --outdir $output_directory/08_pilon_results/${line}  > $output_directory/08_pilon_results/${line}/pilon.log
```


:white_check_mark: [PART-Eleven: CheckM2](#white_check_mark-part-eleven-checkm2)
====

### first you need to install it, the best way to do it on CARC is the following conda environment
```
# make conda environment
conda create -n checkm2 -c conda-forge -c bioconda python=3.8 checkm2

# activate conda environment you just made
conda activate checkm2

# to make sure it is installed
checkm2 -h
```
please do the following steps once and don't use the discovery login node, reminder that you can use the command below to ask for a compute node
```
salloc --ntasks=1 --cpus-per-task=16 --time=01:00:00 --mem=128GB -p main
```
### also if you want to run CheckM2 as bash script, you should activate conda using the following commands
```
module purge
eval "$(conda shell.bash hook)"
conda activate checkm2
```

### you also need to download a database which CheckM2 uses
```
checkm2 database --download --path /path/to/save/database
```


### run CheckM2
```
checkm2 predict -t $SLURM_CPUS_PER_TASK -x input_file_type -i path/to/input.fasta -o path/to/output
# example
checkm2 predict -t 16 -x .fasta -i pilon.fasta -o checkm2_out
```

