# Nanopore assembly and polishing for raw assembly
Nan Hu, July 11, 2022

---

## Software preparation
1. flye
2. porechop
3. trimmomatic
4. nextpolish

## Trim nanopore reads using Porechop
Raw nanopore reads contain sequencing adaptors. Before we use our reads to assemble our genome, we need to use software to trim the adaptors from the reads. Here, we get Porechop to fulfil this demand. The software is already installed so you do not need to install again if you are in HPCC of TTU.

It will be easier to run the scripts when we merge all the reads together. You may go into the directory containing all reads and run a simple `cat` command:
```bash
cat *.fastq > [SpeciesName].merged.fastq
```

Then, we can run the `porechop-runner.py` to trim the reads.

```bash
#!/bin/bash
#SBATCH -J trim
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128

date

python /home/gufeng/Porechop/porechop-runner.py -i ./[SpeciesName].merged.fastq -o ./[SpeciesName].merged.trimmed.fastq --discard_middle --threads 128

```
Before running on any assemblers, we need to have this trimmed `.fastq` file.

## Assemble nanopore reads using Flye
Flye is a one-line command software to assemble nanopore reads. We basically use the defalut settings in assemble our genome. It requires a estimated genome size, which may fetched from closed related species genome.

```bash
#!/bin/bash
#SBATCH -J flye
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128

date

/home/nhu/software/Flye/bin/flye --nano-raw ./[SpeciesName].merged.trimmed.fastq --genome-size 350m --out-dir out_nano --threads 128 
```
> The new version of Flye is somehow unstable in our HPCC system due to RAM issue. Thus, it may abort in the middle of run and report an out of memory error. However, good news is that Flye is capable with resuming from where it aborted. We can simply add `--resume` to the command line above to continue the run.

When the run finishes, the raw assembly will be in `./out_nano/assembly.fasta`. We need to use this file as the input file for polishing.

## Polish assembly using Nextpolish
Nanopore reads have longer length than NGS data but poor accuracy. For polishing raw nanopore assembly, we need to use NGS illumina reads to correct the nanopore assembly. Similarly, the illumina reads needs to be trimmed before using. We will use Trimmomatic to trim the illumina reads.

```bash
#!/bin/bash
#SBATCH -J trim
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 16

java -jar /home/nhu/software/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE -phred33 \
        [SpeciesName].R1.fastq \
        [SpeciesName].R2.fastq \
        [SpeciesName].R1_paired.fastq \
        [SpeciesName].R1_unpaired.fastq \
        [SpeciesName].R2_paired.fastq \
        [SpeciesName].R2_unpaired.fastq \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

```
We will use two paired outputs as the input files for nexpolish step. In this script, [TruSeq3-PE.fa](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/TruSeq3-PE.fa) is provided here as the adaptor sequences.

Then, create a folder for polishing. Within this folder, we need two configuration files for polishing. The first one is `run.cfg`:
```bash
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 128
genome = assembly.fasta #genome file
genome_size = auto
workdir = ./polish1
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwatask = best
```
The second one is `sgs.fofn`, which contains two lines of where the trimmed paired illumina reads located:
```bash
/lustre/scratch/nhu/salix_assembly/nanopore_polish/illumina_reads/Salix-exigua-SE920_S5_L002_R1_001_paired.fastq
/lustre/scratch/nhu/salix_assembly/nanopore_polish/illumina_reads/Salix-exigua-SE920_S5_L002_R2_001_paired.fastq
```

After creating these two files, check all the required file are in the correct location (raw assembly, illumina reads). Then, run submission scripts:
```bash
#!/bin/bash
#SBATCH -J nextpolish
#SBATCH -o log/%x.o%j
#SBATCH -e log/%x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128

source ~/.bashrc

date

../NextPolish/nextPolish run.cfg

date

```

The output will be in `./polish1` folder named `genome.nextpolish.fasta`. We will rerun the polishing step for 3 more rounds. Each time, remember to use the polished genome from previous round as the input genome file for polishing. Also, remember to change the output directory name. Normally, after 4 rounds of polishing, the genome size will not change a lot.
