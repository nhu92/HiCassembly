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



