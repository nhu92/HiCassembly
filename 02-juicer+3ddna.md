# HiC Assembly using juicer + 3D-DNA
Nan Hu, May, 2021

---

## Software preparation
1. [Juicer](https://github.com/aidenlab/juicer)
2. [3D-DNA](https://github.com/aidenlab/3d-dna)
3. last

## Install Juicer and configuration
Juicer does not require to install. In order to use Juicer, we only need to run `git clone` from author's github page. There are some [dependencies](https://github.com/aidenlab/juicer#readme) for this software.

Becasue SLURM job submission system is installed in HPCC of TTU, we need to use Juicer version for SLURM. After we clone the whole Juicer repository, we need to copy `[Juicer Dir]/SLURM/scripts` to our working directory. Then, we should install Juicertools by commands under scripts folder:
```
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
```
> (updated 07/11/2022) Check the newest version of juicer tools. Currently, I use v1.22.01.

However, our job submission system has a few different command settings that does not compatible to Juicer. Thus, there are some files need to be edited. I have files done for this step, so you can replace these files under `scripts` folder. ([juicer.sh](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/juicer.sh), [split_rmdups.awk](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/split_rmdups.awk))

Then, we should build up our directory structure like below, assuming our working directory's name is `HiCassembly`:
```
-- HiCassembly
 |
 |--  scripts             copied from [Juicer Dir]/SLURM/scripts
 |--  restriction_sites   generated from raw nanopore assembly by generate_site_positions.py
 |--  fastq               contains fastq file
 |--  references          contains raw nanopore assembly
 |--  splits              contains fastq file splits
 |--  chrom.sizes         a tablized file has contig names and sizes
```
For `references` folder, we need to copy our raw nanopore assembly fasta file into it.

For `fastq` folder, we need to copy our raw HiC sequencing reads into it. Juicer only treat with unzipped reads. Thus, if input fastq file is `.fastq.gz`, we need to unzip it by:
```bash
gunzip [.fastq.gz file]
```

For `restriction_sites` folder, we need to generate a file include all restrcition enzyme cutting sites on our raw assembled genome. There is a python script under `[Juicer Dir]/misc` named `generate_site_positions.py`. Here is the submission scripts:
```bash
#!/bin/bash
#SBATCH -J rez_size
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

python generate_site_positions.py MboI [Species Name] [raw nanopore assembly fasta file]

```
'MboI' is the restriction enzyme we used to digest our genome for HiC sequencing. Species name is our customed abbreviation of our target species. I used SAEX for Salix exigua and SANI for Salix nigra. After this step, copy the output file back to `restriction_sites` folder.

`chrom.sizes` is a file contains all the contig/scaffold name and a side column of the legnth for each segments. It can be generated from following commands:
```bash
python seqlength.py -f [raw nanopore assembly fasta file] > [output file name]

```
[seqlength.py](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/seqlength.py) is available in miscellaneous folder.

For `splits` folder, we need to split fastq files into pieces for software to run. Here are some commands to generate the splits. This step takes long time to finish. You may need to run a submission scripts or run under `interactive` command.
```bash
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/[.fastq files]
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/[.fastq files]

```
Do this for all fastq files.

## Generate HiC hitmaps by Juicer

After configurating the software, we should be able to run Juicer for our HiC initial assembly. Go to scripts folder then create a submission script like below:
```bash
#!/bin/bash
#SBATCH -J juicer
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 6
#SBATCH -n 128

source ~/.bashrc

bash juicer.sh -d .. -z ../references/[raw nanopore assembly fasta file] \
               -q nocona \ 
               -y ../restriction_sites/[output of restriction enzyme sites] \ 
               -t 128 \ 
               -p ../chrom.sizes \ 
               -l nocona

```
This step takes very very long time to finish. You can check the progress by check the error files in `debug` folder. We will get a file called `merged_nodups.txt` under `aligned` folder as our input file for 3D-DNA software.

## Install 3D-DNA and configuration
3D-DNA also does not required to install while it need some packages to work. Here is a [checklist](https://github.com/aidenlab/3d-dna#readme).

The structure of our working directory is below, assuming we have another working directory named `HiCscaffolding`:
```
-- HiCscaffolding
 |
 |-- 3d-dna            git clone from 3d-dna
 |-- raw_assembly      contains raw nanopore assembly 
 |-- hic               contains 'merged_nodups.txt' from Juicer
 |-- [results folder]  a custome folder to save results
```
Our `merged_nodups.txt` file may contain illegal form of data that causing skipping steps in 3d-dna scaffolding. We need to run following command to remove those lines.
```bash
awk 'NF==16' merged_nodups.txt > new_merged_nodups.txt
mv merged_nodups.txt old_merged_nodups.txt
mv new_merged_nodups.txt merged_nodups.txt
```

Our raw assembly need some indexing by running following scripts:
```bash
#!/bin/bash
#SBATCH -J indexing
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 16

bwa index -a bwtsw [raw assembly] 
samtools faidx [raw assembly] 
```

## Generating HiC scaffolding results (with polishing)
In order to run this software in parallel way, you may load `gnu parallel` module from HPCC. You can check it with `module spider parallel` and find the latest version of `gnu parallel` module to load.

After generating approporate structure of 3d-dna working directory, go to your custome result folder and create following file:
```bash
#!/bin/bash
#SBATCH -J 3ddnar10
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 128

bash ../3d-dna/run-asm-pipeline.sh -r 10 ../raw_assembly/[raw assembly fasta] ../hic/merged_nodups.txt

```
`-r` parameter is the round of polishment (corrections) you want to apply to assembly before scaffolding. Every additional polishment round requires 45 minutes to 1 hour to finish. There will be a `.fasta` file under result directory when it finishes.

## All to all alignment
All to all alignment is applied to examine the assembly quality by comparing assembly to well-done closed species genome assembly, as we assume there was no significant chromosomal rearrangement between closed related species.

Before we apply all to all alignment, we need to select assembled chromosomal-length scaffolds from HiC assemblies in order to avoid short sequences alignment. We can check the length of the assembly by [seqlength.py](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/seqlength.py) also. For example, we expected the number of chromosome in *Salix* would be 19. Thus, to cover major scaffolds and contigs, we can select top 30 longest sequences from `.FINAL.fasta` by:
```bash
grep ">" [.FINAL.fasta] | head -30 | sed 's/>//g' > top30.list
python pick_seq_list.py -f [.FINAL.fasta] -l top30.list -o [output FASTA]
```

[pick_seq_list.py](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/pick_seq_list.py) is in miscellaneous folder. The output will be used as the querry genome file for lastz software. 

We use lastz to run all-to-all alignment by `multiple` mode. Here is the command:
```bash
#!/bin/bash
#SBATCH -J lastz
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64

lastz SUBJECT_GENOME.fasta[multiple] QUERRY_GENOME.fasta--format=MAF --chain --gapped --transition --maxwordcount=4 --exact=100 --step=20 > RESULTS.maf

```
Normally, it will take 40 min to run the alignment for a genome with 300M size. When `.maf` results generated, we can input it into plotting scripts for plotting:
```bash
#!/bin/bash
#SBATCH -J lastz_plot
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 64

/home/nhu/software/last-1133/scripts/last-dotplot RESULT.maf ALL2ALL.png

```

