# HiC Assembly using juicer + 3D-DNA
Nan Hu, May, 2021

---

## Software preparation
1. Juicer
2. 3D-DNA
3. last

## Install Juicer and configuration
Juicer does not require to install. In order to use Juicer, we only need to run `git clone` from author's github page. There are some required environment for this software.

Becasue SLURM job submission system is installed in HPCC of TTU, we need to use Juicer version for SLURM. After we clone the whole Juicer repository, we need to copy `[Juicer Dir]/SLURM/scripts` to our working directory.

However, our job submission system has a few different command settings that does not compatible to Juicer. Thus, there are some files need to be edited. I have files done for this step, so you can replace these files under `scripts` folder.

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

For `fastq` folder, we need to copy our raw HiC sequencing reads into it.

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
grep ">" [raw nanopore assembly fasta file] | sed 's/>//g' > contig.namelist
python seqlength.py -f [raw nanopore assembly fasta file] -l contig.namelist > [output file name]

```
seqlength.py is available in miscellaneous folder.

For `splits` folder, we need to split fastq files into pieces for software to run. Here are some commands to generate the splits. This step takes very long time to finish.
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

bash juicer.sh -d .. -z ../references/[raw nanopore assembly fasta file] -q nocona -y ../restriction_sites/[output of restriction enzyme sites] -t 128 -p ../chrom.sizes -l nocona

```
This step takes very very long time to finish. You can check the progress by check the error files in `debug` folder. We will get a file called `merged_nodups.txt` under `aligned` folder as our input file for 3D-DNA software.

## Install 3D-DNA and configuration
3D-DNA also does not required to install and it need some packages to work. Here is a checklist.

The structure of our working directory is below:
```

```
