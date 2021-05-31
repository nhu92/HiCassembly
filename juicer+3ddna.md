# HiC Assembly using juicer + 3D-DNA
Nan Hu, May, 2021

---

## Software preparation
1. Juicer
2. 3D-DNA
3. last

### Install Juicer and configuration
Due to SLURM job submission system in HPCC of TTU, we need to use Juicer version for SLURM. After we clone the whole Juicer repository, we need to copy `<Juicer Dir>/SLURM/scripts` to our owrking directory.

However, our job submission system has a few different command settings that does not compatible to Juicer. Thus, there are some files need to be edited. I have files done for this step, so you can replace these files under `scripts` folder.

Then, we should build up our directory structure like below, assuming our working directory's name is `HiCassembly`:
```
-- HiCassembly
 |
 |--  scripts             copied from <Juicer Dir>/SLURM/scripts
 |--  restriction_sites
 |--  fastq
 |--  references
 |--  splits
 |--  chrom.sizes
```
For `references` folder, we need to copy our raw nanopore assembly fasta file into it.

For `fastq` folder, we need to copy our raw HiC sequencing reads into it.

For `restriction_sites` folder, we need to generate a file include all restrcition enzyme cutting sites on our raw assembled genome. There is a python script under `<Juicer Dir>/misc` named `generate_site_positions.py`. Here is the submission scripts:
```bash
#!/bin/bash
#SBATCH -J rez_size
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

python generate_site_positions.py MboI <Species Name> <raw nanopore assembly fasta file>

```
'MboI' is the restriction enzyme we used to digest our genome for HiC sequencing. Species name is our customed abbreviation of our target species. I used SAEX for Salix exigua and SANI for Salix nigra.

`chrom.sizes` is a file contains all the contig/scaffold name and a side column of the legnth for each segments. It can be generated from following commands:
```bash
grep ">" <raw nanopore assembly fasta file> | sed 's/>//g' > contig.namelist
python seqlength.py -f <raw nanopore assembly fasta file> -l contig.namelist > <output file name>

```
seqlength.py is available in miscellaneous folder.

For splits folder, we need to split fastq files into pieces for software to run. Here are some commands to generate the splits. This step takes very long time to finish.
```bash
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/<.fastq files>
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/<.fastq files>

```
Do this for all fastq files.

### Install 3D-DNA and configuration

