# Using package 'allHiC' in assembling nanopore assembly with HiC information
Nan Hu / Oct. 18, 2020

---

## Preview
This pipeline is designed for using HiC sequencing data with nanopore assemblies attempting fetching relatively chromosomal level assemblies. Most of ideas refer to [ALLHiC Github](https://github.com/tangerzhang/ALLHiC/wiki).

The entire workflow for nanopore assembly, NGS data polishment, and HiC chromatin signal incorpration are shown below:
![Main Workflow](https://github.com/gudusanjiao/HiCassembly/blob/main/miscellaneous/Workflow.png "Workflow")
Steps below are a pipeline started from raw nanopore aseembly to reach chromosomal level assembly (Yellow box above). It works fine with *Salix nigra* now. There have been issues with *Salix exigua* which we still working on it. 

## Step 1: Index raw nanopore assembly
This step take raw nanopore aseembly fasta file as input to create index file for both alignment and BAM file operations.
> When we call 'raw assembly' here, we refer to assembly using nanopore reads assembled by chosen aseemblers, and then polished by NextPolish using NGS read. Noted here, the 'raw assembly' is not 'raw', it is after 4 rounds of polishing.
If your genome assembly is already indexed previously due to other analyzing process, you may skip this step.
```bash
#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SN_allHiC_index
#$ -q omni
#$ -l h_vmem=5.3G
#$ -o log/01/$JOB_NAME.o$JOB_ID
#$ -e log/01/$JOB_NAME.e$JOB_ID
#$ -pe sm 1
#$ -P quanah

bwa index -a bwtsw SN_flye1_nextpolish4.fasta  
samtools faidx SN_flye1_nextpolish4.fasta
```

## Step 2: Align HiC reads to raw nanopore aseembly
(Occupied needs updating)

## Step 3: Filter alignment
(Occupied needs updating)

## Step 4: Merge BAM files
(Occupied needs updating)

## Step 5: Prune alletic links and weak signals (Optional. Only for polyploids)
(Occupied needs updating)

## Step 6: Partition: Assign contigs into a pre-defined number of groups
(Occupied needs updating)

## Step 7: Rescue: Assign unplaced contigs into partitioned clusters (Optional. Only for ployploids)
(Occupied needs updating)

## Step 8: Extract grouping inforation
(Occupied needs updating)

## Step 9: Optimize directions of contigs
(Occupied needs updating)

## Step 10: Build FASTA file
(Occupied needs updating)

## Step 11: Plot chromatin contacting map
(Occupied needs updating)
