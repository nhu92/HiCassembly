# HiC Assembly using juicer + 3D-DNA
Nan Hu, May, 2021

---

## Software prepare
1. Juicer
2. 3D-DNA
3. last
Due to SLURM job submission system in HPCC of TTU, we need to use Juicer version for SLURM. After we clone the whole Juicer repository, we need to copy `<Juicer Dir>/SLURM/scripts` to our owrking directory.

However, our jib submission system has a few different command settings that does not compatible to Juicer. Thus, there are some files need to be edited. I have files done for this step, so you can replace these files under `scripts` folder.

