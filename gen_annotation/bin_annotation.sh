#!/bin/sh
#SBATCH --time=1-10:15 -n12 -p dque


python full_annotation.py all_bins_CpG_counts_multi_chr19.csv chr19.gene chr19 58617616
