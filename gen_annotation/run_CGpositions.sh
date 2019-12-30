#!/bin/sh
#SBATCH --time=16:00:00 -n24 -p bynode

python run_annotatoin.py chrsize

#for i in {1..22}
#do
#perl ../cgpositionFinder.pl chr$i.fa
#python ~/scripts/jack_bin/bin_CpG_counts.py 100 chr$i.CpG.positions.txt
#done

#for i in X Y
#do
#perl ../cgpositionFinder.pl chr$i.fa
#python ~/scripts/jack_bin/bin_CpG_counts.py 100 chr$i.CpG.positions.txt
#done

