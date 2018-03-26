#!/bin/bash
#SBATCH -J sstacks
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=5g
#SBATCH -t 14-0:0:0
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=15pl16@queensu.ca
#SBATCH -o output_sstacks%j.o

module load stacks/1.46

sstacks -c ./output/batch_1 -s ./output/1007 \
-s ./output/100 \
-s ./output/14 \
-s ./output/15 \
-s ./output/16 \
-s ./output/17 \
-s ./output/18 \
-s ./output/19 \
-s ./output/25 \
-s ./output/26 \
-s ./output/27 \
-s ./output/28 \
-s ./output/30 \
-s ./output/31 \
-s ./output/32 \
-s ./output/33 \
-s ./output/34 \
-s ./output/35 \
-s ./output/36 \
-s ./output/37 \
-s ./output/38 \
-s ./output/39 \
-s ./output/40 \
-s ./output/41 \
-s ./output/42 \
-s ./output/43 \
-s ./output/44 \
-s ./output/45 \
-s ./output/46 \
-s ./output/47 \
-s ./output/48 \
-s ./output/49 \
-s ./output/60 \
-s ./output/61 \
-s ./output/70 \
-s ./output/71 \
-s ./output/723 \
-s ./output/724 \
-s ./output/725 \
-s ./output/726 \
-s ./output/727 \
-s ./output/728 \
-s ./output/729 \
-s ./output/72 \
-s ./output/730 \
-s ./output/731 \
-s ./output/732 \
-s ./output/73 \
-s ./output/741 \
-s ./output/742 \
-s ./output/743 \
-s ./output/753 \
-s ./output/754 \
-s ./output/755 \
-s ./output/756 \
-s ./output/757 \
-s ./output/758 \
-s ./output/75 \
-s ./output/760 \
-s ./output/761 \
-s ./output/762 \
-s ./output/763 \
-s ./output/76 \
-s ./output/776 \
-s ./output/777 \
-s ./output/77 \
-s ./output/78 \
-s ./output/792 \
-s ./output/793 \
-s ./output/794 \
-s ./output/795 \
-s ./output/79 \
-s ./output/80 \
-s ./output/810 \
-s ./output/81 \
-s ./output/82 \
-s ./output/833 \
-s ./output/84 \
-s ./output/850 \
-s ./output/851 \
-s ./output/85 \
-s ./output/868 \
-s ./output/869 \
-s ./output/86 \
-s ./output/87 \
-s ./output/88 \
-s ./output/89 \
-s ./output/90 \
-s ./output/91 \
-s ./output/92 \
-s ./output/93 \
-s ./output/94 \
-s ./output/95 \
-s ./output/96 \
-s ./output/97 \
-s ./output/98 \
-s ./output/99 \
-s ./output/L_627 \
-o ./output -b 1
