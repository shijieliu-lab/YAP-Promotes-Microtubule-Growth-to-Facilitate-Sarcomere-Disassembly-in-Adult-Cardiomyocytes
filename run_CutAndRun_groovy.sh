#!/bin/bash
#BSUB -W 11:00
#BSUB -n 2
#BSUB -R "span[ptile=2]"
#BSUB -M 1000
#BSUB -e %J.err
#BSUB -o %J.out

#module load STAR/2.7.9
#module load bowtie2/2.4.2
#module load python3/3.7.1
#module load sratoolkit/3.0.0
#module load deeptools/2.0.0
#module load samtools/1.9.0

source ~/.bash_profile
bpipe run CutAndRun_workflow_uniq.groovy
