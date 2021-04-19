#! /bin/bash

for i in "lin1-1" "lin1-2" "lin1-3" "lin2-1" "lin2-2" "lin2-3" "lin3-1" "lin3-2" "lin3-3" "rat1-1" "rat1-2" "rat1-3" "rat2-1" "rat2-2" "rat2-3" "rat3-1" "rat3-2" "rat3-3"
do
STAR --genomeDir CrattaniiGenome \
--readFilesIn /$i\_R1.fastq /$i\_R2.fastq \
--outFileNamePrefix CollinsiaRNAalignment/Pass2/$i \
--outSAMtype BAM SortedByCoordinate \
--sjdbFileChrStartEnd CollinsiaRNAalignment/Pass1/lin1-1SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin1-2SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin1-3SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin2-1SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin2-2SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin2-3SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin3-1SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin3-2SJ.out.tab \
CollinsiaRNAalignment/Pass1/lin3-3SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat1-1SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat1-2SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat1-3SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat2-1SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat2-2SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat2-3SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat3-1SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat3-2SJ.out.tab \
CollinsiaRNAalignment/Pass1/rat3-3SJ.out.tab \
--runThreadN 16
done