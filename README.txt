23 Feb 2021 -- by Lauren J. Frazee

This folder contains the following 8 files:
File S1. bash-script-command-line
File S2. CORA.gff
File S3. STARmappingPass1script.sh
File S4. STARmappingPass2script.sh
File S5. featureCounts.R
File S6. CollinsiaRNAreadcountsonlyERGO.txt
File S7. Analyses & select visualizations.R
File S8. bp_cora.tab
 
This is how to use them to recreate our results:
1. Start with 'bash-script-command-line' which will use 'CORA.gff',
'STARmappingPass1script.sh', and 'STARmappingPass1script.sh'. You will also need the
genome and raw RNAseq reads, which are both stored as described in the manuscript
(Frazee et al., submitted to G3).
2. Next, run 'featureCounts.R' which will use the the output (.bam) files from
step 1 and 'CORA.gff' (again) as input. This step will produce
'CollinsiaRNAreadcountsonlyERGO.txt' as output.
3. Lastly, run 'Analyses & select visualizations.R' which use the output file
('CollinsiaRNAreadcountsonlyERGO.txt', provided) from step 2, 'bp_cora.tab', and
'CORA.gff' as input. 