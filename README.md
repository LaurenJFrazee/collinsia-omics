# collinsia-omics
code behind Frazee et al. 2021 paper accepted in journal 'G3: Genes, Genomes, Genetics'

'New genomic resources and comparative analyses reveal differences in floral gene expression in selfing and outcrossing Collinsia sister species'


This repository contains the following 7 files:

bash-script-command-line

STARmappingPass1script.sh

STARmappingPass2script.sh

featureCounts.R

CollinsiaRNAreadcountsonlyERGO.txt

Analyses & select visualizations.R

bp_cora.tab
 

This is how to use them to recreate our results:
1. Start with 'bash-script-command-line' which will use 'CORA.gff' (NOT AVAILABLE HERE ON GITHUB, too large?; contact me for access!), 'STARmappingPass1script.sh', and 'STARmappingPass1script.sh'. You will also need the genome and raw RNAseq reads, which are both stored as described in the manuscript.
2. Next, run 'featureCounts.R' which will use the the output (.bam) files from step 1 and 'CORA.gff' (again NOT AVAILABLE HERE ON GITHUB) as input. This step will produce 'CollinsiaRNAreadcountsonlyERGO.txt' as output.
3. Lastly, run 'Analyses & select visualizations.R' which use the output file ('CollinsiaRNAreadcountsonlyERGO.txt', provided) from step 2, 'bp_cora.tab', and
'CORA.gff' (again NOT AVAILABLE HERE ON GITHUB) as input.




Article summary

We present a high-quality genome assembly for Collinsia rattanii, a self-fertilizing herb. We used flower bud transcriptomes from C. rattanii and its predominantly outcrossing sister species, C. linearis, to explore the genomic basis of mating system and phenotypic evolution in Collinsia, a self-compatible genus. Transcriptional regulation of enzymes involved in pollen formation may influence floral traits that distinguish selfing and outcrossing Collinsia species through pleiotropic functions. These patterns provide clues about parallel evolution in selfing plants.


Abstract

The evolutionary transition from outcross-fertilization to self-fertilization is one of the most common in angiosperms and is often associated with a parallel shift in floral morphological and developmental traits, such as reduced flower size and pollen to ovule ratios, known as the ‘selfing syndrome’. How these convergent phenotypes arise, the extent to which their characteristics are shaped by selection, and the nature of their underlying genetic basis are unsettled questions in evolutionary biology. The genus Collinsia includes seven independent transitions from outcrossing or mixed-mating to high selfing rates accompanied by selfing syndrome traits. Accordingly, Collinsia represents an ideal system for investigating this parallelism, but requires genomic resource development. We present a high quality de novo genome assembly for the highly selfing species C. rattanii. To begin addressing the basis of selfing syndrome developmental shifts, we evaluate and contrast patterns of gene expression from floral transcriptomes across three stages of bud development for C. rattanii and its outcrossing sister species C. linearis. Relative to C. linearis, total gene expression is less variable among individuals and bud stages in C. rattanii. In addition, there is a common pattern among differentially expressed genes: lower expression levels that are more constant across bud development in C. rattanii relative to C. linearis. Transcriptional regulation of enzymes involved in pollen formation specifically in early bud development may influence floral traits that distinguish selfing and outcrossing Collinsia species through pleiotropic functions. Future work will include additional Collinsia outcrossing-selfing species pairs to identify genomic signatures of parallel evolution.


Keywords

Collinsia, RNA-seq, selfing syndrome, pollen, floral development, differential gene expression, DESeq2, dichogamy, evolutionary genomics, Hi-C scaffolding, parallel evolution


Authors

Lauren J. Frazee, Joanna Rifkin, Dinusha C. Maheepala, Alannie-Grace Grant, Stephen Wright, Susan Kalisz, Amy Litt, & Rachel Spigler
