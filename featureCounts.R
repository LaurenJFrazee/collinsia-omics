# Use this code after alignment from within R/RStudio through the HPC cluster

#install & load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")
library(Rsubread)

#Turn .bam files into summarized read count file
CollinsiaRNAreadcounts <- featureCounts(
  isGTFAnnotationFile = TRUE, 
  isPairedEnd = TRUE,   #yes these were paired end reads
  nthreads = 4,   #how much computing power to use
  GTF.featureType = "CDS",   #see manual and structure of GTF versus GFF3 files
  GTF.attrType = "Parent",   #see manual and structure of GTF versus GFF3 files
  annot.ext = "~/CrattaniiGenome/CORA.gff", 
  files = c("~/CollinsiaRNAalignment/Pass2/lin1-1Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin1-2Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin1-3Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin2-1Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin2-2Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin2-3Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin3-1Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin3-2Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/lin3-3Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat1-1Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat1-2Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat1-3Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat2-1Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat2-2Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat2-3Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat3-1Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat3-2Aligned.sortedByCoord.out.bam",
            "~/CollinsiaRNAalignment/Pass2/rat3-3Aligned.sortedByCoord.out.bam"))

#re-format featureCounts output & write to file
write.table(x=data.frame(CollinsiaRNAreadcounts$annotation[,c("GeneID")], CollinsiaRNAreadcounts$counts,
                         stringsAsFactors=FALSE),file="CollinsiaRNAreadcountsonlyERGO.txt",
            quote=FALSE,sep="\t",row.names=FALSE)
