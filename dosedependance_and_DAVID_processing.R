needed <- subset(RESULT,select = c(Geneid, gene, gene_product, H37Rv_gene, logFC_DM1_vs_DM0, 
                                    PValue_DM1_vs_DM0, FDR_DM1_vs_DM0, logFC_DM2_vs_DM0,  PValue_DM2_vs_DM0, FDR_DM2_vs_DM0,
                                    logFC_DM3_vs_DM0,  PValue_DM3_vs_DM0, FDR_DM3_vs_DM0))
needed$diffexpressed <- 'NO'
needed$diffexpressed[needed$logFC_DM1_vs_DM0 > 1 & needed$PValue_DM1_vs_DM0 < 0.05|
                       needed$logFC_DM2_vs_DM0 > 1 & needed$PValue_DM2_vs_DM0 < 0.05|
                       needed$logFC_DM3_vs_DM0 > 1 & needed$PValue_DM3_vs_DM0 < 0.05] <- 'Up'

needed$diffexpressed[needed$logFC_DM1_vs_DM0 < -1 & needed$PValue_DM1_vs_DM0 < 0.05|
                       needed$logFC_DM2_vs_DM0 < -1 & needed$PValue_DM2_vs_DM0 < 0.05|
                       needed$logFC_DM3_vs_DM0 < -1 & needed$PValue_DM3_vs_DM0 < 0.05] <- 'Down'
overexpressed <- needed[needed$diffexpressed == "Up", 1:3]
underexpressed <- needed[needed$diffexpressed == "Down", 1:3]
overexpressed$Geneid

# dose-dependence
dose_overexpressed <- needed[needed$diffexpressed == "Up" & 
                               2*(needed$logFC_DM1_vs_DM0) < needed$logFC_DM2_vs_DM0 & 
                               2*(needed$logFC_DM2_vs_DM0) < 2*(needed$logFC_DM3_vs_DM0), 1:3]
dose_underexpressed <- needed[needed$diffexpressed == "Down" & 
                               2*(needed$logFC_DM1_vs_DM0) > needed$logFC_DM2_vs_DM0 & 
                               2*(needed$logFC_DM2_vs_DM0) > 2*(needed$logFC_DM3_vs_DM0), 1:3]


# Построим график, на котором будет видна зависимость уровня экспресии от концентрации, 
# то есть log2(FC) от номера 1, 2 или 3

# heatmap
for_heatmap <- subset(needed, select = c(logFC_DM1_vs_DM0, 
                                         logFC_DM2_vs_DM0,
                                         logFC_DM3_vs_DM0))
names(for_heatmap)[1] <- "DM1"
names(for_heatmap)[2] <- "DM2"
names(for_heatmap)[3] <- "DM3"
rownames(for_heatmap) <- needed$Geneid
for_heatmap <- for_heatmap[order(for_heatmap$DM1),]
for_heatmap <- na.omit(for_heatmap)
data <- as.matrix(for_heatmap)
heatmap(data)

## select genes for DAVID analisys
write.table(dose_overexpressed$Geneid,"filename.txt",sep="\t",row.names=FALSE)

## count genes but not pseudogenes
write.table(needed$Geneid,"all_locus_tags.txt",sep="\t",row.names=FALSE)

## data for bubbleplot


DM1_UP <- needed$Geneid[needed$logFC_DM1_vs_DM0 > 1 & needed$PValue_DM1_vs_DM0 < 0.05]
DM2_UP <- needed$Geneid[needed$logFC_DM2_vs_DM0 > 1 & needed$PValue_DM2_vs_DM0 < 0.05]
DM3_UP <- needed$Geneid[needed$logFC_DM3_vs_DM0 > 1 & needed$PValue_DM3_vs_DM0 < 0.05]

write.table(DM1_UP,"DM1_Up.txt",sep="\t",row.names=FALSE)
write.table(DM2_UP,"DM2_Up.txt",sep="\t",row.names=FALSE)
write.table(DM3_UP,"DM3_Up.txt",sep="\t",row.names=FALSE)

David_DM1 <- subset(Da_DM1_1,select = c(Term, Count, PValue, Fold.Enrichment, FDR))
David_DM1 <- David_DM1[David_DM1$FDR < 0.2, 1:4]
David_DM2 <- subset(Da_DM2,select = c(Term, Count, PValue, Fold.Enrichment, FDR))
David_DM2 <- David_DM2[David_DM2$FDR < 0.013, 1:4]
David_DM3 <- subset(Da_DM3,select = c(Term, Count, PValue, Fold.Enrichment, FDR))
David_DM3 <- David_DM3[David_DM3$FDR < 0.0242, 1:4]



typeof(David_DM1$FDR)

library(ggplot2)
library(dplyr)
library(gridExtra)

ggplot(David_DM1, aes(x=Fold.Enrichment, y=Term, size = Count, color = PValue)) +
  geom_point(alpha=0.7)

ggplot(David_DM2, aes(x=Fold.Enrichment, y=Term, size = Count, color = PValue)) +
  geom_point(alpha=0.7)

ggplot(David_DM3, aes(x=Fold.Enrichment, y=Term, size = Count, color = PValue)) +
  geom_point(alpha=0.7)
ggarrange(p1, p2, p3, ncol = 3, nrow = 1)


## merge diffexpression with gene names
names(All_locus_genes)[names(All_locus_genes) == "LOCUS_TAG"] <- "Geneid"
overexpressed_genes <- merge(overexpressed, All_locus_genes)
underexpressed_genes <- merge(underexpressed, All_locus_genes)
dose_over_genes <- merge(dose_overexpressed, All_locus_genes)
dose_under_genes <- merge(dose_underexpressed, All_locus_genes)

need_over <- subset(overexpressed_genes, select = c(Geneid, Name))
need_under <- subset(underexpressed_genes, select = c(Geneid, Name))

dose_need_over <- subset(dose_over_genes, select = c(Geneid, Name))
dose_need_under <- subset(dose_under_genes, select = c(Geneid, Name))  

write.table(dose_need_over,"dose_over.txt",sep="\t",row.names=TRUE)
write.table(dose_need_under,"dose_under.txt",sep="\t",row.names=TRUE)
  
  
  
  