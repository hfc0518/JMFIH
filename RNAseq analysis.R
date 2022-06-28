#### set working directory ####
setwd("/Insert path to Jingyi FIH folder/JingyiFIH/")
#### load libraries ####
library(ggplot2)
library(edgeR)
library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
#### load functions ####
save_pheatmap <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svg(filename)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#### read in table of raw counts and annotate with sample names to create annotated counts ####
counts <- read.table("Raw_32 samples_2.txt", header=T, stringsAsFactors=F, sep="\t") #raw counts table from Jingyi 20190312
ref_genes <- read.table("mouse reference genome_biomaRt.txt", sep=',', header = T) #load df containing gene symbols
counts <- merge(counts, ref_genes, by.x = "Geneid",
                by.y = "Ensembl.Gene.ID", all.x = TRUE, all.y = FALSE) #match geneids and corresponding gene symbols
counts[,1] <- counts[,39] #replace original geneids with gene symbols
counts[,1] <- toupper(counts[,1]) #capitalize gene names
head(counts)[1] #check the first few gene names
counts <- counts[!duplicated(counts$Geneid),] #remove rows with duplicated gene names
ind_NT <- which(colnames(counts) %like% "NT") #NT abbv for control group 
ind_TM <- which(colnames(counts) %like% "LPS") #TM abbv for treatment group 
counts_ann <- counts[,c(1,ind_NT,ind_TM)] #our annotated counts table contains gene symbols, NT and TM groups
print(colnames(counts_ann)) #check again the colnames
write.table(counts_ann, "annotatedcounts_NT_LPS_james_20190513.txt", sep = "\t", row.names = F) #save file

##### check between-sample distribution of the annotated counts #####
x <- read.table("annotatedcounts_NT_LPS_james_20190513.txt", header=T, sep= "\t", row.names = 1)
colnames(x)
head(rownames(x))
logx <- log2(x + 1)
logx <- as.data.frame(t(logx))
logx$samples = rownames(logx)
head(logx$samples)
df <- melt(logx, variable_name = "samples")
df$samples <- gsub("Hif1an.Hif1an","KO",df$samples)
df$samples <- gsub("NT", "CTL", df$samples)
df <- data.frame(df,Condition=substr(df$samples,11,12))
df <- subset(df, df$Condition=="WT"|df$Condition=="KO")
ggplot(df,aes(x = value, colour = samples, fill = samples))+
  ylim(c(0, 0.25)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))

#### prepare data for Normalization by trimmed mean of M values (TMM) ####
x <- read.table("annotatedcounts_NT_LPS_james_20190513.txt", header=T, sep= "\t", row.names = 1)
#replace Hif1an.Hif1an with KO annotation
colnames(x) <- gsub("Hif1an.Hif1an","KO",colnames(x))
#subset selected samples from the same condition, i.e., LPS, using %like% from data.table library
ind_LPS <- which(colnames(x) %like% "LPS") #get the columns with LPS-treated samples
x <- x[,ind_LPS] # subset the relevant columns
colnames(x) # double check rownames
head(rownames(x))
group <- factor(c(rep("1",4), rep("2",4))) #two groups - FIH-WT (WT) and FIH-KO (KO)
y <- DGEList(counts = x,group=group) #create a DGEList object for subsequent edgeR analysis
keep <- rowSums(cpm(y)>1) >= 4 #set the threshold for filtering out lowly expressed genes
table(keep) #check how many genes are filtered out
y <- y[keep, , keep.lib.sizes=FALSE] #filter out lowly expressed genes using the foregoing threshold
y <- calcNormFactors(y, method = "TMM") #TMM normalization

#### quality-check of normalized reads - using multi-dimension plot to visualize sample clustering ####
plotMDS(y, col=as.numeric(y$samples$group))
# density plot to visualize between sample distribution of counts
x <- cpm.DGEList(y, normalized.lib.sizes=TRUE, log=FALSE)
logx <- log2(x+1)
logx <- as.data.frame(t(logx))
logx$samples <- rownames(logx)
df <- melt(logx, variable_name = "samples")
df <- data.frame(df,Condition=substr(df$samples,11,12))
df <- subset(df, df$Condition=="WT"|df$Condition=="KO")
ggplot(df,aes(x=value,colour=samples,fill=samples))+
  ylim(c(0, 0.25)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))

#### estimate dispersion of the negative binomial distribution to account for variability between biological replicates ####
design <- model.matrix(~group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design)

#### use exactTest in the edgeR package for DEG analysis ####
et <- exactTest(y)
results_exactTest <- topTags(et, n = nrow(y), sort.by = "none")
et_table = results_exactTest$table
et_table$genes <- rownames(et_table)
toplot <- et_table
toplot$genes <- rownames(toplot)
# set thresholds for positive/negative hits - FDR < 0.05 and absolute logFC > 1
toplot$Threshold <- if_else(toplot$FDR<0.05 & toplot$logFC > 1 ,"pos_hit",
                           if_else(toplot$FDR<0.05 & toplot$logFC < -1 ,"neg_hit","non_hit"))
#### create heatmap figure ####
tmm <- as.data.frame(cpm(y)) #get the cpm
hits <- toplot$genes[which(toplot$Threshold %in% c('pos_hit','neg_hit'))] # get the genenames that we identified as hits
ind <- which(rownames(tmm) %in% hits) 
tmm <- tmm[ind,] # subset the rows that contain the identified hits
# draw heatmap using pheatmap function from 'pheatmap' package and save as svg
hm <- pheatmap((tmm), color = colorspace::diverge_hsv(100),cellwidth =15, cellheight= 15, fontsize_col = 12,fontsize_row = 12,cluster_cols = T, cluster_rows = T, 
         scale = 'row')
hm_name <- paste0("LPS_treated_RNAseq_heatmap_Figure6A",".svg")
# call the save_pheatmap function to save the heatmap as svg
save_pheatmap(hm, hm_name)

#### repeat the steps for the untreated condition ####
#### prepare data for Normalization by trimmed mean of M values (TMM) ####
x <- read.table("annotatedcounts_NT_LPS_james_20190513.txt", header=T, sep= "\t", row.names = 1)
#replace Hif1an.Hif1an with KO annotation
colnames(x) <- gsub("Hif1an.Hif1an","KO",colnames(x))
#subset selected samples from the same condition, i.e., LPS, using %like% from data.table library
ind_NT <- which(colnames(x) %like% "NT") #get the columns with non-treated samples
x <- x[,ind_NT] # subset the relevant columns
colnames(x) # double check rownames
head(rownames(x))
group <- factor(c(rep("1",4), rep("2",4))) #two groups - FIH-WT (WT) and FIH-KO (KO)
y <- DGEList(counts = x,group=group) #create a DGEList object for subsequent edgeR analysis
keep <- rowSums(cpm(y)>1) >= 4 #set the threshold for filtering out lowly expressed genes
table(keep) #check how many genes are filtered out
y <- y[keep, , keep.lib.sizes=FALSE] #filter out lowly expressed genes using the foregoing threshold
y <- calcNormFactors(y, method = "TMM") #TMM normalization
#### quality-check of normalized reads - using multi-dimension plot to visualize sample clustering ####
plotMDS(y, col=as.numeric(y$samples$group))
# density plot to visualize between sample distribution of counts
x <- cpm.DGEList(y, normalized.lib.sizes=TRUE, log=FALSE)
logx <- log2(x+1)
logx <- as.data.frame(t(logx))
logx$samples <- rownames(logx)
df <- melt(logx, variable_name = "samples")
df <- data.frame(df,Condition=substr(df$samples,11,12))
df <- subset(df, df$Condition=="WT"|df$Condition=="KO")
ggplot(df,aes(x=value,colour=samples,fill=samples))+
  ylim(c(0, 0.25)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))

#### estimate dispersion of the negative binomial distribution to account for variability between biological replicates ####
design <- model.matrix(~group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design)

#### use exactTest in the edgeR package for DEG analysis ####
et <- exactTest(y)
results_exactTest <- topTags(et, n = nrow(y), sort.by = "none")
et_table = results_exactTest$table
et_table$genes <- rownames(et_table)
toplot <- et_table
toplot$genes <- rownames(toplot)
# set thresholds for positive/negative hits - FDR < 0.05 and absolute logFC > 1
toplot$Threshold <- if_else(toplot$FDR<0.05 & toplot$logFC > 1 ,"pos_hit",
                            if_else(toplot$FDR<0.05 & toplot$logFC < -1 ,"neg_hit","non_hit"))
#### create heatmap figure ####
tmm <- as.data.frame(cpm(y)) #get the cpm
hits <- toplot$genes[which(toplot$Threshold %in% c('pos_hit','neg_hit'))] # get the genenames that we identified as hits
ind <- which(rownames(tmm) %in% hits) 
tmm <- tmm[ind,] # subset the rows that contain the identified hits
# draw heatmap using pheatmap function from 'pheatmap' package
hm <- pheatmap((tmm), color = colorspace::diverge_hsv(100),cellwidth =15, cellheight= 15, fontsize_col = 12,fontsize_row = 12,cluster_cols = T, cluster_rows = T, 
         scale = 'row')
# save the figure as svg 
hm_name <- paste0("NT_treated_RNAseq_heatmap_SFigure6A",".svg")
# call the save_pheatmap function to save the heatmap as svg
save_pheatmap(hm, hm_name)
