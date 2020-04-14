#Start here if first time running code on this computer
source("http://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("GenomicFeatures")
BiocManager::install("org.Mm.eg.db")
install.packages("readr")
library(GenomicFeatures)
download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.annotation.gtf.gz", "gencode.vM19.annotation.gtf.gz")
BiocManager::install("tximport")
library("tximport")
txdb_ms <- makeTxDbFromGFF("gencode.vM19.annotation.gtf.gz")
saveDb(txdb_ms, file="gencode.vM19.sqlite")
# next time you can just load with loadDb command below (no need to makeTxDb...)

#Start here if code has been run on this computer before (and packages are installed, databases made, etc.)
library("DESeq2")
library("edgeR")
library("tximport")
library("readr")
library(GenomicFeatures)

txdb_ms <- loadDb("/Users/Sammy/Desktop/scripts/gencode.vM19.sqlite") 
columns(txdb_ms)
k <- keys(txdb_ms, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb_ms, k, "GENEID", "TXNAME") #for every gene, tell me the transcripts that are associated with it
head(tx2gene)
dim(tx2gene)

dir <- file.path("/Volumes/LaCie/SequencingData/H2BE/old_brains/quants")
list.files(dir)
samplenames <- list.files(dir) #this is the directory with my Salmon outputs and I wirte their names into a vector 
samplenames
#samplenames_reorder <- samplenames[c(6,7,5,1,2,3,4)] #from the list of files, I select the ones I want for the untreated and LPs comparison (from the list of file names I selected above), and put them in the correct order
#samplenames_reorder

# for differential expression between control and one condition alone 
samplenames_subset <- samplenames[c(5:7, 2:4)] #from the list of files, I select the ones I want for the LPS alone and LPS-RMM comparison (from the list of file names I selected above), and put them in the correct order
samplenames_subset
#Write out the treatment conditions: (the levels tells the analysis which conditions to set as a control to compare to)
#Treatment <- factor(c("Control", "Control", "Control", "Dot1Li_4h","Dot1Li_4h","Dot1Li_4h"), levels=c("Control","Dot1Li_4h"))
#Alternatively can 'repeat' a name to avoid having to retype:
Treatment <- factor(c(rep("WT",3),rep("KO",3)), levels=c("WT","KO"))
Treatment
colData <- data.frame(samplenames_subset, Treatment)
colData #this is a table that has the exact file names of the Salmon output, and the other columns are descriptors of the experimental design that will be important for the DESeq2 analysis later on

#Now we can build a vector which points to our quantification files using this column of coldata. We use names to name this vector with the run IDs as well.
files <- file.path(dir,colData$samplenames_subset,"quant.sf")
names(files) <- colData$samplenames_subset
head(files,6)

#use Tximport to read in files correctly. dim gives dimension readout. Should be the number of lines and the number of samples
txi<- tximport(files, type="salmon", 
               tx2gene = tx2gene, 
               countsFromAbundance = "lengthScaledTPM", 
               ignoreAfterBar = TRUE)
dim(txi$abundance)
dim(txi$counts)
dim(txi$length)
head(txi$counts)
# At this point, we have a list object called txi that contains abundance, counts, and length information from 
# the aggregated Salmon quant.sf files. As the code is currectly written, the counts are TPMs that have been
# scaled to both the length of the gene (avg length of the detected transcripts for a given gene) and the size 
# of the library. Thus, we now ought to be able to compare counts of a gene across the 6 libraries that are 
# present in these dataset. My understanding is that we would like to look at counts BY CONDITION so I 
# attempting to separate the 6 columns (libraries) by mouse background (WT or KO). Then, we can get a group 
# average count for each gene. 

## NOTE: if running just the Deseq pipeline, skip to making the DeseqDataSet on line 76 


# covert the counts matrix from tximport into a data frame for easier usage.
library(tidyverse)
cts <- as_tibble(txi$counts, rownames = NA)
cts <- rownames_to_column(cts)

# Rename the column names to something more manageable and make sure to store the key 
# somewhere in case you want to know which data set a specific number came from

colnames(cts) <- c("GeneID", "WT.1", "WT.2", "WT.3", "KO.1", "KO.2", "KO.3")
head(cts)

# Gather the columns into a tidy format for easier manipulation. 
cts_tidy <- gather(cts, key = "sample", value = "count", -GeneID)
head(cts_tidy)

# Add a factor column that denotes which group, WT or KO, the measured count is from.
cts_tidy <- cts_tidy %>%
   mutate(condition = as.factor(grepl("KO.*", sample))) 

levels(cts_tidy$condition) <- c("WT", "KO")

# Now, we can perform calculations on different groupings of data In this case, we want the mean counts 
# per gene per condition, so we group on those variables and add a column, mean. Make sure to ungroup at the
# end or later use of functions may act weird. 

cts_tidy.2 <- cts_tidy %>%
 group_by(condition, GeneID) %>%
  mutate(mean = mean(count)) %>%
  ungroup()

# I think the next thing to do is to move the mean count data into a wider
# format so that the df has two columns WT and KO, underneath which you see
# the count (averaged across the three corresponding samples) per gene.

cts_tidy.2 <- cts_tidy.2 %>%
  pivot_wider(names_from = condition, values_from = mean) %>%
  mutate(condition = as.factor(grepl("KO.*", sample)))

levels(cts_tidy.2$condition) <- c("WT", "KO")
  
head(cts_tidy.2)
tail(cts_tidy.2)
cts_tidy.2 <- select(cts_tidy.2, -sample, -count)
# I'm not sure this is the best way to proceed, but I am now splitting the df
# into two separate dfs, one for WT and one for KO. I'm also selecting only the 
# necessary rows so that when the tables are ultimately joined, it won't be a mess

cts_split <- split(x = cts_tidy.2, f = cts_tidy.2$condition)

# We need to cut off the decimals as we do in the Deseq pipeline. I just used the same procedure
# There's definitely a way to do this in one line but that's a later issue

cts_split$WT$GeneID <- make.unique(substr(cts_split$WT$GeneID, 1, 18))
cts_split$KO$GeneID <- make.unique(substr(cts_split$KO$GeneID, 1, 18))
head(cts_split)
#Now, we will build a DESeqDataSet from the matrices in tx

dds <- DESeqDataSetFromTximport(txi, colData, ~ Treatment)
str(dds)
#My favorite of these transformation is the vst, mostly because it is very fast, and provides transformed (nearly log-scale) data which is robust to many problems associated with log-transformed data (for more details, see the DESeq2 workflow or vignette ).
#blind=FALSE refers to the fact that we will use the design in estimating the global scale of biological variability, but not directly in the transformation:
vst <- vst(dds, blind=FALSE)
plotPCA(vst, intgroup = "Treatment")

#We will chop off the version number of the gene IDs, so that we can better look up their annotation information later.
#However, we have a few genes which would have duplicated gene IDs after chopping off the version number, so in order to proceed we have to also use make.unique to indicate that some genes are duplicated. (It might be worth looking into why we have multiple versions of genes with the same base ID coming from our annotation.)
head(dds)
table(duplicated(substr(rownames(dds),1,18)))
rownames(dds) <- make.unique(substr(rownames(dds),1,18)) 
head(dds)

#Now we can run our differential expression pipeline. First, it is sometimes convenient to remove genes where all the samples have very small counts. It's less of an issue for the statistical methods, and mostly just wasted computation, as it is not possible for these genes to exhibit statistical significance for differential expression. 
# Here we count how many genes (out of those with at least a single count) have 3 samples with a count of 10 or more:
dds <- dds[rowSums(counts(dds)) > 0,]
keep_dds <- rowSums(counts(dds) >= 1) >= 3
table(keep_dds)
dds_over1 <- dds[keep_dds,] #filter them out

dds_over1 <- DESeq(dds_over1)
resultsNames(dds_over1)
ResName <- resultsNames(dds_over1)
ResName
ResName_input <- ResName[2]
ResName_input
res_dds_over1 <- results(dds_over1, name = ResName_input)
head(res_dds_over1)
res_dds_over1.sort <- res_dds_over1[order(res_dds_over1$pvalue),]

#plotMA(res_dds_over1, ylim=c(-5,5))

summary(res_dds_over1)

#add gene names (symbols) to deseq results file
#modify your deseq results (res) table to take off numbers after decimal point to allow for matching to this database, needed to install org.Mm.eg.db previously for this to work
#geneIDs <- substr(rownames(res_dds_over1), 1, 18)
# running mapIDs: collect gene symbols for the ensembl names in your geneID list
#library(org.Mm.eg.db)
#gene_symbols <- mapIds(org.Mm.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#add gene symbols as a new column to your res file
#res_dds_over1$GeneSymbol <- gene_symbols

par("mar")
par(mar=c(1,1,1,1))
#make volcano plot
with(res_dds_over1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))
# Add colored points: red if padj<0.05. (Other options are for orange of log2FC>1, green if both)
with(subset(res_dds_over1, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
#with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
res_dds_over1 <- as.data.frame(res_dds_over1)
#res_dds_over1 <- res_dds_over1[ ,c(7, 1:6)]
#head(res_dds_over1)
#row.names(res_dds_over1)
#a <- row.names(res_dds_over1)
#namesTrim <- gsub("\\|.*", "", a)
#row.names(res_dds_over1) <- namesTrim
#head(res_dds_over1)

## Now, rather than writing out the res_dds_over1 file as usual, we will pull in the grouped
## mean count data from out cts_split list. 

library(biobroom)
res_over1_tidy <- tidy.DESeqResults(res_dds_over1)
head(res_over1_tidy)

# I'm now going to try to join the mean count data to this above df. I am also removing/adding so
# that we're left with the columns that are most useful.

res_cts_tidy <- res_over1_tidy %>% 
  left_join(cts_split$WT, by = c("gene" = "GeneID")) %>% 
  left_join(cts_split$KO, by = c("gene" = "GeneID"), suffix = c("_WT", "_KO")) %>% 
  select(-baseMean, log2FC = estimate, -condition_KO, -condition_WT, -KO_WT, -WT_KO,
         WT_cts = WT_WT, KO_cts = KO_KO) %>%
  pivot_longer(cols = WT_cts:KO_cts, 
               names_to = "condition", 
               names_ptypes = list(
                 condition = factor()),
               values_to = "count")

# I had some trouble getting the gene names into the table earlier for some reason, so I'm putting
# them in here. It's likely that they can be added to res_dds_over1 as they are in the original
# pipeline
geneIDs <- substr(res_cts_tidy$gene, 1, 18)
library(org.Mm.eg.db)
gene_symbols <- mapIds(org.Mm.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#add gene symbols as a new column to your res file
res_cts_tidy$GeneSymbol <- gene_symbols
res_cts_tidy <- res_cts_tidy[ ,c(1,9,2:8)]

## Here I am adding a label of quartiles based on where the log2FC falls in the distribution
res_cts_tidy <- res_cts_tidy %>% mutate(quartile = ntile(log2FC, 4))
head(res_cts_tidy)

# Now, the dataset is ready for EDA. Since it's in a tidy format, ggplot2 will work well for many
# different types of visualizations for this data. 

write.csv((res_cts_tidy),
          file="/Volumes/LaCie/SequencingData/H2BE/old_brains/DESeq/Sam_analysis/1.H2beKOvsWT_deseqWithGeneExpr.csv")

## Generate a histogram of WT expression values to see if there is an elbow

res_cts_tidy %>%
  filter(condition == "WT_cts") %>%
  ggplot(aes(count)) +
    geom_histogram(binwidth = 0.05) +
    scale_x_log10() +
    xlab("WT expression count distribution") +
    geom_vline(xintercept = 5, color = "red") 

## Get distributions of expression levels by condition
res_cts_tidy %>%
  ggplot(aes(condition, count)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") 
  

## Trying another method to also show t.test results as well
## We want a two-sample t-test

res_cts_tidy %>% identify_ou

comparison <- list(c("WT_cts", "KO_cts"))
p <- ggboxplot(res_cts_tidy, 
            x = "condition", y = "count",
            color = "condition", 
            palette = "npg") +
  stat_compare_means(method = "t.test", 
                     label.x = 1.3) +
  scale_y_log10()

p

## Get the distribution of KO Log2FC vs WT

res_cts_tidy %>%
  ggplot() +
  geom_density(aes(x = log2FC, bw = 2.2))

## Get distribution of log2FCs, binned by quartiles 

res_cts_tidy %>%
  filter(condition == "WT_cts",
         count > 0) %>%
  ggplot(aes(x = log2FC, y = count, group = quartile)) +
    geom_boxplot() +
  facet_wrap(. ~ as.factor(quartile)) +
  scale_y_log10()
  
res_cts_tidy %>%
  filter(condition == "WT_cts",
         count > 0) %>%
  ggplot(aes(x = log2FC, y = count, group = quartile)) +
  geom_boxplot(varwidth = TRUE) +
  scale_y_log10()
  
## THis one gives a nice visual for all datapoints
res_cts_tidy %>%
  group_by(quartile) %>%
  filter(condition == "WT_cts",
         count > 0) %>%
  ggplot(aes(x = as.factor(quartile), y = count)) +
  geom_boxplot() +
  scale_y_log10() +
  ggtitle("WT counts vs KO Log2FC") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(caption = "Binned as quartiles of KO Log2FC vs WT") 

# Showing just data points with p.adjusted < 0.05
res_cts_tidy %>%
  filter(condition == "WT_cts",
         count > 0,
         p.adjusted < 0.05) %>%
  ggplot(aes(x = as.factor(quartile), y = count)) +
  geom_boxplot() +
  scale_y_log10() +
  ggtitle("WT counts vs Significant KO Log2FC") +
  theme(plot.title = element_text(hjust = 0.5))

## This code makes a plot that has the WT count distribution for significant UP and DOWN DE genes
## as well as the distribution of ALL WT COUNTS as a comparison. 
## First we make a temporary data frame to hold all WT counts with a minimum expression over 5
temp <- res_cts_tidy %>%
  filter(condition == "WT_cts",
         count > 0) 

## Then, we can begin a new filter on the original data to remove H2be, get counts over 5 and
## significant adjusted p values
## I'm plotting that data first and then in the second geom_boxplot I'm adding the full WT count
## distribution on top

res_cts_tidy %>%
  filter(GeneSymbol != "Hist2h2be",
         condition == "WT_cts",
         count > 0,
         p.adjusted < 0.05) %>%
  ggplot(aes(x = as.factor(quartile), y = count)) +
  geom_boxplot() +
  geom_boxplot(data = temp, aes(x = condition, y = count)) +
  scale_y_log10() +
  ggtitle("WT counts vs Significant KO Log2FC")+
  labs(caption = "WT Counts > 5; Signif Down (1), Up (4) < 0.05 p.adjust") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())


# Here I'm going to try to do the same thing but with unnamed genes removed 
# (i.e. GeneSymbol = NA)
temp2 <- res_cts_tidy %>%
          filter(condition == "WT_cts",
          count > 0,
          !is.na(GeneSymbol)) 

res_cts_tidy %>%
  filter(GeneSymbol != "Hist2h2be",
         !is.na(GeneSymbol),
         condition == "WT_cts",
         count > 0,
         p.adjusted < 0.05) %>%
  ggplot(aes(x = as.factor(quartile), y = count)) +
  geom_boxplot() +
  geom_boxplot(data = temp2, aes(x = condition, y = count)) +
  scale_y_log10() +
  ggtitle("WT counts vs Significant KO Log2FC_noNAs")+
  labs(caption = "WT Counts > 5; Signif Down (1), Up (4) < 0.05 p.adjust") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())




