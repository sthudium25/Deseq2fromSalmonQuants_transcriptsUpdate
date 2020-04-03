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
samplenames_subset <- samplenames[c(2:7)] #from the list of files, I select the ones I want for the LPS alone and LPS-RMM comparison (from the list of file names I selected above), and put them in the correct order
samplenames_subset
#Write out the treatment conditions: (the levels tells the analysis which conditions to set as a control to compare to)
#Treatment <- factor(c("Control", "Control", "Control", "Dot1Li_4h","Dot1Li_4h","Dot1Li_4h"), levels=c("Control","Dot1Li_4h"))
#Alternatively can 'repeat' a name to avoid having to retype:
Treatment <- factor(c(rep("KO",3),rep("WT",3)), levels=c("KO","WT"))
Treatment
colData <- data.frame(samplenames_subset, Treatment)
colData #this is a table that has the exact file names of the Salmon output, and the other columns are descriptors of the experimental design that will be important for the DESeq2 analysis later on

#Now we can build a vector which points to our quantification files using this column of coldata. We use names to name this vector with the run IDs as well.
files <- file.path(dir,colData$samplenames_subset,"quant.sf")
names(files) <- colData$samplenames_subset
head(files,6)

#use Tximport to read in files correctly. dim gives dimension readout. Should be the number of lines and the number of samples
txi<- tximport(files, type="salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreAfterBar = TRUE)
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
cts <- as_tibble(txi$counts, rownames = NA)
cts <- rownames_to_column(cts)

# Rename the column names to something more manageable and make sure to store the key 
# somewhere in case you want to know which data set a specific number came from

colnames(cts) <- c("GeneID", "KO.1", "KO.2", "KO.3", "WT.1", "WT.2", "WT.3")
head(cts)

# Gather the columns into a tidy format for easier manipulation. 
cts_tidy <- gather(cts, key = "sample", value = "count", -GeneID)
head(cts_tidy)

# Add a factor column that denotes which group, WT or KO, the measured count is from.
cts_tidy <- cts_tidy %>%
   mutate(condition = as.factor(grepl("KO.*", sample))) 

# There's probably a better way around this but in order to label the conditions with WT and KO, 
# I've loaded the plyr package whcih has this handy revalue function. However, we then need to immediately
# remove plyr from the library because it interferes with dplyr (which I found out the hard way)
library(plyr)
cts_tidy$condition <- revalue(cts_tidy$condition, c("FALSE" = "WT", "TRUE" = "KO"))
head(cts_tidy)
detach(package:plyr)

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
  pivot_wider(names_from = condition, values_from = mean)

# I'm not sure this is the best way to proceed, but I am now splitting the df
# into two separate dfs, one for WT and one for KO. I'm also selecting only the 
# necessary rows so that when the tables are ultimately joined, it won't be a mess

cts_split <- split(x = cts_tidy.2, f = cts_tidy.2$condition)

# We need to cut off the decimals as we do in the Deseq pipeline. I just used the same procedure
# There's definitely a way to do this in one line but that's a later issue

cts_split$WT$GeneID <- make.unique(substr(cts_split$WT$GeneID, 1, 18))
cts_split$KO$GeneID <- make.unique(substr(cts_split$KO$GeneID, 1, 18))

#Now, we will build a DESeqDataSet from the matrices in tx

dds <- DESeqDataSetFromTximport(txi, colData, ~ Treatment)
str(dds)
#My favorite of these transformation is the vst, mostly because it is very fast, and provides transformed (nearly log-scale) data which is robust to many problems associated with log-transformed data (for more details, see the DESeq2 workflow or vignette ).
#blind=FALSE refers to the fact that we will use the design in estimating the global scale of biological variability, but not directly in the transformation:
vst <- vst(dds, blind=FALSE)
plotPCA(vst, "Treatment")

#We will chop off the version number of the gene IDs, so that we can better look up their annotation information later.
#However, we have a few genes which would have duplicated gene IDs after chopping off the version number, so in order to proceed we have to also use make.unique to indicate that some genes are duplicated. (It might be worth looking into why we have multiple versions of genes with the same base ID coming from our annotation.)
head(dds)
table(duplicated(substr(rownames(dds),1,18)))
rownames(dds) <- make.unique(substr(rownames(dds),1,18)) 
head(dds)
?map2


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
geneIDs <- substr(rownames(res_dds_over1), 1, 18)
# running mapIDs: collect gene symbols for the ensembl names in your geneID list
library(org.Mm.eg.db)
gene_symbols <- mapIds(org.Mm.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#add gene symbols as a new column to your res file
res_dds_over1$GeneSymbol <- gene_symbols

par("mar")
par(mar=c(1,1,1,1))
#make volcano plot
with(res_dds_over1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))
# Add colored points: red if padj<0.05. (Other options are for orange of log2FC>1, green if both)
with(subset(res_dds_over1, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
#with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
res_dds_over1 <- as.data.frame(res_dds_over1)
res_dds_over1 <- res_dds_over1[ ,c(7, 1:6)]
head(res_dds_over1)
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

# I'm now going to try to join the mean count data to this above df.

res_cts_tidy <- res_over1_tidy %>% 
  left_join(cts_split$WT, by = c("gene" = "GeneID")) %>% 
  left_join(cts_split$KO, by = c("gene" = "GeneID"), suffix = c("_WT", "_KO")) %>% 
  select(-baseMean, log2FC = estimate, -sample_WT, -sample_KO, -condition_KO, -condition_WT)

head(res_cts_tidy)

write.csv((res_cts_tidy),
          file="/Volumes/LaCie/SequencingData/H2BE/old_brains/DESeq/Sam_analysis/1.H2beOldBrains_WTvKO_Deseq_transcript.2_UseingTxOut.csv")



