# --- Stage 1: Install & Load Libraries ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

# --- Stage 2: Set Directory & Create Metadata ---
# Use the exact path to your folder
directory <- "C:/Users/Lola/Desktop/mapping/counts"

# Create a table linking filenames to conditions
sample_files <- c("BH_rep1.counts", "BH_rep2.counts", "BH_rep3.counts", 
                  "Serum_rep1.counts", "Serum_rep2.counts", "Serum_rep3.counts")

sample_names <- c("BH_1", "BH_2", "BH_3", "Serum_1", "Serum_2", "Serum_3")
sample_conditions <- c("BH", "BH", "BH", "Serum", "Serum", "Serum")

sample_table <- data.frame(sampleName = sample_names,
                           fileName = sample_files,
                           condition = sample_conditions)

# Set "BH" as the reference (the baseline we compare AGAINST)
sample_table$condition <- factor(sample_table$condition, levels = c("BH", "Serum"))

# --- Stage 3: Import into DESeq2 ---
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                  directory = directory,
                                  design = ~ condition)

# --- Stage 4: The Analysis ---
# This runs normalization, dispersion estimation, and the statistical tests
dds <- DESeq(dds)

# Get the results
res <- results(dds, contrast=c("condition", "Serum", "BH"))

# Filter for Significance (Padj < 0.05)
res_sig <- res[which(res$padj < 0.05), ]

# See how many genes were up/down regulated
summary(res)

# Install ggrepel if you want nice labels later, but for now:
res_df <- as.data.frame(res)

# Add a column for "Significant" to color the dots
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("black", "red")) +
  theme_minimal() +
  labs(title="Volcano Plot: Serum vs BH",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") # This marks the 0.05 threshold

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

### EGGNOG ###

# Load the annotation file (skip the header lines starting with #)
egg_anno <- read.delim("C:/Users/Lola/Desktop/eggnog/Galaxy3-[eggNOG Mapper on dataset 1_ annotations].tabular", 
                       skip = 4, header = FALSE)

# The columns are usually: V1 = query (ID), V10 = GOs
# Let's rename them for clarity
colnames(egg_anno)[c(1, 10)] <- c("locus_tag", "go_terms")
egg_anno <- egg_anno[, c("locus_tag", "go_terms")]

# Some genes have no GO terms (they have a '-'). Let's remove them.
egg_anno <- egg_anno[egg_anno$go_terms != "-", ]

# Assuming 'res_sig' is your table of significant genes from DESeq2
# and it has a column called 'locus_tag'
final_map <- merge(res_sig, egg_anno, by = "locus_tag")

# Because one gene can have many GO terms separated by commas, 
# we need to "split" them so each GO term has its own row
library(tidyr)
final_map_long <- separate_rows(final_map, go_terms, sep = ",")

colnames(as.data.frame(res_sig))
colnames(egg_anno)

# 1. Convert DESeq2 results to a data frame and move row names to 'locus_tag'
res_df <- as.data.frame(res_sig)
res_df$locus_tag <- rownames(res_df)

# 2. Now the merge will work because both tables have a "locus_tag" column
final_map <- merge(res_df, egg_anno, by = "locus_tag")

# 3. Check the first few rows to make sure they combined
head(final_map)

library(tidyr)

# Split the comma-separated GO terms into separate rows
final_long <- separate_rows(final_map, go_terms, sep = ",")

# Clean up any accidental spaces
final_long$go_terms <- trimws(final_long$go_terms)

# Look at the result - you should see the same gene appearing multiple times
head(final_long)

install("rrvgo")
library(rrvgo)

# 1. Create a named vector of scores (using -log10 of the adjusted p-value)
# This tells rrvgo which GO terms are the most significant
genes_to_plot <- final_long$go_terms
gene_scores <- setNames(-log10(final_long$padj), final_long$go_terms)

# 2. Calculate the Similarity Matrix
simMatrix <- calculateSimMatrix(genes_to_plot,
                                orgdb="org.Hs.eg.db", # Using human as a hierarchy reference
                                ont="BP", 
                                method="Rel")

# 3. Reduce the redundancy (grouping similar GO terms together)
reducedTerms <- reduceSimMatrix(simMatrix,
                                gene_scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

# 1. Replace p-values of 0 with the smallest possible number R can handle
# This prevents -log10(0) from becoming Infinity
final_long$padj[final_long$padj == 0] <- .Machine$double.xmin

# 2. Create the scores again
gene_scores <- setNames(-log10(final_long$padj), final_long$go_terms)

# 3. Remove any NA scores just in case
gene_scores <- gene_scores[!is.na(gene_scores)]

# 4. Re-run the reduction (using the cleaned scores)
reducedTerms <- reduceSimMatrix(simMatrix,
                                gene_scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

# 5. Check if you actually have terms to plot
# If this number is 0 or 1, a treemap won't work well
nrow(reducedTerms)

# 6. Try the plot again
treemapPlot(reducedTerms)

# 1. Re-calculate the reduction with a higher threshold (clustering more tightly)
reducedTerms_tight <- reduceSimMatrix(simMatrix,
                                      gene_scores,
                                      threshold=0.9, # Try 0.9 first
                                      orgdb="org.Hs.eg.db")

# 2. Check the new count - ideally we want this under 100 for a pretty treemap
nrow(reducedTerms_tight)

# 3. Try the plot again
treemapPlot(reducedTerms_tight)

# Sort by score and take the top 50 rows
reducedTerms_top <- reducedTerms_tight[order(reducedTerms_tight$score, decreasing = TRUE), ]
reducedTerms_top <- head(reducedTerms_top, 50)

# Plot the top 50
treemapPlot(reducedTerms_top)

scatterPlot(simMatrix, reducedTerms_tight)
