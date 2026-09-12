# =============================================================================
# 16S Microbiome Diversity - R pipeline (dada2 + phyloseq + vegan)
# Reads -> ASVs -> taxonomy -> alpha/beta diversity -> PCoA -> PERMANOVA.
# Dataset: dada2 tutorial (MiSeq SOP, mouse gut), Early vs Late timepoints.
# Companion QIIME2 pipeline: microbiome_16s_qiime2.sh (Moving Pictures, body sites).
# =============================================================================
# Run in RStudio: Session -> Set Working Directory -> To Source File Location -> Source.
# First run installs Bioconductor packages + downloads reads and SILVA (~minutes).
# =============================================================================

## ---- Setup ----
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("dada2", "phyloseq")) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
library(dada2); library(phyloseq); library(vegan); library(ggplot2)
dir.create("data_R", showWarnings = FALSE); dir.create("results_R", showWarnings = FALSE)

## ---- 1. Get reads + SILVA reference ----
options(timeout = 1200)
if (!dir.exists("data_R/MiSeq_SOP")) {
  download.file("https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip", "data_R/sop.zip", mode = "wb")
  unzip("data_R/sop.zip", exdir = "data_R")
}
if (!file.exists("data_R/silva_nr99_v138.1_train_set.fa.gz"))
  download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
                "data_R/silva_nr99_v138.1_train_set.fa.gz", mode = "wb")
path <- "data_R/MiSeq_SOP"
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
cat("samples:", length(fnFs), "\n")

## ---- 2. Quality filter + trim ----
# Truncate low-quality read tails and cap expected errors so the denoiser sees clean reads.
filtFs <- file.path("data_R/filt", paste0(sample.names, "_F.fastq.gz"))
filtRs <- file.path("data_R/filt", paste0(sample.names, "_R.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160),
                     maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, multithread = FALSE)
print(head(out))

## ---- 3. Learn errors + denoise to ASVs ----
# dada2 learns the run's error model, then infers exact sequence variants (single-nucleotide resolution).
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

## ---- 4. Merge pairs -> ASV table -> remove chimeras ----
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
seqtab <- makeSequenceTable(merged)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
cat("ASVs:", ncol(seqtab.nochim), "| chimera-free reads kept:",
    round(sum(seqtab.nochim) / sum(seqtab), 3), "\n")

## ---- 5. Assign taxonomy (SILVA) ----
taxa <- assignTaxonomy(seqtab.nochim, "data_R/silva_nr99_v138.1_train_set.fa.gz", multithread = FALSE)

## ---- 6. Build phyloseq object (ASV table + taxonomy + metadata) ----
# Group = Early/Late (day <100 vs >=100), the dysbiosis timepoints of this mouse-gut study.
day <- as.integer(sub("^F3D", "", sample.names)); day[is.na(day)] <- 0
meta <- data.frame(row.names = sample.names, When = factor(ifelse(day < 100, "Early", "Late")))
rownames(seqtab.nochim) <- sapply(strsplit(basename(rownames(seqtab.nochim)), "_"), `[`, 1)  # match meta
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa), sample_data(meta))
ps <- prune_samples(sample_sums(ps) > 0, ps)
cat("phyloseq:", nsamples(ps), "samples x", ntaxa(ps), "ASVs\n")

## ---- 7. Alpha diversity (within-sample: Shannon) ----
alpha <- estimate_richness(ps, measures = c("Observed", "Shannon"))
alpha$When <- sample_data(ps)$When
write.csv(alpha, "results_R/alpha_diversity.csv")
ggsave("results_R/alpha_diversity.png",
       plot_richness(ps, x = "When", measures = c("Observed", "Shannon")) + geom_boxplot(),
       width = 7, height = 4, dpi = 150)

## ---- 8. Beta diversity -> PCoA -> PERMANOVA ----
# Relative abundance (16S is compositional); Bray-Curtis dissimilarity; PCoA to see, PERMANOVA to test.
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ord <- ordinate(ps.rel, method = "PCoA", distance = "bray")
ggsave("results_R/pcoa.png",
       plot_ordination(ps.rel, ord, color = "When") + geom_point(size = 3) +
         ggtitle("PCoA (Bray-Curtis) - mouse gut, Early vs Late"),
       width = 6, height = 5, dpi = 150)
bc <- phyloseq::distance(ps.rel, method = "bray")
perm <- adonis2(bc ~ When, data = data.frame(sample_data(ps.rel)))
print(perm)
capture.output(perm, file = "results_R/permanova.txt")

## ---- Interpretation ----
# Alpha (Shannon) per sample; beta (Bray-Curtis) between samples; PCoA visualises group separation;
# PERMANOVA tests whether Early vs Late community centroids differ (R^2 = variance explained, p from
# permutations). Expectation: the known post-weaning ("Late") shift structures the gut community, so
# Early and Late separate on the PCoA with a significant PERMANOVA. Caveats: compositional data; small
# n; 16S resolves to genus (not species/function); PERMANOVA can confound location vs dispersion.
