#!/usr/bin/env bash
# =============================================================================
# 16S Microbiome Diversity - QIIME2 pipeline (Moving Pictures, human body sites)
# Reads -> demux -> ASVs (DADA2) -> tree -> diversity -> PERMANOVA -> taxonomy.
# Companion R pipeline: microbiome_16s.R (dada2 + phyloseq, MiSeq SOP).
# =============================================================================
# SETUP (one time): install the QIIME2 amplicon distribution, then `conda activate qiime2`.
#   (If the env solve fails on deblur/sortmerna, strip those lines from the .yml:
#    sed -i -E '/deblur|sortmerna/d' <the-conda>.yml   -- deblur is an unused denoiser.)
# Run: conda activate qiime2 && bash microbiome_16s_qiime2.sh
# View any .qzv at https://view.qiime2.org
# =============================================================================
set -euo pipefail
mkdir -p data_q results_q
cd data_q

## ---- 1. Get data (multiplexed reads + barcodes + metadata) ----
[ -f sample-metadata.tsv ] || wget -O sample-metadata.tsv \
  "https://data.qiime2.org/2024.10/tutorials/moving-pictures/sample_metadata.tsv"
mkdir -p emp-single-end-sequences
[ -f emp-single-end-sequences/barcodes.fastq.gz ] || wget -O emp-single-end-sequences/barcodes.fastq.gz \
  "https://data.qiime2.org/2024.10/tutorials/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz"
[ -f emp-single-end-sequences/sequences.fastq.gz ] || wget -O emp-single-end-sequences/sequences.fastq.gz \
  "https://data.qiime2.org/2024.10/tutorials/moving-pictures/emp-single-end-sequences/sequences.fastq.gz"

## ---- 2. Import + demultiplex (barcodes -> per-sample reads) ----
qiime tools import --type EMPSingleEndSequences \
  --input-path emp-single-end-sequences --output-path emp-single-end-sequences.qza
qiime demux emp-single --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv --m-barcodes-column barcode-sequence \
  --o-per-sample-sequences demux.qza --o-error-correction-details demux-details.qza
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

## ---- 3. Denoise to ASVs (DADA2) ----
# --p-trunc-len chosen from demux.qzv quality profile (120 for this V4 single-end data).
qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 --p-trunc-len 120 \
  --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza

## ---- 4. Phylogenetic tree (for UniFrac / Faith's PD) ----
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

## ---- 5. Core diversity metrics (alpha + beta + PCoA), rarefied to even depth ----
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza \
  --p-sampling-depth 1103 --m-metadata-file sample-metadata.tsv --output-dir ../results_q/core-metrics

## ---- 6. Group significance (PERMANOVA) on body-site ----
qiime diversity beta-group-significance \
  --i-distance-matrix ../results_q/core-metrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv --m-metadata-column body-site \
  --o-visualization ../results_q/unweighted-unifrac-body-site.qzv --p-pairwise
qiime diversity alpha-group-significance \
  --i-alpha-diversity ../results_q/core-metrics/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization ../results_q/shannon-group-significance.qzv

## ---- 7. Taxonomy (reference-based vsearch; classifier-free, version-robust) ----
# NOTE: the pre-trained gg/silva sklearn classifiers were moved off data.qiime2.org for 2024.10,
# so we classify by SEARCH against the SILVA 515F/806R reference (no version-pinned model needed).
[ -f silva-seqs-515-806.qza ] || wget -O silva-seqs-515-806.qza \
  "https://data.qiime2.org/2024.10/common/silva-138-99-seqs-515-806.qza"
[ -f silva-tax-515-806.qza ] || wget -O silva-tax-515-806.qza \
  "https://data.qiime2.org/2024.10/common/silva-138-99-tax-515-806.qza"
qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seqs.qza \
  --i-reference-reads silva-seqs-515-806.qza \
  --i-reference-taxonomy silva-tax-515-806.qza \
  --p-threads 4 \
  --o-classification taxonomy.qza --o-search-results vsearch-hits.qza
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv --o-visualization ../results_q/taxa-bar-plots.qzv

## ---- 8. Export ASV table to a plain TSV ----
qiime tools export --input-path table.qza --output-path ../results_q/exported-table
biom convert -i ../results_q/exported-table/feature-table.biom \
  -o ../results_q/asv_table.tsv --to-tsv

## ---- Interpretation ----
# core-metrics gives alpha (Shannon, Faith PD) + beta (Bray-Curtis, weighted/unweighted UniFrac) + PCoA.
# Expectation: human body sites host distinct communities -> tight PCoA clusters by body-site, a
# significant PERMANOVA (unweighted UniFrac), gut the most diverse. The taxa barplot shows the
# biomarker taxa per site. Caveats: compositional data; rarefaction discards depth; 16S -> genus, not
# function; PERMANOVA can confound location vs dispersion.
