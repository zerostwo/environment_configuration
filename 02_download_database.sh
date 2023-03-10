#!/bin/bash
# -----------------------------------------------------------------------------
# Basic configuration
WORKING_DIR=$1

# -----------------------------------------------------------------------------
# 1. Genome data
# Data download from: https://www.gencodegenes.org/
# Human
GENOME_DIR=${WORKING_DIR}/genome/homo_sapiens/GRCh38.p13
mkdir -p $GENOME_DIR
cd $GENOME_DIR
axel -n 10 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.p13.genome.fa.gz
axel -n 10 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
axel -n 10 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gff3.gz
gzip -d ./*

# Mouse
GENOME_DIR=${WORKING_DIR}/genome/mus_musculus/GRCm39
mkdir -p $GENOME_DIR
cd $GENOME_DIR
axel -n 10 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.genome.fa.gz
axel -n 10 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gtf.gz
axel -n 10 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gff3.gz
gzip -d ./*

# -----------------------------------------------------------------------------
# 2. 10x Genomics references
# Data download from: https://www.10xgenomics.com/software
# Cell Ranger references - 2020-A (July 7, 2020)
REFDATA=${WORKING_DIR}/references
mkdir -p $REFDATA
cd $REFDATA
# Human
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# Mouse
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

# Cell Ranger ATAC GRCh38 Reference - 2020-A-2.0.0 (May 3, 2021)
# Human
curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
# Mouse
curl -O https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz

tar -xzvf ./*tar.gz

# -----------------------------------------------------------------------------
# 3. cisTarget resources
# Data download from: https://resources.aertslab.org/cistarget/
# Human
HUMAN_DIR=${WORKING_DIR}/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/
mkdir -p $HUMAN_DIR
cd $HUMAN_DIR
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather

# Mouse
MOUSE_DIR=${WORKING_DIR}/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/
mkdir -p $MOUSE_DIR
cd $MOUSE_DIR
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
# drosophila_melanogaster
FLY_DIR=${WORKING_DIR}/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc9nr/region_based/
mkdir -p $FLY_DIR
cd $FLY_DIR
wget https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc9nr/region_based/dm6-regions-11species.mc9nr.regions_vs_motifs.rankings.feather

# Motif2TF annotations
MOTIF2TF_DIR=${WORKING_DIR}/cistarget/motif2tf
mkdir -p $MOTIF2TF_DIR
cd $MOTIF2TF_DIR
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v8-nr.flybase-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.flybase-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.chicken-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

# TF lists
TF_LISTS=${WORKING_DIR}/cistarget/tf_lists
mkdir -p $TF_LISTS
cd $TF_LISTS
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_dmel.txt
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt


