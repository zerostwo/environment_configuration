#!/bin/bash
# -----------------------------------------------------------------------------
# Basic configuration
WORKING_DIR=$1

# -----------------------------------------------------------------------------
# R
mamba create -p ${WORKING_DIR}/r/4.4.1 -c conda-forge r-base=4.4.1 r-seurat r-devtools r-tidyverse -y
# -----------------------------------------------------------------------------
# Download data
# mamba create -p ${WORKING_DIR}/sratoolkit/2.5.7 -c daler sratoolkit=2.5.7 -y
mamba create -p ${WORKING_DIR}/parallel-fastq-dump/0.6.7 -c bioconda parallel-fastq-dump=0.6.7 -y
mamba create -p ${WORKING_DIR}/entrez-direct/24.0 -c bioconda entrez-direct=24.0 -y
# -----------------------------------------------------------------------------
# Quality control
mamba create -p ${WORKING_DIR}/fastp/1.0.1 -c bioconda fastp=1.0.1 -y
mamba create -p ${WORKING_DIR}/falco/1.2.5 -c bioconda falco=1.2.5 -y
mamba create -p ${WORKING_DIR}/fastqc/0.12.1 -c bioconda fastqc=0.12.1 -y
mamba create -p ${WORKING_DIR}/fastq-pair/1.0 -c bioconda fastq-pair=1.0 -y
mamba create -p ${WORKING_DIR}/fastq-screen/0.16.0 -c bioconda fastq-screen=0.16.0 -y
mamba create -p ${WORKING_DIR}/trim-galore/0.6.10 -c bioconda trim-galore=0.6.10 -y
mamba create -p ${WORKING_DIR}/trimmomatic/0.40 -c bioconda trimmomatic=0.40 -y
# -----------------------------------------------------------------------------
# Mapping
mamba create -p ${WORKING_DIR}/hisat2/2.2.1 -c bioconda hisat2=2.2.1 -y
mamba create -p ${WORKING_DIR}/star/2.7.11b -c bioconda star=2.7.11b -y
mamba create -p ${WORKING_DIR}/bowtie2/2.5.4 -c bioconda bowtie2=2.5.4 -y
mamba create -p ${WORKING_DIR}/bwa/0.7.19 -c bioconda bwa=0.7.19 -y
mamba create -p ${WORKING_DIR}/bbmap/39.52 -c bioconda bbmap=39.52 -y
# -----------------------------------------------------------------------------
# Quantification
mamba create -p ${WORKING_DIR}/subread/2.1.1 -c bioconda subread=2.1.1 -y
mamba create -p ${WORKING_DIR}/htseq/2.0.9 -c bioconda htseq=2.0.9 -y
mamba create -p ${WORKING_DIR}/cufflinks/2.2.1 -c bioconda cufflinks=2.2.1 -y
mamba create -p ${WORKING_DIR}/stringtie/3.0.3 -c bioconda stringtie=3.0.3 -y
mamba create -p ${WORKING_DIR}/salmon/1.10.3 -c bioconda salmon=1.10.3 -y
mamba create -p ${WORKING_DIR}/kallisto/0.51.1 -c bioconda kallisto=0.51.1 -y
mamba create -p ${WORKING_DIR}/rsem/1.3.3 -c bioconda rsem=1.3.3 -y
# -----------------------------------------------------------------------------
# Alternative splicing
mamba create -p ${WORKING_DIR}/rmats/4.3.0 -c bioconda rmats=4.3.0 -y
mamba create -p ${WORKING_DIR}/rmats2sashimiplot/3.0.0 -c bioconda rmats2sashimiplot=3.0.0 -y
mamba create -p ${WORKING_DIR}/suppa/2.4 -c bioconda suppa=2.4 -y
# -----------------------------------------------------------------------------
# Others
mamba create -p ${WORKING_DIR}/multiqc/1.33 -c bioconda multiqc=1.33 -y
mamba create -p ${WORKING_DIR}/bedtools/2.31.1 -c bioconda bedtools=2.31.1 -y
mamba create -p ${WORKING_DIR}/samtools/1.23 -c bioconda samtools=1.23 -y
mamba create -p ${WORKING_DIR}/snakemake/9.14.5 -c bioconda snakemake=9.14.5 -y
mamba create -p ${WORKING_DIR}/sambamba/1.0.1 -c bioconda sambamba=1.0.1 -y
mamba create -p ${WORKING_DIR}/picard/3.4.0 -c bioconda picard=3.4.0 -y
mamba create -p ${WORKING_DIR}/umi_tools/1.1.6 -c bioconda umi_tools=1.1.6 -y
mamba create -p ${WORKING_DIR}/environment-modules/5.1.1 -c bioconda environment-modules=5.1.1 -y
# -----------------------------------------------------------------------------
# MeRIP-seq
mamba create -p ${WORKING_DIR}/homer/5.1 -c bioconda homer=5.1 -y
mamba create -p ${WORKING_DIR}/macs2/2.2.9.1 -c bioconda macs2=2.2.9.1 -y
mamba create -p ${WORKING_DIR}/meme/5.5.9 -c bioconda meme=5.5.9 -y
mamba create -p ${WORKING_DIR}/piranha/1.2.1 -c bioconda piranha=1.2.1 -y
# -----------------------------------------------------------------------------
# ChIP-seq
mamba create -p ${WORKING_DIR}/deeptools/3.5.6 -c bioconda deeptools=3.5.6 -y
# -----------------------------------------------------------------------------
# DNA methylation
# mamba create -p ${WORKING_DIR}/bismark/0.25.1 -c bioconda bismark=0.25.1 -y
# -----------------------------------------------------------------------------
# SNP
# mamba create -p ${WORKING_DIR}/gatk4/4.6.2.0 -c bioconda gatk4=4.6.2.0 -y
# -----------------------------------------------------------------------------
# Nanopore sequencing
# mamba create -p ${WORKING_DIR}/nanopolish/0.14.0 -c bioconda nanopolish=0.14.0 -y
# -----------------------------------------------------------------------------
# Genome
conda create -p ${WORKING_DIR}/blast/2.17.0 -c bioconda blast=2.17.0 -y
conda create -p ${WORKING_DIR}/rmblast/2.14.1 -c bioconda rmblast=2.14.1 -y
# -----------------------------------------------------------------------------
# Structure
mamba create -p ${WORKING_DIR}/autodock-vina/1.1.2 -c bioconda autodock-vina=1.1.2 -y
mamba create -p ${WORKING_DIR}/gromacs/2025.4 -c bioconda gromacs=2025.4 -y
mamba create -p ${WORKING_DIR}/openbabel/3.1.1 -c bioconda openbabel=3.1.1 -y
# -----------------------------------------------------------------------------
# Single-cell RNA-seq
mamba create -p ${WORKING_DIR}/velocyto.py/0.17.17 -c bioconda velocyto.py=0.17.17 -y
mamba create -p ${WORKING_DIR}/cellrank/1.5.1 -c bioconda cellrank=1.5.1 -y
mamba create -p ${WORKING_DIR}/scanpy/1.9.5 -c conda-forge scanpy=1.9.5 -y
mamba create -p ${WORKING_DIR}/scrublet/0.2.3 -c bioconda scrublet=0.2.3 -y
mamba create -p ${WORKING_DIR}/celltypist/1.6.1 -c bioconda celltypist=1.6.1 -y
# pySCENIC
PYSCENIC=${WORKING_DIR}/pyscenic/0.12.1
mamba create -p $PYSCENIC python=3.10 -y
${WORKING_DIR}/pyscenic/0.12.1/bin/pip install pyscenic==0.12.1
# CellPhoneDB
CELLPHONEDB=${WORKING_DIR}/cellphonedb/5.0.1
mamba create -p $CELLPHONEDB python=3.8 -y
${WORKING_DIR}/cellphonedb/5.0.1/bin/pip install cellphonedb==5.0.1
# CellBender
CELLBENDER=${WORKING_DIR}/cellbender/0.3.0
mamba create -p $CELLBENDER python=3.7 pytorch torchvision pytables -c conda-forge -y
cd $CELLBENDER
git clone https://github.com/broadinstitute/CellBender.git
${WORKING_DIR}/cellbender/0.3.0/bin/pip install -e CellBender
