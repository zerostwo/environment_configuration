#!/bin/bash
# -----------------------------------------------------------------------------
# Basic configuration
WORKING_DIR=$1

# -----------------------------------------------------------------------------
# R
conda create -p ${WORKING_DIR}/r/4.2.2 -c conda-forge r-base=4.2.2 r-seurat r-devtools r-tidyverse r-hdf5r gcc libcxx -y
# -----------------------------------------------------------------------------
# Download data
conda create -p ${WORKING_DIR}/sratoolkit/2.5.7 -c bioconda sratoolkit=2.5.7 -y
conda create -p ${WORKING_DIR}/parallel-fastq-dump/0.6.7 -c bioconda parallel-fastq-dump=0.6.7 -y
# -----------------------------------------------------------------------------
# Quality control
conda create -p ${WORKING_DIR}/fastqc/0.11.9 -c bioconda fastqc=0.11.9 -y
conda create -p ${WORKING_DIR}/fastp/0.23.2 -c bioconda fastp=0.23.2 -y
conda create -p ${WORKING_DIR}/fastq-pair/1.0 -c bioconda fastq-pair=1.0 -y
conda create -p ${WORKING_DIR}/fastq-screen/0.15.2 -c bioconda fastq-screen=0.15.2 -y
conda create -p ${WORKING_DIR}/trim-galore/0.6.10 -c bioconda trim-galore=0.6.10 -y
conda create -p ${WORKING_DIR}/trimmomatic/0.39 -c bioconda trimmomatic=0.39 -y
# -----------------------------------------------------------------------------
# Mapping
conda create -p ${WORKING_DIR}/hisat2/0.2.2 -c bioconda hisat2=2.2.1 -y
conda create -p ${WORKING_DIR}/tophat/2.1.1 -c bioconda tophat=2.1.1 -y
conda create -p ${WORKING_DIR}/bowtie2/2.5.1 -c bioconda bowtie2=2.5.1 -y
conda create -p ${WORKING_DIR}/star/2.7.10b -c bioconda star=2.7.10b -y
conda create -p ${WORKING_DIR}/bwa/0.7.17 -c bioconda bwa=0.7.17 -y
conda create -p ${WORKING_DIR}/bbmap/39.01 -c bioconda bbmap=39.01 -y
# -----------------------------------------------------------------------------
# Quantification
conda create -p ${WORKING_DIR}/subread/2.0.3 -c bioconda subread=2.0.3 -y
conda create -p ${WORKING_DIR}/htseq/2.0.2 -c bioconda htseq=2.0.2 -y
conda create -p ${WORKING_DIR}/salmon/1.9.0 -c bioconda salmon=1.9.0 -y
conda create -p ${WORKING_DIR}/kallisto/0.48.0 -c bioconda kallisto=0.48.0 -y
conda create -p ${WORKING_DIR}/rsem/1.3.3 -c bioconda rsem=1.3.3 -y
conda create -p ${WORKING_DIR}/cufflinks/2.2.1 -c bioconda cufflinks=2.2.1 -y
conda create -p ${WORKING_DIR}/stringtie/2.2.1 -c bioconda stringtie=2.2.1 -y
# -----------------------------------------------------------------------------
# Alternative splicing
conda create -p ${WORKING_DIR}/rmats/4.1.2 -c bioconda rmats=4.1.2 -y
conda create -p ${WORKING_DIR}/rmats2sashimiplot/2.0.4 -c bioconda rmats2sashimiplot=2.0.4 -y
conda create -p ${WORKING_DIR}/suppa/2.3 -c bioconda suppa=2.3 -y
# -----------------------------------------------------------------------------
# Others
conda create -p ${WORKING_DIR}/multiqc/1.14 -c bioconda multiqc=1.14 -y
conda create -p ${WORKING_DIR}/bedtools/2.30.0 -c bioconda bedtools=2.30.0 -y
conda create -p ${WORKING_DIR}/snakemake/7.21.0 -c bioconda snakemake=7.21.0 -y
conda create -p ${WORKING_DIR}/samtools/1.16.1 -c bioconda samtools=1.16.1 -y
conda create -p ${WORKING_DIR}/sambamba/1.0 -c bioconda sambamba=1.0 -y
conda create -p ${WORKING_DIR}/picard/2.27.5 -c bioconda picard=2.27.5 -y
# conda create -p ${WORKING_DIR}/umi_tools/1.1.2 -c bioconda umi_tools=1.1.2 -y
conda create -p ${WORKING_DIR}/environment-modules/5.1.1 -c bioconda environment-modules=5.1.1 -y
# -----------------------------------------------------------------------------
# MeRIP-seq
conda create -p ${WORKING_DIR}/homer/4.11 -c bioconda homer=4.11 -y
conda create -p ${WORKING_DIR}/macs2/2.2.7.1 -c bioconda macs2=2.2.7.1 -y
conda create -p ${WORKING_DIR}/meme/5.4.1 -c bioconda meme=5.4.1 -y
conda create -p ${WORKING_DIR}/piranha/1.2.1 -c bioconda piranha=1.2.1 -y
# -----------------------------------------------------------------------------
# ChIP-seq
conda create -p ${WORKING_DIR}/deeptools/3.5.1 -c bioconda deeptools=3.5.1 -y
# -----------------------------------------------------------------------------
# DNA methylation
conda create -p ${WORKING_DIR}/bismark/0.24.0 -c bioconda bismark=0.24.0 -y
# -----------------------------------------------------------------------------
# SNP
conda create -p ${WORKING_DIR}/gatk4/4.3.0.0 -c bioconda gatk4=4.3.0.0 -y
# -----------------------------------------------------------------------------
# Nanopore sequencing
conda create -p ${WORKING_DIR}/nanopolish/0.14.0 -c bioconda nanopolish=0.14.0 -y
# -----------------------------------------------------------------------------
# Genome
conda create -p ${WORKING_DIR}/blast/2.13.0 -c bioconda blast=2.13.0 -y
conda create -p ${WORKING_DIR}/rmblast/2.11.0 -c bioconda rmblast=2.11.0 -y
# -----------------------------------------------------------------------------
# Structure
conda create -p ${WORKING_DIR}/autodock-vina/1.1.2 -c bioconda autodock-vina=1.1.2 -y
conda create -p ${WORKING_DIR}/gromacs/2021.3 -c bioconda gromacs=2021.3 -y
conda create -p ${WORKING_DIR}/openbabel/3.1.1 -c bioconda openbabel=3.1.1 -y
# -----------------------------------------------------------------------------
# Single-cell RNA-seq
conda create -p ${WORKING_DIR}/velocyto.py/0.17.17 -c bioconda velocyto.py=0.17.17 -y
conda create -p ${WORKING_DIR}/cellrank/1.5.1 -c bioconda cellrank=1.5.1 -y
conda create -p ${WORKING_DIR}/scanpy/1.7.2 -c bioconda scanpy=1.7.2 -y
conda create -p ${WORKING_DIR}/scrublet/0.2.3 -c bioconda scrublet=0.2.3 -y
conda create -p ${WORKING_DIR}/celltypist/1.3.0 -c bioconda celltypist=1.3.0 -y
# pySCENIC
PYSCENIC=${WORKING_DIR}/pyscenic/0.12.1
conda create -p $PYSCENIC python=3.10 -y
${WORKING_DIR}/pyscenic/0.12.1/bin/pip install pyscenic==0.12.1
# CellPhoneDB
CELLPHONEDB=${WORKING_DIR}/cellphonedb/3.1.0
conda create -p $CELLPHONEDB python=3.7 -y
${WORKING_DIR}/cellphonedb/3.1.0/bin/pip install cellphonedb==3.1.0
# CellBender
CELLBENDER=${WORKING_DIR}/cellbender/0.2.2
conda create -p $CELLBENDER python=3.7 pytorch torchvision pytables -c conda-forge -y
cd $CELLBENDER
git clone https://github.com/broadinstitute/CellBender.git
${WORKING_DIR}/cellbender/0.2.2/bin/pip install -e CellBender
