#!/bin/bash
# -----------------------------------------------------------------------------
# Basic configuration
WORKING_DIR=$1
GENOME_PATH=$2
GTF_PATH=$3
SPECIES=$4
N_CORES=$5

hisat2_build_path=/opt/hisat2/2.2.1/bin/hisat2-build
bowtie2_build_path=/opt/bowtie2/2.5.2/bin/bowtie2-build
kallisto_path=/opt/kallisto/0.50.0/bin/kallisto
salmon_path=/opt/salmon/1.10.2/bin/salmon
bismark_genome_preparation_path=/opt/bismark/0.24.2/bin/bismark_genome_preparation
star_path=/opt/star/2.7.11a/bin/STAR

# -----------------------------------------------------------------------------
# hisat2
INDEX_PATH=${WORKING_DIR}/hisat2/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
$hisat2_build_path -p $N_CORES $GENOME_PATH genome

# -----------------------------------------------------------------------------
# bowtie2
INDEX_PATH=${WORKING_DIR}/bowtie2/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
$bowtie2_build_path --threads $N_CORES $GENOME_PATH genome

# -----------------------------------------------------------------------------
# kallisto
INDEX_PATH=${WORKING_DIR}/kallisto/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
$kallisto_path index $GENOME_PATH -i genome

# -----------------------------------------------------------------------------
# salmon
INDEX_PATH=${WORKING_DIR}/salmon/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
$salmon_path index --threads $N_CORES -t $GENOME_PATH -i genome

# -----------------------------------------------------------------------------
# bismark
INDEX_PATH=${WORKING_DIR}/bismark/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
ln -s $GENOME_PATH ./
$bismark_genome_preparation_path \
	--parallel $N_CORES \
	--path_to_aligner $bowtie2_build_path \
	--verbose $INDEX_PATH

# -----------------------------------------------------------------------------
# star
INDEX_PATH=${WORKING_DIR}/star/$SPECIES
$star_path \
  --runThreadN $N_CORES \
  --runMode genomeGenerate \
  --genomeDir $INDEX_PATH \
  --genomeFastaFiles $GENOME_PATH \
  --sjdbGTFfile $GTF_PATH \
  --sjdbOverhang 149
