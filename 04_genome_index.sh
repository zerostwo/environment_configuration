#!/bin/bash
# -----------------------------------------------------------------------------
# Basic configuration
WORKING_DIR=$1
GENOME_PATH=$2
GTF_PATH=$3
SPECIES=$4
N_CORES=$5
BOWTIE2_PATH=$6
# -----------------------------------------------------------------------------
# hisat2
module load hisat2/0.2.2
INDEX_PATH=${WORKING_DIR}/hisat2/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
hisat2-build -p $N_CORES $GENOME_PATH genome

# -----------------------------------------------------------------------------
# bowtie2
module load bowtie2/2.5.1
INDEX_PATH=${WORKING_DIR}/bowtie2/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
bowtie2-build --threads $N_CORES $GENOME_PATH genome

# -----------------------------------------------------------------------------
# kallisto
module load kallisto/0.48.0
INDEX_PATH=${WORKING_DIR}/kallisto/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
kallisto index $GENOME_PATH -i genome

# -----------------------------------------------------------------------------
# salmon
module load salmon/1.9.0
INDEX_PATH=${WORKING_DIR}/salmon/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
salmon index --threads $N_CORES -t $GENOME_PATH -i genome

# -----------------------------------------------------------------------------
# bismark
module load bismark/0.24.0
INDEX_PATH=${WORKING_DIR}/bismark/$SPECIES
mkdir -p $INDEX_PATH
cd $INDEX_PATH
ln -s $GENOME_PATH ./
bismark_genome_preparation \
	--parallel $N_CORES \
	--path_to_aligner $BOWTIE2_PATH \
	--verbose $INDEX_PATH

# -----------------------------------------------------------------------------
# star
module load star/2.7.10b
INDEX_PATH=${WORKING_DIR}/star/$SPECIES
STAR \
  --runThreadN $N_CORES \
  --runMode genomeGenerate \
  --genomeDir $INDEX_PATH \
  --genomeFastaFiles $GENOME_PATH \
  --sjdbGTFfile $GTF_PATH \
  --sjdbOverhang 149
