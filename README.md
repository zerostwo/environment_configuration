# Environment configuration

📜Environment configuration script for bioinformatics analysis server.


```bash
mkdir logs
SOFTWARE_INSTALL_PATH=/your/software/install/path
nohup bash ./01_install_software.sh $SOFTWARE_INSTALL_PATH > logs/01_install_software.log &

DATABASE_DOWNLOAD_PATH=/your/database/download/path
nohup bash ./02_download_database.sh $DATABASE_DOWNLOAD_PATH > logs/02_download_database.log &

SOFTWARE_PATH=/opt
MODULEFILES_PATH=/opt/modules/5.0.1/modulefiles
bash ./03_modulefiles_generate.sh $SOFTWARE_PATH $MODULEFILES_PATH
```

```bash
# Human
INDEX_PATH=/DATA/public/index
GENOME_PATH=/DATA/public/genome/homo_sapiens/GRCh38.p13/GRCh38.p13.genome.fa
GTF_PATH=/DATA/public/genome/homo_sapiens/GRCh38.p13/gencode.v42.annotation.gtf
SPECIES=homo_sapiens/GRCh38.p13
N_CORES=20
BOWTIE2_PATH=/opt/bowtie2/2.5.1/bin

nohup bash ./04_genome_index.sh $INDEX_PATH $GENOME_PATH $GTF_PATH $SPECIES $N_CORES $BOWTIE2_PATH > logs/04_human_genome_index.log &

# Mouse
INDEX_PATH=/DATA/public/index
GENOME_PATH=/DATA/public/genome/mus_musculus/GRCm39/GRCm39.genome.fa
GTF_PATH=/DATA/public/genome/mus_musculus/GRCm39/gencode.vM31.annotation.gtf
SPECIES=mus_musculus/GRCm39
N_CORES=20
BOWTIE2_PATH=/opt/bowtie2/2.5.1/bin

nohup bash ./04_genome_index.sh $INDEX_PATH $GENOME_PATH $GTF_PATH $SPECIES $N_CORES $BOWTIE2_PATH > logs/04_mouse_genome_index.log &
```