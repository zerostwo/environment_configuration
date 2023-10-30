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
# Human GRCh38.p14
mkdir -p logs
INDEX_PATH=/datapool/index
GENOME_PATH=/datapool/genome/homo_sapiens/GRCh38.p14/GRCh38.p14.genome.fa
GTF_PATH=/datapool/genome/homo_sapiens/GRCh38.p14/gencode.v44.annotation.gtf
SPECIES=homo_sapiens/GRCh38.p14
N_CORES=10

nohup bash ./04_genome_index.sh $INDEX_PATH $GENOME_PATH $GTF_PATH $SPECIES $N_CORES > logs/04_human_GRCh38.p14_genome_index.log &

# Mouse GRCm38.p6 mm10
INDEX_PATH=/datapool/index
GENOME_PATH=/datapool/genome/mus_musculus/GRCm38.p6/GRCm38.p6.genome.fa
GTF_PATH=/datapool/genome/mus_musculus/GRCm38.p6/gencode.vM25.annotation.gtf
SPECIES=mus_musculus/GRCm38.p6
N_CORES=10

nohup bash ./04_genome_index.sh $INDEX_PATH $GENOME_PATH $GTF_PATH $SPECIES $N_CORES > logs/04_mouse_GRCm38.p6_genome_index.log &

# Mouse GRCm39
INDEX_PATH=/datapool/index
GENOME_PATH=/datapool/genome/mus_musculus/GRCm39/GRCm39.genome.fa
GTF_PATH=/datapool/genome/mus_musculus/GRCm39/gencode.vM33.annotation.gtf
SPECIES=mus_musculus/GRCm39
N_CORES=10

nohup bash ./04_genome_index.sh $INDEX_PATH $GENOME_PATH $GTF_PATH $SPECIES $N_CORES > logs/04_mouse_GRCm39_genome_index.log &
```
