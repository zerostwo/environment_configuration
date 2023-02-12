# Environment configuration

📜Environment configuration script for bioinformatics analysis server.


```bash
mkdir logs
SOFTWARE_INSTALL_PATH=/your/software/install/path
nohup bash ./01_install_software.sh $SOFTWARE_INSTALL_PATH > logs/01_install_software.log &

DATABASE_DOWNLOAD_PATH=/your/database/download/path
nohup bash ./02_download_database.sh $DATABASE_DOWNLOAD_PATH > logs/02_download_database.log &
```