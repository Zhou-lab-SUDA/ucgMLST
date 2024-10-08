## 1.2 Installation

**For Chinese users who are unable to access Google Drive, please contact us for the download link of pre-build database**

```shell
# step1. create a conda environment and install the dependencies
conda create -n ucgMLST
conda activate ucgMLST
conda install -c conda-forge -c bioconda python=3.8 pigz bindash diamond hmmer minimap2 iqtree ete3 click scikit-learn pandas pyarrow

# step2. download the ucgMLST package
INSTALL_PATH=/path/to/install/ucgMLST
mkdir $INSTALL_PATH && cd $INSTALL_PATH
git clone https://github.com/Zhou-lab-SUDA/ucgMLST.git
chmod +x $INSTALL_PATH/ucgMLST/modules/*.py

# step3. download our pre-build database (about 5.8 GB)
# For Chinese users, please contact us for the download link if you were unable to access Google Drive
# google drive shared link: https://drive.google.com/file/d/1z8QZ4YUGBYBMZH0L8OOhOtu7_zyn4CJN/view?usp=drive_link
wget -c -O ucgMLST_database.tar.gz "https://drive.usercontent.google.com/download?id=1z8QZ4YUGBYBMZH0L8OOhOtu7_zyn4CJN&export=download&authuser=0&confirm=t"
DB_LOCATION=/path/to/database/location
mkdir $DB_LOCATION && cd $DB_LOCATION
tar zxf ucgMLST_database.tar.gz

# step4 (optional): add ucgMLST to Environmental Variables 
echo "export PATH=\$PATH:$INSTALL_PATH/module" >> ~/.bashrc
```

## 1.3 Quick Start

```shell
# 0. activate the environment before use ucgMLST
conda activate ucgMLST

# 1. profile taxon from metagenomic data, use file path of genoQuery.py if you didn't do installation step4 
genoQuery.py -q read1.fq.gz -q read2.fq.gz -d /path/to/database -o result_dir

# 2. get taxon profiles table in tab seprate format from one or more result_dir of genoQuery.py
genoCompile.py -o profile_compare result_dir1 result_dir2
```
