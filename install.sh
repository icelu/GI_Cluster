#!/bin/bash

cd GI_Cluster
# Revise INSTALLDIR accordingly
prefix=/mnt/projects/lub/workspace/GI_Cluster/program/
INSTALLDIR=$prefix

# Sample commands to download and organize the database files.
# Note that the names of folders and files must be in accord with those in GI_Feature.sh.
mkdir db
cd db

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
mkdir db/MOB
gunzip Pfam-A.hmm.gz Pfam-A.hmm.dat.gz
mv Pfam-A.hmm* db/MOB

wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz
mkdir -p db/RNA/CMs
mv Rfam.tar.gz db/RNA
mv db/RNA/*.cm db/RNA/CMs
cat *.cm > db/RNA/CMs/Rfam.cm

wget http://phast.wishartlab.com/phage_finder/DB/prophage_virus.db
mkdir db/PHAST
mv prophage_virus.db db/PHAST

wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
mkdir db/VFDB 
gunzip VFDB_setB_pro.fas.gz
mv VFDB_setB_pro.fas db/VFDB 

mkdir db/CARD 
mv card-data.tar.bz2 db/CARD 
tar xvfj card-data.tar.bz2

wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz   
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv
mkdir db/COG
gunzip prot2003-2014.fa.gz
mv prot2003-2014.fa cog2003-2014.csv db/COG



#===================================================
## Install required software and packages

mkdir program 
cd program

git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal
make install INSTALLDIR=$prefix

tar xvfz CodonWSourceCode_1_4_4.tar.gz
cd codonW/
./codonWinstall all

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
tar xvfz ncbi-blast-2.7.1+-x64-linux.tar.gz
# Add program/ncbi-blast-2.7.1+/bin to environment variable $PATH

wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz
tar xvfz hmmer-3.2.1.tar.gz
cd hmmer-3.2.1
./configure --prefix $prefix
make
make check               
make install             

wget http://circos.ca/distribution/circos-0.69-6.tgz
tar xvfz circos-0.69-6.tgz
# Add program/circos-0.69-6/bin to environment variable $PATH

wget http://eddylab.org/infernal/infernal-1.1.2.tar.gz
tar xvfz infernal-1.1.2.tar.gz
cd infernal-1.1.2
./configure --prefix $prefix
make
make check               
make install   

wget http://trna.ucsc.edu/software/trnascan-se-2.0.0.tar.gz
tar xvfz trnascan-se-2.0.0.tar.gz
cd tRNAscan-SE-2.0
./configure --prefix $prefix
make
make check               
make install

cd program/bin
wget http://wwwabi.snv.jussieu.fr/public/RepSeek/linux/repseek

tar xvf COGsoft.201204.tar
cd COGmakehash/
make
mv COGmakehash ../bin/
cd ../COGreadblast/
make
mv COGreadblast ../bin/
cd ../COGcognitor
make
mv COGcognitor ../bin/
cd ../COGtriangles
make
mv COGtriangles ../bin/
cd ../COGlse
make
mv COGlse ../bin/
cd ../

cd program/bin
wget https://github.com/icelu/GI_Prediction/blob/master/GI_SVM/scripts/jenks2.py
wget https://github.com/icelu/GI_Prediction/blob/master/GI_SVM/scripts/GI_SVM.py
 
wget ftp://ftp.sanger.ac.uk/pub/resources/software/alien_hunter/alien_hunter.tar.gz
# Install according to the documents 

# Remember to add program/bin to environment variable $PATH

# There are issues in installing the latest version of Ckmeans.1d.dp
wget https://cran.r-project.org/src/contrib/Archive/Ckmeans.1d.dp/Ckmeans.1d.dp_3.4.6-6.tar.gz
R CMD INSTALL Ckmeans.1d.dp_3.4.6-6.tar.gz
 
pip install scipy

