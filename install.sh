#!/bin/bash
apt-get install -y samtools bcftools bwa bedtools r-base
Rscript -e 'install.packages(c("tidyverse", "ggplot2", "cowplot", "plotly", "htmlwidgets", "vcfR"), repos="https://cloud.r-project.org")'

# Get script dir, posix version
a="/$0"; a=${a%/*}; a=${a:-.}; a=${a#/}/; sdir=$(cd $a; pwd)
echo 'export stilPCR="${sdir}/process.sh"' >> ~/.bash_profile

# primer3
sudo apt-get install -y build-essential g++ cmake git-all
git clone https://github.com/primer3-org/primer3.git primer3
cd primer3/src
make
make install
cd

#seqfold
python3 -m pip install seqfold