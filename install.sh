#!/bin/bash
apt-get install -y samtools bcftools bwa bedtools r-base
Rscript -e 'install.packages(c("tidyverse", "ggplot2", "cowplot", "plotly", "htmlwidgets", "vcfR"), repos="https://cloud.r-project.org")'

# Get script dir, posix version
a="/$0"; a=${a%/*}; a=${a:-.}; a=${a#/}/; sdir=$(cd $a; pwd)
echo 'export stilPCR="${sdir}/process.sh"' >> ~/.bash_profile