# hissPCR

Analysis of Illumina amplicon sequencing produced using hisPCR (link to protocols.io). The program will trim off the primer sequences (so they do not interfere with variant calling), produce numerous plots of the data, and calls variants. Outputs will be output in the current working directory.

## process.sh
### Process fastq and create plots/ call variants
>bash process.sh \
>  --R1 "test_data/read_R1_001.fastq.gz" \
>  --R2 "test_data/read_R2_001.fastq.gz" \
>  --ref "refs/BDQ_duplex.fasta" \
>  --primers "refs/primers.tsv"


| Flag              | Description                                                       |
|-------------------|-------------------------------------------------------------------|
| -t\|--threads     | number of threads (def = 1)                                       |
| -n\|--sample_name | sample name (def = ${R1/R1*/})                                    |
| -1\|--R1          | path to read 1 (required)                                         |
| -2\|--R2          | path to read 2 (required)                                         |
| -f\|--ref         | path to reference fasta (required)                                |
| -p\|--primers     | path to file containin primer sequences (optinal but recommended, require samtools ≥v1.4) |
| -d\|--script_dir  | path to script dir (def = posix calculated)                       |
| -o\|--out_dir     | path to output dir (def = cwd)                                    |


#### Primers bed file
This file is highly recommended as otherwise primers will be used to call variants (or the lack of variants)
The file must contain thee columns, in the first, the chromosome name matching the name in the reference fasta file, second, the start of where this primer binds to the reference sequence, and third, when the binding ends.
for example

>$cat primers.bed

```
Rv0678	0	20
Rv0678	421	441
Rv0678	459	479
Rv0678	967	987
```


## Required software
These can be installed on linux systems by running 
>./install.sh

- [samtools >=1.4](http://www.htslib.org/download/)
- [bcftools >=1.4](http://www.htslib.org/download/)
- [bwa >= 0.7.17](https://sourceforge.net/projects/bio-bwa/files/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)


R packages
- vcfR
- tidyverse
- ggplot2
- cowplot
- plotly
- htmlwidgets
