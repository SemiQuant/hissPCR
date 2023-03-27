#!/bin/bash
usage () { #echo -e to read \n
  echo "
Usage Options
  -t|--threads = number of threads (def = 1)
  -n|--sample_name = sample name (def = ${R1/R1*/})
  -1|--R1 = path to read 1 (required)
  -2|--R2 = path to read 2 (required)
  -f|--ref = path to reference fasta (required)
  -p|--primers = path to file containin primer sequences (optinal but recommended)
  -d|--scrip_dir = path to script dir (def = posix calculated)
  -o|--out_dir = path to output dir (def = cwd)
"
}

declare_globals () {
    # if same thing twice will take second one
    while [[ "$#" -gt 0 ]]
    do
        case $1 in
        -t|--threads)
        threads="$2"
        ;;
        -n|--sample_name)
        sample="$2"
        ;;
        -1|--R1)
        R1="$2"
        ;;
        -2|--R2)
        R2="$2"
        ;;
        -f|--ref)
        ref="$2"
        ;;
        -p|--primers)
        primers="$2"
        ;;
        -d|--script_dir)
        scrip_dir="$2"
        ;;
        -o|--out_dir)
        out_dir="$2"
        ;;
    esac
        shift
    done
}

declare_globals "$@"
sample=${sample:-${R1/R1*/}}
threads=${threads:-1}
primers="$primers" # requires samtools >=1.4
# Get script dir, posix version
a="/$0"; a=${a%/*}; a=${a:-.}; a=${a#/}/; sdir=$(cd $a; pwd)
script_path="${script_dir:-$sdir}"
out_dir="${out_dir:-$pwd}"
cd "$out_dir"
TRIM="${script_path}/Trimmomatic-0.39/trimmomatic-0.39.jar"

#check if programs installed
command -v bcftools >/dev/null 2>&1 || { echo >&2 "I require bcftools but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v R >/dev/null 2>&1 || { echo >&2 "I require R but it's not installed. Aborting."; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools but it's not installed. Aborting."; exit 1; }
if [ ! -f "$TRIM" ]; then echo "$TRIM not found!"; exit 1; fi



# bwa index
if [ -z "$ref" ]
then
  echo "I need a reference fasta."
  exit 1
fi

if [ -e "${ref}.sa" ] #check if indexed alread
then
  bwa index "$ref"
else
  echo "Found bwa index for $ref?"
fi



java -jar "$TRIM" PE \
  -threads $threads \
  "$R1" "$R2" \
  "${R1/.fastq.gz/.trimmed.fastq.gz}" "${R1/.fastq.gz/.trimmedUmpaired.fastq.gz}" \
  "${R2/.fastq.gz/.trimmed.fastq.gz}" "${R2/.fastq.gz/.trimmedUmpaired.fastq.gz}" \
  ILLUMINACLIP:"${script_path}/Trimmomatic-0.39/adapters/TruSeq3-PE.fa":2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36


R1="${R1/.fastq.gz/.trimmed.fastq.gz}"
R2="${R2/.fastq.gz/.trimmed.fastq.gz}"

bam="${sample}.bam"

bwa mem -t $threads "$ref" "$R1" "$R2" > "${sample}.sam"
samtools view -bS "${sample}.sam" | samtools sort - > "$bam"


if [ ! -z "$primers" ]
then
  samtools ampliconclip --tolerance 8 -u --rejects-file "${sample}_NOclipped.bam" -b "$primers" --hard-clip "$bam" | samtools fixmate -@ $threads -O bam - "${sample}_clipped.tmp.bam"
  samtools sort -@ $threads "${sample}_clipped.tmp.bam" -o "${sample}_clipped.bam"
  bam="${sample}_clipped.bam"
else
  echo "NOT CLIPPING PRIMERS!!"
fi

samtools index "$bam"
samtools flagstat "$bam" > "${sample}.flagstat.txt"
# samtools depth -a "${sample}.bam" -o "${sample}.depth.tsv"
bedtools genomecov -dz -ibam "$bam" > "${sample}.depth.tsv"

# samtools view -F 4 "$bam" | gawk '{ if ( and($2, 16) == 0 ) { strand="+" } else { strand="-" }; print $3 "\t" $4 "\t" $4 + length($10) "\t"  strand }' > "${sample}.counts.tsv"

# bedtools genomecov -bg -strand + -ibam "$bam" | awk 'BEGIN{FS=OFS="\t+"}{print $0 OFS }' > "${sample}.counts.tsv"
# samtools view -b -F 4 -f 16 "$bam" | bedtools genomecov -bg -ibam - | awk 'BEGIN{FS=OFS="\t-"}{print $0 OFS }' > "${sample}.counts.tsv"
# samtools view -b -F 4 -f 32 "$bam" | bedtools genomecov -bg -ibam - | awk 'BEGIN{FS=OFS="\t+"}{print $0 OFS }' >> "${sample}.counts.tsv"
samtools view -F 4 -f 16 "$bam" | gawk '{ print $3 "\t" $4 "\t" $4 + length($10) "\t-"}' > "${sample}.counts.tsv"
samtools view -F 4 -f 32 "$bam" | gawk '{ print $3 "\t" $4 "\t" $4 + length($10) "\t+"}' >> "${sample}.counts.tsv"

Rscript "${script_path}/stillPCR_plots.R" "${sample}.counts.tsv" "${sample}.depth.tsv"


# samtools ampliconstats -@ $threads "$primers" "${sample}_clipped.bam" -o "${sample}_astats.txt"
# plot-ampliconstats -size 1200,900 mydata "${sample}_astats.txt"

bcftools mpileup --fasta-ref "$ref" --max-depth 999999999 \
  --ignore-RG \
  --min-BQ 30 \
  --threads $threads \
  --excl-flags "UNMAP,SECONDARY,QCFAIL" \
  --max-idepth 9999999999 \
  -Ou \
  --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
  "$bam" |
  bcftools call -m --keep-masked-ref --threads $threads \
  --novel-rate 1e-100 \
  --pval-threshold 1.0 \
  --prior 0 \
  --keep-alts > "${sample}.vcf"

  
Rscript "${script_path}/variantAnalysis.R" "${sample}.vcf" "$ref"
