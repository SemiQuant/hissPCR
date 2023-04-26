#!/bin/bash

usage () {
  echo "
Usage Options
  -n|--name = name for amplicons (default = hisPCR_nameless)
  -c|--seq_cycles = Illuina sequecnign cycles (read length; default = 250)
  -h|--calc_hairpin = should the structure of the inner amplicon be caluclated, require seqfold and take significant time. Alternatively see NuPack and DNAmelt.
  -start|--start = start position for region of intrest on template (REQUIRED)
  -end|--end = end position for region of intrest on template (REQUIRED)
  -t|--template = template sequence (REQUIRED)
  -o|--out_dir = path to output directory
  -s|--Script_dir = path to script directory

# Required Software
  primer3
  seqfold (if calc_hairpin set)
"
}

declare_globals () {
    # if same thing twice will take second one
    while [[ "$#" -gt 0 ]]
    do
        case $1 in
        -d|--Script_dir)
        Script_dir="$2"
        ;;
        -c|--seq_cycles)
        seq_cycles="$2"
        ;;
        -n|--name)
        id="$2"
        ;;
        -h|--calc_hairpin)
        calc_hairpin="Y"
        ;;
        -t|--template)
        template="$2"
        ;;
        -s|--start)
        start="$2"
        ;;
        -e|--end)
        end="$2"
        ;;
        -o|--out_dir)
        out_dir="$2"
        ;;
    esac
        shift
    done
}

declare_globals "$@"

command -v primer3_core >/dev/null 2>&1 || { echo >&2 "I require primer3_core but it's not installed. Aborting."; exit 1; }
if [[ $calc_hairpin == "Y" ]]
then
    command -v seqfold >/dev/null 2>&1 || { echo >&2 "I require seqfold but it's not installed. Aborting."; exit 1; }
fi

if [ -z ${start+x} ]
then
    echo 'Please set the start position!'
    exit 1
fi

if [ -z ${end+x} ]
then 
    echo 'Please set the end position!'
    exit 1
fi

if [ -z ${template+x} ]
then 
    echo 'Please set the template sequence!'
    exit 1
fi

if [ ! -z ${out_dir+x} ]
then 
    cd "$out_dir"
fi


# set defults
Script_dir_tmp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Script_dir="${Script_dir:-$Script_dir_tmp}"
seq_cycles=${seq_cycles:-250}
id="${id:-hisPCR_nameless}"
calc_hairpin=${calc_hairpin:-N}
size_max=$(expr $seq_cycles \* 4 - 60)
size_min=$(expr $end - $start)


if [ $size_min -gt $size_max ]
then
    echo "Teamplate length too long for $id"
fi


get_primers () {
  echo "id=${2}
PRIMER_NUM_RETURN=1
SEQUENCE_ID=${2}
PRIMER_TASK=generic
SEQUENCE_TARGET=$3,$4
SEQUENCE_TEMPLATE=${1}
PRIMER_PRODUCT_SIZE_RANGE=${4}-${5}
PRIMER_MIN_SIZE=15
PRIMER_OPT_SIZE=20
PRIMER_MAX_SIZE=30
PRIMER_MIN_TM=58
PRIMER_OPT_TM=60
PRIMER_MAX_TM=65
PRIMER_GC_CLAMP=0
PRIMER_MAX_GC=80
PRIMER_OPT_GC_PERCENT=50
PRIMER_MIN_GC=30
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_MAX_HAIRPIN_TH=47.0
PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.0
PRIMER_MAX_END_STABILITY=100.0
PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=-1
PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=-1
PRIMER_EXPLAIN_FLAG=0
PRIMER_LIBERAL_BASE=0
PRIMER_FIRST_BASE_INDEX=0
PRIMER_MAX_TEMPLATE_MISPRIMING=-1.00
PRIMER_MAX_TEMPLATE_MISPRIMING_TH=-1.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=-1.00
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=-1.00
PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=0
=" > "${id}.primer.tmp"
}

#PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/semiquant/bioinfomatics/primer3-2.4.0/src/primer3_config/

get_primers "${template}" "${id}" "$start" "$size_min" "$size_max"
primer3_core --format_output < "${id}.primer.tmp" | grep -A 2 "OLIGO" | sed "s/LEFT PRIMER/${id}_F1/g" | sed "s/RIGHT PRIMER/${id}_R2/g" | awk '$9="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"$9' > "${id}.primers.tsv"
primer_found=$(wc -l < "${id}.primers.tsv")
rm "${id}.primer.tmp"



if [ $primer_found -gt 0 ]
then
    start=$(expr $seq_cycles \* 2 + $start - 15)
    get_primers "${template}" "${id}" "$start" 30 150
    primer3_core --format_output < "${id}.primer.tmp" | grep -A 2 "OLIGO" | sed "s/LEFT PRIMER/${id}_F2/g" | sed "s/RIGHT PRIMER/${id}_R1/g" | awk '$9="GACTGGAGTTCAGACGTGTGCTCTTCCGATCT"$9' > "${id}_inner.primers.tsv"
    primer_found=$(wc -l < "${id}_inner.primers.tsv")
    rm "${id}.primer.tmp"

    if [ $primer_found -gt 0 ]
    then
        echo -e "OLIGO\tstart\tlen\ttm\tgc%\tany_th\t3'_th\thairpin\tseq" > "${id}_hissPCR_primers.tsv"
        tail -n +2 "${id}.primers.tsv" >> "${id}_hissPCR_primers.tsv"
        tail -n +2 "${id}_inner.primers.tsv" >> "${id}_hissPCR_primers.tsv"
        if [[ $calc_hairpin == "Y" ]]
        then
            tail -n 2 "${id}_hissPCR_primers.tsv" | read -A arr
            strt=$arr[2]
            tail -n 1 "${id}_hissPCR_primers.tsv" | read -A arr
            end=$arr[2]
            inner_amplicon=${template:$strt:$end}
            seqfold "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGACATCTCGTATGCCGTCTTCTGCTTG${inner_amplicon}CAAGCAGAAGACGGCATACGAGATGTCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC" --celcius 60 --dot-bracket > "${id}_hissPCR_innerHairpin.txt"
        fi
        rm "${id}.primers.tsv" "${id}_inner.primers.tsv"
    else
        echo "could not find inner priemr set,for $id please relax conditions"
    fi
else 
    echo "could not find inner priemr set,for $id please relax conditions"
fi






