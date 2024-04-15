#!/bin/bash

set -e
print_help() {
  echo "Summary: Recognize neoantigens from IPA events in tumor."
  echo "Usage:   $0 -a <IPAFinder_anno.txt> -f <bam_fq_input.txt> -n <Normal_proteome.fa> -g <refseq.gtf> -G <genome.fa> -o <output_dir> -optitype_script <OptiTypePipeline.py> -optitype_config <optitype.config>"
  echo
  echo "Options:"
  printf "%-35s %-s\n"  "-a, -annotated_finder"  "Annotation file containing intron and exon information from IPAFinder."
  printf "%-35s %-s\n"  "-b, -bam_hla_input"  "File containing BAM file paths and related HLA typing information."
  printf "%-35s %-s\n"  "-f, -bam_fq_input"  "File containing BAM file paths and related FASTQ paths."
  printf "%-35s %-s\n"  "-n, -normal_proteome"  "Amino acid sequences in the normal proteome."
  printf "%-35s %-s\n"  "-g, -annotated_gtf"  "GTF file with annotated transcripts (RefSeq)."
  printf "%-35s %-s\n"  "-G, -genome_fasta"  "Genome in fasta format."
  printf "%-35s %-s\n"  "-t, --rank_threshold"  "Set the rank threshold for high binding peptides. [default=2]"
  printf "%-35s %-s\n"  "-m, --matched_normal"  "Specify whether there are matched normal samples. If False, IPA-derived peptides commonly found in normal samples are used. [default=True]"
  printf "%-35s %-s\n"  "-o, -output_dir"  "Output directory."
  printf "%-35s %-s\n"  "-optitype_script"  "The path to OptiTypePipeline.py."
  printf "%-35s %-s\n"  "-optitype_config"  "Configuration file of OptiType."
  printf "%-35s %-s\n"  "-h, -help"  "Print this help menu."
}

args_count=$#
if [ $args_count -eq 0 ]; then
  echo "No arguments provided. Exiting..."
  print_help
  exit 1
fi

ARGS=`getopt -a -o a:b:f:n:g:G:t:m:o:h  --long annotated_finder:,bam_hla_input:,bam_fq_input:,normal_proteome:,annotated_gtf:,genome_fasta:,rank_threshold:,matched_normal:,output_dir:,optitype_script:,optitype_config:,help -n "$0" -- "$@"`
[ $? -ne 0 ] && exit 1
echo ARGS=[$ARGS]
eval set -- "${ARGS}"


# set default
rank_threshold=2
matched_normal=True
save_dir=DIPAN_result

# set parameter
while true
do
    case $1 in
        -a|--annotated_finder)
          IPAFinder_anno=$2
          shift 2
          ;;
        -b|--bam_hla_input)
          bam_hla_input=$2
          shift 2
          ;;        
        -f|--bam_fq_input)
          bam_fq_input=$2
          shift 2
          ;;
        -n|--normal_proteome)
          normal_proteome_pt=$2
          shift 2
          ;;
        -g|--annotated_gtf)
          annotated_gtf_pt=$2
          shift 2
          ;;
        -G|--genome_fasta)
          genome_fasta_pt=$2
          shift 2
          ;;
        -t|--rank_threshold)
          rank_threshold=$2
          shift 2
          ;;
        -m|--matched_normal)
          [[ $2 == True ]] || [[ $2 == False ]] || { echo Error: matched_normal must be Ture or False! ; exit 1 ; }
          matched_normal=$2
          shift 2
          ;;
        -o|--output_dir)
          save_dir=$2
          shift 2
          ;;
        --optitype_script)
          optitype_script_pt=$2
          shift 2
          ;;
        --optitype_config)
          optitype_config_pt=$2
          shift 2
          ;;
        -h|--help)
          print_help
          exit 1
          ;;
        --)
          shift
          break
          ;;
        *)
          echo "Internal error!"
          exit 1
          ;;
    esac
done   

# check input
[ -f "$IPAFinder_anno" ] || { echo "Error: Invalid annotated_finder!"; exit 1; }
[ -f "$normal_proteome_pt" ] || { echo "Error: Invalid normal_proteome!"; exit 1; }
[ -f "$annotated_gtf_pt" ] || { echo "Error: Invalid annotated_gtf!"; exit 1; }
[ -f "$genome_fasta_pt" ] || { echo "Error: Invalid genome_fasta!"; exit 1; }
if [ -f "$bam_fq_input" ]  && [ -f "$optitype_script_pt" ] && [ -f "$optitype_config_pt" ] && [ -z "$bam_hla_input" ] ; then
    echo "OptiType will be used to calculate HLA typing."
    [ -x $optitype_script_pt ] || { echo "$optitype_script_pt is not executable!"; exit 1; }
elif [ -f "$bam_hla_input" ] && [ -z "$bam_fq_input" ]  && [ -z "$optitype_script_pt" ] && [ -z "$optitype_config_pt" ]; then
    echo "Custom HLA typing."
else
    echo "Error: Please provide either (FASTQ paths + Optitype paths) or HLA class I typing information, but not both!"
    exit 1
fi


# set default
process=5
# Run
[[ -d $save_dir ]] || mkdir $save_dir

########## Detect IPA ###########
echo -e `date '+%Y.%m.%d %H:%M'` Detecting IPA ... 
if [[ -f "$bam_hla_input" ]]; then
    # tumor samples
    ipa_tumor_dir=$save_dir/1.IPA_result/tumor
    [[ -d $ipa_tumor_dir ]] || mkdir -p $ipa_tumor_dir
    tumor_ipa_pts=""
    declare -A hla_dict
    while IFS=$'\t' read -r tp pt hla_types ; do
        sample=${pt##*/}
        sample=${sample%%.*}
        echo -e `date '+%Y.%m.%d %H:%M'` $sample
        sample_dir=${ipa_tumor_dir}/${sample}
        [[ -d $sample_dir ]] || mkdir $sample_dir
        tmp_pt=${sample_dir}/tmp.txt
        echo condition1=$pt > $tmp_pt
        python IPAFinder_DetectIPA.py -b $tmp_pt -anno $IPAFinder_anno -p $process -o ${sample_dir}/IPAFinder_IPUI.txt
        tumor_ipa_pts="${tumor_ipa_pts} ${sample_dir}/IPAFinder_IPUI.txt"  
        rm $tmp_pt
        hla_types=$(echo $hla_types | sed 's/*//g')
        hla_dict[$sample]=${hla_types}
    done < <(grep "Tumor" $bam_hla_input)

    # normal samples
    if [[ $matched_normal == True ]] ; then
        ipa_normal_dir=$save_dir/1.IPA_result/normal
        [[ -d $ipa_normal_dir ]] || mkdir -p $ipa_normal_dir
        normal_ipa_pts=""
        while IFS=$'\t' read -r tp pt ; do
            sample=${pt##*/}
            sample=${sample%%.*}
            echo -e `date '+%Y.%m.%d %H:%M'` $sample
            sample_dir=${ipa_normal_dir}/${sample}
            [[ -d $sample_dir ]] || mkdir $sample_dir
            tmp_pt=${sample_dir}/tmp.txt
            echo condition2=$pt > $tmp_pt
            python IPAFinder_DetectIPA.py -b $tmp_pt -anno $IPAFinder_anno -p $process -o ${sample_dir}/IPAFinder_IPUI.txt
            normal_ipa_pts="${normal_ipa_pts} ${sample_dir}/IPAFinder_IPUI.txt"
            rm $tmp_pt
        done < <(grep "Normal" $bam_hla_input)
     else
         echo "There are No matched nromal sample."
     fi
elif [[ -f $bam_fq_input ]]; then
    # tumor samples
    ipa_tumor_dir=$save_dir/1.IPA_result/tumor
    [[ -d $ipa_tumor_dir ]] || mkdir -p $ipa_tumor_dir
    tumor_ipa_pts=""
    declare -A fq_dict
    while IFS=$'\t' read -r tp pt fq_pts ; do
        sample=${pt##*/}
        sample=${sample%%.*}
        echo -e `date '+%Y.%m.%d %H:%M'` $sample
        sample_dir=${ipa_tumor_dir}/${sample}
        [[ -d $sample_dir ]] || mkdir $sample_dir
        tmp_pt=${sample_dir}/tmp.txt
        echo condition1=$pt > $tmp_pt
        python IPAFinder_DetectIPA.py -b $tmp_pt -anno $IPAFinder_anno -p $process -o ${sample_dir}/IPAFinder_IPUI.txt
        tumor_ipa_pts="${tumor_ipa_pts} ${sample_dir}/IPAFinder_IPUI.txt"  
        rm $tmp_pt
        fq_dict[$sample]=${fq_pts}
    done < <(grep "Tumor" $bam_fq_input)

    # normal samples
    if [[ $matched_normal == True ]] ; then
        ipa_normal_dir=$save_dir/1.IPA_result/normal
        [[ -d $ipa_normal_dir ]] || mkdir -p $ipa_normal_dir
        normal_ipa_pts=""
        while IFS=$'\t' read -r tp pt ; do
            sample=${pt##*/}
            sample=${sample%%.*}
            echo -e `date '+%Y.%m.%d %H:%M'` $sample
            sample_dir=${ipa_normal_dir}/${sample}
            [[ -d $sample_dir ]] || mkdir $sample_dir
            tmp_pt=${sample_dir}/tmp.txt
            echo condition2=$pt > $tmp_pt
            python IPAFinder_DetectIPA.py -b $tmp_pt -anno $IPAFinder_anno -p $process -o ${sample_dir}/IPAFinder_IPUI.txt
            normal_ipa_pts="${normal_ipa_pts} ${sample_dir}/IPAFinder_IPUI.txt"
            rm $tmp_pt
        done < <(grep "Normal" $bam_fq_input)
    else
        echo "There are No matched nromal sample."
    fi
fi
########### Detect IPA ###########


########### Filter IPA ###########
echo -e `date '+%Y.%m.%d %H:%M'` Filtering IPA ... 
filter_ipa_dir=$save_dir/2.IPA_filter
[[ -d $filter_ipa_dir ]] || mkdir $filter_ipa_dir
# merge normal
if [[ $matched_normal == True ]] ; then
    merge_normal_ipa_pt=${filter_ipa_dir}/normal.ipa.txt
    awk -v OFS="_" 'FNR>1 {print $1,$2}' $normal_ipa_pts | sort | uniq > $merge_normal_ipa_pt
    regex=""
    for ipa in `cat $merge_normal_ipa_pt`
    do
        regex+="${ipa}$\|"
    done
    regex=$(echo $regex | sed 's/..$//')
else
    regex="^$"
fi
# tumor samples
filter_tumor_ipa_pts=""
for pt in $tumor_ipa_pts
do
    sample=${pt%/*}
    sample=${sample##*/}
    sample_dir=${filter_ipa_dir}/${sample}
    [[ -d $sample_dir ]] || mkdir $sample_dir
    result_pt=${sample_dir}/IPAFinder_IPUI.txt
    awk '{print $0,$1"_"$2}' $pt | grep -v $regex | awk -v OFS="\t" '{print $1,$2,$3,$4,$5}' > $result_pt
    filter_tumor_ipa_pts="${filter_tumor_ipa_pts} ${result_pt}"
done
########### Filter IPA ###########


########### Assembel transcript and extract sequences ###########
echo -e `date '+%Y.%m.%d %H:%M'` Assembel transcript and extract sequences ... 
gtf_dir=$save_dir/3.IPA_gtf
[[ -d $gtf_dir ]] || mkdir $gtf_dir
for pt in $filter_tumor_ipa_pts
do
    sample=${pt%/*}
    sample=${sample##*/}
    echo -e `date '+%Y.%m.%d %H:%M'` Processing $sample ...
    sample_dir=${gtf_dir}/${sample}
    [[ -d $sample_dir ]] || mkdir $sample_dir
    gtf_file=${sample_dir}/${sample}.gtf
    # assemble transcript
    python build_transcript.py \
        --finder_pa_file $pt \
        --annotated_gtf $annotated_gtf_pt \
        --fa_path $genome_fasta_pt \
        --save_file $gtf_file
    # extract sequence
    python extract_sequence.py \
        --gtf_pt $gtf_file \
        --genome_fa_pt $genome_fasta_pt \
        --save_dir $sample_dir
done
########### Assembel transcript and extract sequences ###########


########### HLA typing and netMHC ###########
echo -e `date '+%Y.%m.%d %H:%M'` HLA typing and netMHC ... 
mhc_dir=$save_dir/4.IPA_netMHC
[[ -d $mhc_dir ]] || mkdir $mhc_dir

if [[ -f "$bam_fq_input" ]]; then
    # HLA typing
    declare -A hla_dict
    for tumor_sample in `ls $gtf_dir`
    do
        echo -e `date '+%Y.%m.%d %H:%M'` Processing $tumor_sample ...
        sample_dir=${mhc_dir}/${tumor_sample}
        [[ -d $sample_dir ]] || mkdir $sample_dir
        fq_file_info=${fq_dict[$tumor_sample]}
        fq_file_number=$(echo $fq_file_info | awk -F "," '{print NF}')
        if [ $fq_file_number -eq 2 ]
        then
            fq1_file=$(echo $fq_file_info | awk -F "," '{print $1}')
            [ -f "$fq1_file" ] || { echo "Error: Invalid fastq path!"; exit 1; }
            fq2_file=$(echo $fq_file_info | awk -F "," '{print $2}')
            [ -f "$fq2_file" ] || { echo "Error: Invalid fastq path!"; exit 1; }
            $optitype_script_pt --input $fq1_file $fq2_file -v --rna -p $tumor_sample --outdir $sample_dir -c $optitype_config_pt
        elif [ $fq_file_number -eq 1 ]
        then
            fq_file=$(echo $fq_file_info | awk -F "," '{print $1}')
            [ -f "$fq_file" ] || { echo "Error: Invalid fastq path!"; exit 1; }
            $optitype_script_pt --input $fq_file -v --rna -p $tumor_sample --outdir $sample_dir -c $optitype_config_pt
        fi    
        type_result_pt=${sample_dir}/${tumor_sample}_result.tsv
        hla_types=$(awk -v OFS="\n" 'NR>1 {print $2,$3,$4,$5,$6,$7}' $type_result_pt | sort | uniq | sed 's/*//g' | awk '{print "HLA-"$0}' | tr "\n" "," | sed 's/.$//')
        hla_dict[$tumor_sample]=$hla_types
    done
fi

# netMHC
for tumor_sample in `ls $gtf_dir`
do
    echo -e `date '+%Y.%m.%d %H:%M'` Processing $tumor_sample
    sample_dir=${mhc_dir}/${tumor_sample}
    [[ -d $sample_dir ]] || mkdir $sample_dir
    for kmer in {8..11};do
        kmer_fa_pt=${gtf_dir}/$tumor_sample/extend.${kmer}.fa
        hla_type_info=${hla_dict[$tumor_sample]}
        fname=${sample_dir}/netMHC.${kmer}.result.txt
        netMHCpan -a $hla_type_info -f $kmer_fa_pt -l $kmer > ${fname}
        echo -e `date '+%Y.%m.%d %H:%M'` Save $fname
    done   
done

# filter
for tumor_sample in `ls $mhc_dir`
do
    tumor_sample_result_dir=${mhc_dir}/${tumor_sample}
    filter_result_file=${save_dir}/${tumor_sample}.neoantigen.txt
    netMHC_result_pts=$(find $tumor_sample_result_dir -name "netMHC.*txt")
    ipa_info=${save_dir}/3.IPA_gtf/${tumor_sample}/${tumor_sample}.IPA_info.txt
    # header
    awk '/BindLevel/' $netMHC_result_pts | awk 'NR==1 {print}' > $filter_result_file
    # filter peptides
    awk '!/#/ && $2~/HLA-/' $netMHC_result_pts | awk -v rt=$rank_threshold '{if($13<rt) print $0}' | awk '{gsub(/<= *SB/, "<=SB"); gsub(/<= *WB/, "<=WB"); print}' >> $filter_result_file
    # filter with normal proteome
    python filter_normal_proteome.py \
        --normal_proteome_fa $normal_proteome_pt \
        --peptide_file $filter_result_file \
        --ipa_annotation_file $ipa_info \
        --matched_normal $matched_normal \
        --save_file $filter_result_file
done

# delete intermediate file
rm -r $save_dir/1.IPA_result
rm -r $save_dir/2.IPA_filter
rm -r $save_dir/3.IPA_gtf
rm -r $save_dir/4.IPA_netMHC
