# DIPAN

## Description
DIPAN is designed to identify neoantigens originating from intronic polyadenylation (IPA) events detected in tumor transcriptomes. These IPA-derived neoantigens have the potential to be presented by the MHC I molecules.

## Installation
DIPAN incorporates Python and Bash scripts. To set up the required environment, use the provided `environment.yaml` file to create a conda virtual environment.

1. Clone the repository:
```
git clone https://github.com/YY-TMU/DIPAN.git
```

2. Create the environment:
```
conda env create -f environment.yaml
conda activate DIPAN
```

3. netMHCpan and OptiType

netMHCpan can only be acquired through the official website and it should be added to the environment variable.
- [netMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/)

OptiType relies on Python 2.7. Due to compatibility issues with other scripts, it should be manually installed according to the instrument.
- [OptiType](https://github.com/FRED-2/OptiType)

## DIPAN options
The following options are avaliable:
- --annotated_finder/-a: Annotation file containing intron and exon information from IPAFinder.
- --bam_hla_input/-b: File containing BAM file paths and related HLA typing information.
- --bam_fq_input/-f: File containing BAM file paths and related FASTQ paths.
- --normal_proteome/-n: Amino acid sequences in the normal proteome.
- --annotated_gtf/-g: GTF file with annotated transcripts (RefSeq).
- --genome_fasta/-G: Genome in fasta format.
- --rank_threshold/-t: Set the rank threshold for high binding peptides. [default=2]
- --matched_normal/-m: Specify whether there are matched normal samples. If False, IPA-derived peptides found in normal samples are used. [default=True]
- --output_dir/-o: Output directory.
- --optitype_script: The path to OptiTypePipeline.py.
- --optitype_config: Configuration file of OptiType.


## Usage
### 1.Input
We provided the annotation file of IPAFinder for GRCh38, along with normal proteome amino acid sequences in `annotated_file` directory, and we recommend that users utilize it directly. The annotation file for GRCh38 of RefSeq can be downloaded from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/archive/hg38/ncbiRefSeq/109.20211119/hg38.109.20211119.ncbiRefSeq.gtf.gz). It should be noted that the annotation file of IPAFinder must match the GTF file; if the GTF file is changed, the annotation file should be adjusted accordingly for [IPAFinder](https://github.com/ZhaozzReal/IPAFinder).

DIPAN offers two options for users. If HLA-I typing information is unavailable, users should provide `bam_fq_input`, `optitype_script` and `optitype_config`. OptiType is used to calculate HLA typing. `bam_fq_input` should include paths to BAM files and their corresponding FASTQ files. Alternatively, if HLA typing is already known, users should provide `bam_hla_input`, which includes paths to BAM files along with their associated HLA typing information.

DIPAN could be tested using [recommended files](https://zenodo.org/records/10970002).


### 2. Unknown HLA typing
```
DIPAN.sh -a <IPAFinder_anno.txt> -f <bam_fq_input.txt> -n <Normal_proteome.fa> -g <refseq.gtf> -G <genome.fa> -o <output_dir> -optitype_script <OptiTypePipeline.py> -optitype_config <optitype.config>
```
bam_fq_input.txt contains paths of BAM file and related FASTQ file, as shown below:
```
Tumor  /path/tumor.sorted.bam  /path/tumor_1.fq,/path/tumor_2.fq
Normal /path/normal.sorted.bam   
```

### 3. Known HLA typing
```
DIPAN.sh -a <IPAFinder_anno.txt> -b <bam_hla_input.txt> -n <Normal_proteome.fa> -g <refseq.gtf> -G <genome.fa> -o <output_dir>
```
bam_hla_input.txt contains paths of BAM file and related HLA typing information, as shown below:
```
Tumor  /path/tumor.sorted.bam  HLA-A*01:01,HLA-B*44:02,HLA-C*06:02
Normal /path/normal.sorted.bam   
```
We collected 730 normal samples from TCGA and curated a list of IPA-derived peptides found in these normal samples. This list, along with the normal human proteome provided, can serve as a control when normal RNA-seq data are unavailable. When matched samples are missing, set matched_normal to False.

### 4. Output
The final output includes filtered neoantigens from tumor samples, and the output consists of the following columns:
Column | Description
------ | -----------
SYMBOL | Gene symbol
Terminal_exon | Genomic location of corresponding terminal exon of IPA isoform
IPAtype | Type of terminal exon (Skipped or Composite)
IPUI | Abundance of IPA events
HLA | HLA-I typing
Peptide | Amino acid sequence of the potential IPA-derived neoantigens
%Rank | Rank of the predicted binding score compared to a set of random natural peptides. This measure is not affected by inherent bias of certain molecules towards higher or lower mean predicted affinities.

