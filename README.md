# DIPAN

## Description
DIPAN searches for neoantigens derived from intronic polyadenylation events found in tumor transcriptomes. We show that the neoantigens could be presented by the MHC I complex using mass spectrometry data.

## Installation
This pipeline includes both Python and Bash scripts. A conda virtual environment can be created using the provided `environment.yaml` file.

1. Clone the repository:
```
git clone https://github.com/YY-TMU/DIPAN.git
```

2. Create the environment:
```
conda env create -f environment.yaml
conda activate DIPAN
```

3.Due to permission and compatibility issues, the following software needs to be installed manually and added to the PATH environment variable.
- [netMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/)
- [OptiType](https://github.com/FRED-2/OptiType)



## Usage
### Identify neoantigens
In the following link, genome file for GRCh38 of RefSeq could be downloaded.
- [Human (hg38)](https://hgdownload.soe.ucsc.edu/goldenPath/archive/hg38/ncbiRefSeq/109.20211119/hg38.109.20211119.ncbiRefSeq.gtf.gz)

In the following link, test file could be downloaded.
- [Test data](https://zenodo.org/records/10970002)

The following options are all required:
- --annotated_finder/-a: Annotation file containing intron and exon information.
- --bam_hla_input/-b: File containing BAM file paths and related HLA typing information.
- --bam_fq_input/-f: File containing BAM file paths and related FASTQ paths.
- --normal_proteome/-n: Amino acid sequences in the normal proteome.
- --annotated_gtf/-g: GTF file with annotated transcripts.
- --genome_fasta/-G: Genome in fasta format.
- --rank_threshold/-t: Set the rank threshold for high binding peptides.
- --matched_normal/-m: Specify whether there are matched normal samples. If False, IPA-derived peptides commonly found in TCGA normal samples are used.
- --output_dir/-o: Output directory.
- --optitype_script: The path to OptiTypePipeline.py.
- --optitype_config: Configuration file of OptiType.

##Usage
We present the annotation file of IPAFinder for hg38 and normal proteome amino acid sequences sourced from UniProt in `annotated_file` directory, and we suggest users utilize it directly.
### 1. HLA typing is calculated by OptiType.
```
DIPAN.sh -a <IPAFinder_anno.txt> -f <bam_fq_input.txt> -n <Normal_proteome.fa> -g <refseq.gtf> -G <genome.fa> -o <output_dir> -optitype_script <OptiTypePipeline.py> -optitype_config <optitype.config>
```
bam_fq_input.txt contains paths of BAM file and related FASTQ file, as shown below:
```
Tumor  /path/tumor.sorted.bam  /path/tumor_1.fq,/path/tumor_2.fq
Normal /path/normal.sorted.bam   
```
### 2. HLA typing could be provided by the user.
```
DIPAN.sh -a <IPAFinder_anno.txt> -b <bam_hla_input.txt> -n <Normal_proteome.fa> -g <refseq.gtf> -G <genome.fa> -o <output_dir>
```
bam_hla_input.txt contains paths of BAM file and related HLA typing information, as shown below:
```
Tumor  /path/tumor.sorted.bam  HLA-A*01:01,HLA-B*44:02,HLA-C*06:02
Normal /path/normal.sorted.bam   
```
We collected 730 normal samples from TCGA and curated a list of IPA-derived peptides found in these normal samples. This list, along with the normal human proteome provided, can serve as a control when normal RNA-seq data are unavailable. When matched samples are missing, set matched_normal to False.

The final output includes filtered neoantigens from tumor samples, and the output consists of the following columns:
Column | Description
------ | -----------
SYMBOL | gene symbol
Terminal_exon | genomic location of corresponding terminal exon of IPA isoform
IPAtype | type of terminal exon (Skipped or Composite)
IPUI | abundance of IPA events
HLA | HLA-I typing
Peptide | Amino acid sequence of the potential IPA-derived neoantigens
Core | The minimal 9 amino acid binding core directly in contact with the MHC.
%Rank | Rank of the predicted binding score compared to a set of random natural peptides. This measure is not affected by inherent bias of certain molecules towards higher or lower mean predicted affinities.



