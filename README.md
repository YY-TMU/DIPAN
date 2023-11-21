# DIPAN

## Description
DIPAN searches for neoantigens derived from intronic polyadenylation events found in tumor transcriptomes. We show that the neoantigens could be presented by the MHC I complex using mass spectrometry data.

## Installation
This pipeline includes both Python and Bash scripts. The installation depends on the environment, taking about 10 minutes.

1.Install the following programs before running this pipeline, and the required programs could be found in `requirements.txt` file. 
- python (3.6+, required packages HTSeq, numpy, pandas, tqdm, pyfasta, scipy)
- [IPAFinder](https://github.com/ZhaozzReal/IPAFinder)
- [netMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/)
- [OptiType](https://github.com/FRED-2/OptiType)

2.Clone the repository and change directory.
```
git clone https://github.com/YY-TMU/DIPAN.git
cd DIPAN
```

## Usage
### Identify neoantigens
In the following link, genome file for GRCh38 of RefSeq could be downloaded.
- [Human (hg38)](https://hgdownload.soe.ucsc.edu/goldenPath/archive/hg38/ncbiRefSeq/109.20211119/hg38.109.20211119.ncbiRefSeq.gtf.gz)

The following options are all required:
- --annotated_finder/-a: Annotation file containing intron and exon information.
- --bam_files/-b: File containing BAM file paths between normal and tumor samples.
- --fq_files/-f: File containing FASTQ file paths in tumor samples.
- --normal_proteome/-n: Amino acid sequences in the normal proteome.
- --annotated_gtf/-g: GTF file with annotated transcripts.
- --genome_fasta/-G: Genome in fasta format.
- --output_dir/-o: Output directory.
- --optitype_script: The path to OptiTypePipeline.py.
- --optitype_config: Configuration file of OptiType.

**Command**
```
DIPAN.sh -a <IPAFinder_anno.txt> -b <bamfiles.txt> -f <fqfiles.txt> -n <normal_proteome> -g <refseq.gtf> -G <genome.fa> -o <output_dir> -optitype_script <OptiTypePipeline.py> -optitype_config <optitype.config>
```
We present the annotation file of IPAFinder for hg38 and normal proteome amino acid sequences sourced from UniProt in `annotated_file` directory, and we suggest users utilize it directly.

bamfiles.txt contains paths between two conditions, as shown below:
```
condition1=/path/tumor1.bam,/path/tumor2.bam 
condition2=/path/normal1.bam,/path/normal2.bam
```

fqfiles.txt contains FASTQ file paths in tumor samples, as shown below:
```
tumor1=/path/tumor1_1.fastq, /path/tumor1_2.fastq
tumor2=/path/tumor2_1.fastq, /path/tumor2_2.fastq
```

The final output includes filtered neoantigens from tumor samples, and the output consists of the following columns:
- Pos: Residue number (starting from 0) of the peptide in the protein sequence.
- HLA: Specified MHC molecule / Allele name.
- Peptide: Amino acid sequence of the potential ligand.
- Core: The minimal 9 amino acid binding core directly in contact with the MHC.
- Of: The starting position of the Core within the Peptide (if > 0, the method predicts a N-terminal protrusion).
- Gp: Position of the deletion, if any.
- Gl: Length of the deletion, if any.
- Ip: Position of the insertion, if any.
- Il: Length of the insertion, if any.
- Icore: Interaction core. This is the sequence of the binding core including eventual insertions of deletions.
- Identity: Protein identifier, i.e. the name of the FASTA entry.
- Score: The raw prediction score.
- %Rank: Rank of the predicted binding score compared to a set of random natural peptides. This measure is not affected by inherent bias of certain molecules towards higher or lower mean predicted affinities. Strong binders are defined as having %rank<0.5, and weak binders with %rank<2. We advise to select candidate binders based on %Rank rather than Score
- BindLevel: (SB: Strong Binder, WB: Weak Binder). The peptide will be identified as a strong binder if the %Rank is below the specified threshold for the strong binders (by default, 0.5%). The peptide will be identified as a weak binder if the %Rank is above the threshold of the strong binders but below the specified threshold for the weak binders (by default, 2%).


