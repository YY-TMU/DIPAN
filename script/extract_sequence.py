import re
import os
import sys
import math
import HTSeq
import logging
import argparse
import pandas as pd
from collections import defaultdict
os.environ['NUMEXPR_MAX_THREADS'] = '4'

# set log
logging.basicConfig(format='%(asctime)s: %(message)s',
                    datefmt = '%Y-%m-%d %H:%M')
log = logging.getLogger()
log.setLevel(logging.INFO)

def creatDir(novel_dir):
    if not os.path.exists(novel_dir):
        os.makedirs(novel_dir)
    return


def ivFromGtf(gtf_file_path):
    # read gtf
    gtf_lines=HTSeq.GFF_Reader(gtf_file_path)
    iv_list = []
    novel_transcript_dic = dict()
    for gtf_line in gtf_lines:
        chrom = gtf_line.iv.chrom
        if not re.search("[M|_]",chrom):
            transcript_id = gtf_line.attr['transcript_id'] 
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            gene_id = gtf_line.attr['gene_id']
            gene_id = ":".join([chrom,gene_id]) 
            gene_name = gtf_line.attr['gene_name']
            start = gtf_line.iv.start + 1
            end = gtf_line.iv.end
            strand = gtf_line.iv.strand
            iv_type = gtf_line.type
            frame = gtf_line.frame
            iv_list.append([chrom,start,end,gene_id,transcript_id,strand,iv_type,gene_name,frame])
            if (iv_type == "transcript") and (transcript_id.startswith("IPA")):
                transcript_id_long = gtf_line.attr['transcript_id_long'] 
                novel_transcript_dic[transcript_id_long] = transcript_id
    # dataframe
    columns = ["chrom","start","end","gene_id","transcript_id","strand","type","gene_name","frame"]
    iv_df = pd.DataFrame.from_records(iv_list)
    iv_df.columns = columns
    return iv_df,novel_transcript_dic


def skipTe(skip_transcripts,transcript_dic,iv_df):
    skip_te_aa_dic = dict()
    for skip_transcript in skip_transcripts:
        skip_te_st,skip_te_end = map(int,skip_transcript.split("|")[3].split("_"))
        skip_transcript = transcript_dic[skip_transcript]
        novel_transcript_df = iv_df.loc[iv_df["transcript_id"]==skip_transcript,]
        strand = list(novel_transcript_df["strand"])[0]
        novel_cds_df = novel_transcript_df.loc[novel_transcript_df["type"]=="CDS",]
        if (sum(novel_cds_df["end"] - novel_cds_df["start"]) + novel_cds_df.shape[0]) % 3 == 0:
            cds_coornidates = list(zip(novel_cds_df["start"],novel_cds_df["end"]))
            cds_coornidates = sorted(cds_coornidates,key = lambda x:x[0])
            if strand == "+":
                last_cds_st,last_cds_end = cds_coornidates[-1]
                if (last_cds_st >= skip_te_st) and (last_cds_end >= skip_te_st):
                    te_exon_aa_number = math.ceil((last_cds_end-skip_te_st+1)/3)
                    skip_te_aa_dic[skip_transcript] = te_exon_aa_number
            elif strand == "-":
                last_cds_st,last_cds_end = cds_coornidates[0]
                if (last_cds_end <= skip_te_end) and (skip_te_end >= last_cds_st):
                    te_exon_aa_number = math.ceil((skip_te_end-last_cds_st+1)/3)
                    skip_te_aa_dic[skip_transcript] = te_exon_aa_number
    return skip_te_aa_dic


def compositeTe(composite_transcripts,transcript_dic,iv_df):
    # processing composite transcripts 
    composite_te_aa_dic = dict()
    for composite_transcript in composite_transcripts:
        ref_transcript = composite_transcript.split("|")[1]
        ref_transcript_df = iv_df.loc[iv_df["transcript_id"]==ref_transcript,]
        ref_exon_df = ref_transcript_df.loc[ref_transcript_df["type"]=="exon",]    
        ref_exon_coordinates = list(zip(ref_exon_df["start"],ref_exon_df["end"]))
        composite_te_st,composite_te_end = map(int,composite_transcript.split("|")[3].split("_"))
        composite_transcript = transcript_dic[composite_transcript]
        novel_transcript_df = iv_df.loc[iv_df["transcript_id"]==composite_transcript,]
        strand = list(novel_transcript_df["strand"])[0]
        novel_cds_df = novel_transcript_df.loc[novel_transcript_df["type"]=="CDS",]
        if (sum(novel_cds_df["end"] - novel_cds_df["start"]) + novel_cds_df.shape[0]) % 3 == 0:
            cds_coornidates = list(zip(novel_cds_df["start"],novel_cds_df["end"]))
            cds_coornidates = sorted(cds_coornidates,key = lambda x:x[0]) 
            if strand == "+":
                ref_exon_ls = [x for x in ref_exon_coordinates if x[0]==composite_te_st]
                assert len(ref_exon_ls) == 1
                ref_exon = ref_exon_ls[0]
                last_cds_st,last_cds_end = cds_coornidates[-1]
                if (last_cds_st >= composite_te_st) and (last_cds_end >= ref_exon[1] + 1):
                    te_exon_aa_number = math.ceil((last_cds_end - ref_exon[1])/3)
                    composite_te_aa_dic[composite_transcript] = te_exon_aa_number
            elif strand == "-":
                ref_exon_ls = [x for x in ref_exon_coordinates if x[1]==composite_te_end]
                assert len(ref_exon_ls) == 1
                ref_exon = ref_exon_ls[0]
                last_cds_st,last_cds_end = cds_coornidates[0]
                if (last_cds_end <= composite_te_end) and (last_cds_st <= ref_exon[0]-1):
                    te_exon_aa_number = math.ceil((ref_exon[0] - last_cds_st)/3)
                    composite_te_aa_dic[composite_transcript] = te_exon_aa_number
    return composite_te_aa_dic


def extractSeq(cds_fa,te_aa_dic):
    # processing fasta file
    transcript_aa_seq_dic = {}
    with open(cds_fa,"r") as fa:
        for x in fa.readlines():
            if x.startswith(">"):
                transcript = x.strip()
            else:
                transcript_aa_seq_dic.setdefault(transcript,[]).append(x.strip())
    name_trans = list(transcript_aa_seq_dic.keys())
    name_trans = {x.replace(">","").split(" ")[0]:x for x in name_trans if x.startswith(">IPA")}            
    transcript_aa_seq_dic =  {k.replace(">","").split(" ")[0]:"".join(v) 
                              for k,v in transcript_aa_seq_dic.items()
                              if k.startswith(">IPA")}
    # sequences for kmer
    kmer_seq_dic = defaultdict(dict)
    for kmer in range(8,12):
        for novel_transcript in te_aa_dic.keys():
            novel_transcript_fa = transcript_aa_seq_dic[novel_transcript]
            need_len = te_aa_dic[novel_transcript] + kmer - 1
            kmer_seq_dic[kmer][novel_transcript] = novel_transcript_fa[-need_len:]
    return kmer_seq_dic,name_trans


def writeCnt(fname,ls_cnt):
    with open(fname,"w") as wt:
        wt.write("\n".join(ls_cnt) + "\n")
    return


def create_parser(name):
    p = argparse.ArgumentParser(prog=name,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                description='Extract sequences.')
    # input
    g = p.add_argument_group('input')
    g.add_argument(
        '--gtf_pt',
        type = str,
        help='GTF file including novel isoforms.')
    g.add_argument(
        '--genome_fa_pt',
        type = str,
        help='Genome sequence in fasta format.')
    # output
    g = p.add_argument_group('output')
    g.add_argument(
        '--save_dir',
        type = str,
        help='Output directory.')
    return p

# run
# parse args
args = sys.argv
parser = create_parser(args[0])
args = parser.parse_args(args[1:])
# variable
merge_gtf_pt = args.gtf_pt
genome_fa_pt = args.genome_fa_pt
save_dir = args.save_dir
sample = ".".join(merge_gtf_pt.split("/")[-1].split(".")[:-1])
cds_fa = "".join([sample,".CDS.AA.fa"])
cds_fa = os.path.join(save_dir,cds_fa)
# make directory
creatDir(save_dir)
# extract cds sequence from GTF
cmd = "gffread %s -g %s -y %s" % (merge_gtf_pt,genome_fa_pt,cds_fa)
os.system(cmd)
log.info("Save %s" % cds_fa)
# extract transcripts
iv_df,novel_transcripts = ivFromGtf(merge_gtf_pt)
skip_transcripts = [x for x in novel_transcripts if ("Skip" in x) and ("novel_XM" in x)]
skip_te_aa_dic = skipTe(skip_transcripts,novel_transcripts,iv_df)
composite_transcripts = [x for x in novel_transcripts if ("Composite" in x) and ("novel_XM" in x)]
composite_te_aa_dic = compositeTe(composite_transcripts,novel_transcripts,iv_df)
te_aa_dic = {**skip_te_aa_dic,**composite_te_aa_dic}
# extend seq
kmer_seq_dic,name_trans = extractSeq(cds_fa,te_aa_dic)
for kmer,transcript_fa in kmer_seq_dic.items():
    kmer_ls = []
    for transcript,fa in transcript_fa.items():
        # processing transcript
        kmer_ls.append(name_trans[transcript])
        # process fa
        kmer_ls.extend(re.findall("[A-Za-z]{1,60}",fa))
    fname = ".".join(["extend",str(kmer),"fa"])
    fname = os.path.join(save_dir,fname)
    writeCnt(fname,kmer_ls)
    log.info("Save %s" % fname)

