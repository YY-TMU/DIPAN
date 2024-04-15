#!/usr/bin/env python

import os
import re
import sys
import HTSeq
import logging
import argparse
import numpy as np
import pandas as pd
from pyfasta import Fasta
os.environ['NUMEXPR_MAX_THREADS'] = '4'

# set log
logging.basicConfig(format='%(asctime)s: %(message)s',
                    datefmt = '%Y-%m-%d %H:%M')
log = logging.getLogger()
log.setLevel(logging.INFO)

# function
def exonFromGtf(gtf_file_path,
                genes):
    # read gtf
    gtf_file=HTSeq.GFF_Reader(gtf_file_path)
    # exon list
    exon_list = []
    for gtf_line in gtf_file:
        if gtf_line.type == "exon":   
            gene_name = gtf_line.attr['gene_name']
            if gene_name in genes:
                chrom = gtf_line.iv.chrom
                if not re.search("[M|_]",chrom):
                    transcript_id = gtf_line.attr['transcript_id']
                    if re.match("NM_",transcript_id): 
                        if not chrom.startswith("chr"):
                            chrom = "chr" + chrom
                        start = gtf_line.iv.start
                        end = gtf_line.iv.end                    
                        strand = gtf_line.iv.strand
                        exon_list.append([chrom,start,end,gene_name,transcript_id,strand])        
    # exon df
    columns = ["chrom","start","end","gene_name","transcript_id","strand"]
    exon_df = pd.DataFrame.from_records(exon_list)
    exon_df.columns = columns
    return exon_df

def addID(ipa_infos):
    merge_ipa_infos = []
    assert len(ipa_infos) < 10 ** 7
    for idx,ipa_info in enumerate(ipa_infos):
        trans_idx = str(idx + 1)
        trans_idx = trans_idx.zfill(7)
        ipa_info.append(trans_idx)
        merge_ipa_infos.append(ipa_info)
    return merge_ipa_infos


def handleComposite(rt_pa_composite,exon_df):
    composite_infos = []
    group = "Composite"
    for x in rt_pa_composite.itertuples():
        _,gene_name,_,predict_te_region,_,IPUI = list(x)
        chrom,start,end = re.search("(chr.*?):(\d+)-(\d+)",predict_te_region).groups()
        gene_id = ":".join([chrom,gene_name])
        start,end = int(start),int(end)
        gene_exon_df = exon_df.loc[(exon_df["gene_name"]==gene_name) & (exon_df["chrom"]==chrom),]
        strand = list(set(gene_exon_df["strand"]))
        assert len(strand) == 1
        strand = strand[0]
        if strand == "+":
            template_exon_df = gene_exon_df.loc[gene_exon_df["end"]==start,]
            if not template_exon_df.empty:
                exon_coordinates = list(set(zip(template_exon_df["start"],template_exon_df["end"])))
                if len(exon_coordinates) > 1:
                    max_index = (template_exon_df["end"] - template_exon_df["start"]).idxmax()
                    template_exon_df = template_exon_df.loc[[max_index],]
                    exon_coordinates = list(set(zip(template_exon_df["start"],template_exon_df["end"])))
                template_transcript = list(set(template_exon_df["transcript_id"]))[0]
                # te region
                te_start = exon_coordinates[0][0]+1
                te_end = end
                te_region = ":".join(map(str,[chrom,te_start,te_end,strand]))
                # up exon
                find_up_exon_df = gene_exon_df.loc[gene_exon_df["transcript_id"]==template_transcript,]
                find_up_exon_df = find_up_exon_df.loc[find_up_exon_df["start"]<(te_start-1),]
                if find_up_exon_df.empty:
                    up_site = "first_exon"
                else:
                    up_exon_coordinates = list(find_up_exon_df["end"])
                    up_site = ":".join([chrom,str(max(up_exon_coordinates)),strand])
                composite_infos.append([te_region,up_site,gene_id,template_transcript,predict_te_region,group,IPUI])
            else:
                raise Exception("Check " + predict_te_region)
        elif strand == "-":
            template_exon_df = gene_exon_df.loc[gene_exon_df["start"]==end,]
            if not template_exon_df.empty:
                exon_coordinates = list(set(zip(template_exon_df["start"],template_exon_df["end"])))
                if len(exon_coordinates) > 1:
                    max_index = (template_exon_df["end"] - template_exon_df["start"]).idxmax()
                    template_exon_df = template_exon_df.loc[[max_index],]
                    exon_coordinates = list(set(zip(template_exon_df["start"],template_exon_df["end"])))
                template_transcript = list(set(template_exon_df["transcript_id"]))[0]
                # te region
                te_start = start + 1
                te_end = exon_coordinates[0][1]
                te_region = ":".join(map(str,[chrom,te_start,te_end,strand]))
                # up exon
                find_up_exon_df = gene_exon_df.loc[gene_exon_df["transcript_id"]==template_transcript,]
                find_up_exon_df = find_up_exon_df.loc[find_up_exon_df["start"]>te_end,]
                if find_up_exon_df.empty:
                    up_site = "first_exon"
                else:
                    up_exon_coordinates = list(find_up_exon_df["start"])
                    up_site = ":".join([chrom,str(min(up_exon_coordinates)+1),strand])
                composite_infos.append([te_region,up_site,gene_id,template_transcript,predict_te_region,group,IPUI])
            else:
                raise Exception("Check " + predict_te_region)
    composite_infos =  addID(composite_infos)   
    return composite_infos


def handleSkip(rt_pa_skip,exon_df):
    skip_infos = []
    group = "Skipped"
    template_transcript = "."
    for x in rt_pa_skip.itertuples():
        _,gene_name,_,predict_te_region,_,IPUI = list(x)
        chrom,start,end = re.search("(chr.*?):(\d+)-(\d+)",predict_te_region).groups()
        start,end = int(start),int(end)
        gene_id = ":".join([chrom,gene_name])
        gene_exon_df = exon_df.loc[(exon_df["gene_name"]==gene_name) & (exon_df["chrom"]==chrom),]
        strand = list(set(gene_exon_df["strand"]))
        assert len(strand) == 1
        strand = strand[0]
        te_region = ":".join(map(str,[chrom,start+1,end,strand]))
        if strand == "+":
            # up exon
            find_up_exon_df = gene_exon_df.loc[gene_exon_df["end"] < start,]
            up_exon_coordinates = list(find_up_exon_df["end"])
            up_site = ":".join([chrom,str(max(up_exon_coordinates)),strand])
        elif strand == "-":
            # up exon
            find_up_exon_df = gene_exon_df.loc[gene_exon_df["start"] > end,]
            up_exon_coordinates = list(find_up_exon_df["start"])
            up_site = ":".join([chrom,str(min(up_exon_coordinates)+1),strand])
        skip_infos.append([te_region,up_site,gene_id,template_transcript,predict_te_region,group,IPUI])
    skip_infos =  addID(skip_infos)
    return skip_infos


def extractUpInfo(finder_pa_file,gtf_path):
    rt_pa = pd.read_table(finder_pa_file)
    all_genes = list(set(rt_pa["SYMBOL"]))
    rt_pa_skip = rt_pa.loc[rt_pa["IPAtype"]=="Skipped",]
    rt_pa_composite = rt_pa.loc[rt_pa["IPAtype"]=="Composite",]
    exon_df = exonFromGtf(gtf_path,all_genes)
    composite_infos = handleComposite(rt_pa_composite,exon_df)
    skip_infos = handleSkip(rt_pa_skip,exon_df)
    total_infos = composite_infos + skip_infos
    total_info_df = pd.DataFrame(total_infos,
                                 columns = ["region1","up_site","gene_id","transcript_id","region2","group","IPUI","IPA_ID"])
    return total_info_df


def splitTerminal(terminal_exon_df):
    novel_terminal_genes = list(set(terminal_exon_df["gene_id"]))
    # split te
    rt_te_skip = terminal_exon_df.loc[terminal_exon_df["group"] == "Skipped",]
    rt_te_composite = terminal_exon_df.loc[terminal_exon_df["group"] == "Composite",]
    return novel_terminal_genes,rt_te_skip,rt_te_composite


def ivFromGtf(gtf_file_path,novel_terminal_genes):
    # read gtf
    gtf_lines=HTSeq.GFF_Reader(gtf_file_path)
    iv_list = []
    for gtf_line in gtf_lines:
        chrom = gtf_line.iv.chrom
        if not re.search("[M|_]",chrom):
            transcript_id = gtf_line.attr['transcript_id']
            if re.match("NM_",transcript_id): 
                if not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                gene_id = gtf_line.attr['gene_id']
                gene_id = ":".join([chrom,gene_id]) 
                gene_name = gtf_line.attr['gene_name']
                if gene_id in novel_terminal_genes:
                    start = gtf_line.iv.start + 1
                    end = gtf_line.iv.end
                    strand = gtf_line.iv.strand
                    iv_type = gtf_line.type
                    frame = gtf_line.frame
                    iv_list.append([chrom,start,end,gene_id,transcript_id,strand,iv_type,gene_name,frame])
    # dataframe
    columns = ["chrom","start","end","gene_id","transcript_id","strand","type","gene_name","frame"]
    iv_df = pd.DataFrame.from_records(iv_list)
    iv_df.columns = columns
    return iv_df


class novelTranscript:
    # object of transcript
    def __init__(self,
                 gene_id,
                 gene_name,
                 transcript_id,
                 exon,
                 start_codon,
                 stop_codon,
                 five_UTR,
                 strand,
                 annotation,
                 cds,
                 IPA_region,
                 IPUI,
                 IPA_ID):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.exon = exon
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.five_UTR = five_UTR
        self.strand = strand
        self.annotation = annotation
        self.cds = cds
        self.IPA_region = IPA_region
        self.IPUI = IPUI
        self.IPA_ID = IPA_ID
        self.three_UTR = None
                        
    def to_dict(self):
        return {"gene_id":self.gene_id,
                "gene_name":self.gene_name,
                "transcript_id":self.transcript_id,
                "exon":tuple(self.exon),
                "start_codon":self.start_codon,
                "stop_codon":self.stop_codon,
                "five_UTR":tuple(self.five_UTR),
                "strand":self.strand,
                "annotation":self.annotation,
                "CDS":self.cds,
                "IPA_region":self.IPA_region,
                "IPUI":self.IPUI,
                "IPA_ID":self.IPA_ID,
                "three_UTR":self.three_UTR}

    
def compositeTranscript(rt_te_composite):
    composite_transcript_dic = dict()
    for x in rt_te_composite.itertuples():
        region = x.region1
        chrom,te_start,te_end,strand = region.split(":")
        te_start,te_end = int(te_start),int(te_end)
        gene_id = x.gene_id
        transcript_id = x.transcript_id
        possible_up_site = x.up_site
        IPA_region = x.region2
        IPUI = x.IPUI
        IPA_ID = x.IPA_ID
        gene_info = gtf_iv_df.loc[gtf_iv_df["gene_id"] == gene_id,]
        if possible_up_site == "first_exon":
            novel_transcript_exon = [(te_start,te_end)]
            transcript_info_df = gene_info.loc[gene_info["transcript_id"]==transcript_id,]
            gene_name = list(set(transcript_info_df["gene_name"]))
            assert len(gene_name) == 1
            gene_name = gene_name[0]
            # start codon
            start_codon = transcript_info_df.loc[transcript_info_df["type"] == "start_codon"]
            start_codon_frame = dict(zip(list(zip(start_codon["start"],start_codon["end"])),start_codon["frame"]))
            start_codon = tuple(zip(start_codon["start"],start_codon["end"]))
            start_codon_one_dimen = list(set(np.array(start_codon).flatten()))
            # stop_codon 
            stop_codon = transcript_info_df.loc[transcript_info_df["type"] == "stop_codon"]
            stop_codon_frame = dict(zip(list(zip(stop_codon["start"],stop_codon["end"])),stop_codon["frame"]))
            stop_codon = tuple(zip(stop_codon["start"],stop_codon["end"]))
            stop_codon_one_dimen = list(set(np.array(stop_codon).flatten()))
            # 5'UTR
            five_UTR = transcript_info_df.loc[transcript_info_df["type"] == "5UTR"]
            five_UTR = list(zip(five_UTR["start"],five_UTR["end"]))
            # cds region
            cds_region_df = transcript_info_df.loc[transcript_info_df["type"] == "CDS"]
            cds_region = dict(zip(list(zip(cds_region_df["start"],cds_region_df["end"])),cds_region_df["frame"]))
            annotation = None
            if strand == "+":    
                if min(start_codon_one_dimen) > te_end:
                    annotation = "non_coding"  
            elif strand == "-":
                if max(start_codon_one_dimen) < te_start:
                    annotation = "non_coding"
            novel_transcript = novelTranscript(gene_id = gene_id,
                                               gene_name = gene_name,
                                               transcript_id = transcript_id,
                                               exon = novel_transcript_exon,
                                               start_codon = start_codon_frame,
                                               stop_codon = stop_codon_frame,
                                               five_UTR = five_UTR,
                                               strand = strand,
                                               annotation = annotation,
                                               cds = cds_region,
                                               IPA_region = IPA_region,
                                               IPUI = IPUI,
                                               IPA_ID = IPA_ID)
            composite_transcript_dic.setdefault(region,[]).append(novel_transcript)
        else:
            up_site = int(possible_up_site.split(":")[1])
            exon_info = gtf_iv_df.loc[gtf_iv_df["gene_id"] == gene_id,]
            gene_name = list(set(exon_info["gene_name"]))
            assert len(gene_name) == 1
            gene_name = gene_name[0]
            if strand == "+":
                possible_transcripts = exon_info.loc[exon_info["end"] == up_site,]
                possible_transcripts = possible_transcripts.loc[possible_transcripts["type"] == "exon",]
                transcript_up_exon_dic = {x.transcript_id:(x.start,x.end) for x in possible_transcripts.itertuples()}
                filtered_transcripts = list(transcript_up_exon_dic.keys())
                
                if len(filtered_transcripts) == 1:
                    filtered_transcript = filtered_transcripts[0]
                elif len(filtered_transcripts) > 1:
                    filtered_transcript_df = exon_info.loc[gtf_iv_df["transcript_id"].isin(filtered_transcripts),]
                    filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["type"]=="exon",]
                    filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["start"]<=te_start,]
                    filtered_transcript_df["len"] = filtered_transcript_df["end"] - filtered_transcript_df["start"] + 1
                    filtered_transcript_len_df = filtered_transcript_df.loc[:,["transcript_id","len"]].groupby(by=["transcript_id"]).sum()
                    filtered_transcript = filtered_transcript_len_df["len"].idxmax()
                                                               
                # up exon
                up_exon = transcript_up_exon_dic[filtered_transcript]
                # exon of transcript
                transcript_exon_info = exon_info.loc[(exon_info["type"] == "exon") & (exon_info["transcript_id"] == filtered_transcript)]
                transcript_exon_info = transcript_exon_info.sort_values(by = ["chrom","start"])
                exon_coordinate = list(zip(transcript_exon_info["start"],transcript_exon_info["end"]))
                exon_coordinate = sorted(exon_coordinate,key = lambda x:x[0])
                # location
                up_exon_index = exon_coordinate.index(up_exon)
                if len(exon_coordinate) > (up_exon_index + 2):
                    template_exon = exon_coordinate[up_exon_index + 1]
                    if template_exon[0] == te_start:
                        if exon_coordinate[-1][1] > te_end:
                            novel_transcript_info = transcript_exon_info.loc[transcript_exon_info["end"] <= up_site,]
                            novel_transcript_exon = list(zip(novel_transcript_info["start"],
                                                             novel_transcript_info["end"]))
                            novel_transcript_exon.append((te_start,te_end))
                            # info of transcript
                            transcript_info = exon_info.loc[(exon_info["transcript_id"] == filtered_transcript)]
                            # start codon
                            start_codon = transcript_info.loc[transcript_info["type"] == "start_codon"]
                            start_codon_frame = dict(zip(list(zip(start_codon["start"],start_codon["end"])),start_codon["frame"]))
                            start_codon = tuple(zip(start_codon["start"],start_codon["end"]))
                            start_codon_one_dimen = list(set(np.array(start_codon).flatten()))
                            # stop_codon 
                            stop_codon = transcript_info.loc[transcript_info["type"] == "stop_codon"]
                            stop_codon_frame = dict(zip(list(zip(stop_codon["start"],stop_codon["end"])),stop_codon["frame"]))
                            stop_codon = tuple(zip(stop_codon["start"],stop_codon["end"]))
                            stop_codon_one_dimen = list(set(np.array(stop_codon).flatten()))
                            # 5'UTR
                            five_UTR = transcript_info.loc[transcript_info["type"] == "5UTR"]
                            five_UTR = list(zip(five_UTR["start"],five_UTR["end"]))
                            # cds region
                            cds_region_df = transcript_info.loc[transcript_info["type"] == "CDS"]
                            cds_region = dict(zip(list(zip(cds_region_df["start"],cds_region_df["end"])),cds_region_df["frame"]))
                            # annotation
                            annotation = None
                            if min(start_codon_one_dimen) > te_end:
                                annotation = "non_coding"
                            elif max(stop_codon_one_dimen) < te_start:
                                annotation = "protein_coding"                                
                            novel_transcript = novelTranscript(gene_id = gene_id,
                                                               gene_name = gene_name,
                                                               transcript_id = filtered_transcript,
                                                               exon = novel_transcript_exon,
                                                               start_codon = start_codon_frame,
                                                               stop_codon = stop_codon_frame,
                                                               five_UTR = five_UTR,
                                                               strand = strand,
                                                               annotation = annotation,
                                                               cds = cds_region,
                                                               IPA_region = IPA_region,
                                                               IPUI = IPUI,
                                                               IPA_ID = IPA_ID)
                            composite_transcript_dic.setdefault(region,[]).append(novel_transcript)                    
            elif strand == "-":
                possible_transcripts = exon_info.loc[exon_info["start"] == up_site,]
                possible_transcripts = possible_transcripts.loc[possible_transcripts["type"] == "exon",]
                transcript_up_exon_dic = {x.transcript_id:(x.start,x.end) for x in possible_transcripts.itertuples()}
                filtered_transcripts = list(transcript_up_exon_dic.keys())
                
                if len(filtered_transcripts) == 1:
                    filtered_transcript = filtered_transcripts[0]
                elif len(filtered_transcripts) > 1:
                    filtered_transcript_df = exon_info.loc[gtf_iv_df["transcript_id"].isin(filtered_transcripts),]
                    filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["type"]=="exon",]
                    filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["end"]>=te_end,]
                    filtered_transcript_df["len"] = filtered_transcript_df["end"] - filtered_transcript_df["start"] + 1
                    filtered_transcript_len_df = filtered_transcript_df.loc[:,["transcript_id","len"]].groupby(by=["transcript_id"]).sum()
                    filtered_transcript = filtered_transcript_len_df["len"].idxmax() 
                # up exon
                up_exon = transcript_up_exon_dic[filtered_transcript]
                # exon of transcript
                transcript_exon_info = exon_info.loc[(exon_info["type"] == "exon") & (exon_info["transcript_id"] == filtered_transcript)]
                transcript_exon_info = transcript_exon_info.sort_values(by = ["chrom","start"],
                                                                        ascending = False)
                exon_coordinate = list(zip(transcript_exon_info["start"],transcript_exon_info["end"]))
                exon_coordinate = sorted(exon_coordinate,key = lambda x:x[0],reverse = True)
                # location
                up_exon_index = exon_coordinate.index(up_exon)
                if len(exon_coordinate) > (up_exon_index + 2):
                    template_exon = exon_coordinate[up_exon_index+1]
                    if template_exon[1] == te_end:
                        if exon_coordinate[-1][0] < te_start:
                            novel_transcript_info = transcript_exon_info.loc[transcript_exon_info["start"] >= up_site,]
                            novel_transcript_exon = list(zip(novel_transcript_info["start"],
                                                             novel_transcript_info["end"]))
                            novel_transcript_exon.append((te_start,te_end))
                            # info of transcript
                            transcript_info = exon_info.loc[(exon_info["transcript_id"] == filtered_transcript)]
                            # start codon
                            start_codon = transcript_info.loc[transcript_info["type"] == "start_codon"]
                            start_codon_frame = dict(zip(list(zip(start_codon["start"],start_codon["end"])),start_codon["frame"]))
                            start_codon = tuple(zip(start_codon["start"],start_codon["end"]))
                            start_codon_one_dimen = list(set(np.array(start_codon).flatten()))
                            # stop_codon 
                            stop_codon = transcript_info.loc[transcript_info["type"] == "stop_codon"]
                            stop_codon_frame = dict(zip(list(zip(stop_codon["start"],stop_codon["end"])),stop_codon["frame"]))
                            stop_codon = tuple(zip(stop_codon["start"],stop_codon["end"]))
                            stop_codon_one_dimen = list(set(np.array(stop_codon).flatten()))
                            # 5'UTR
                            five_UTR = transcript_info.loc[transcript_info["type"] == "5UTR"]
                            five_UTR = list(zip(five_UTR["start"],five_UTR["end"]))
                            # cds region
                            cds_region_df = transcript_info.loc[transcript_info["type"] == "CDS"]
                            cds_region = dict(zip(list(zip(cds_region_df["start"],cds_region_df["end"])),cds_region_df["frame"]))
                            # annotation
                            annotation = None
                            if max(start_codon_one_dimen) < te_start:
                                annotation = "non_coding"
                            elif min(stop_codon_one_dimen) > te_end:
                                annotation = "protein_coding"                               
                            novel_transcript = novelTranscript(gene_id = gene_id,
                                                               gene_name = gene_name,
                                                               transcript_id = filtered_transcript,
                                                               exon = novel_transcript_exon,
                                                               start_codon = start_codon_frame,
                                                               stop_codon= stop_codon_frame,
                                                               five_UTR = five_UTR,
                                                               strand = strand,
                                                               annotation = annotation,
                                                               cds = cds_region,
                                                               IPA_region = IPA_region,
                                                               IPUI = IPUI,
                                                               IPA_ID = IPA_ID)
                            composite_transcript_dic.setdefault(region,[]).append(novel_transcript)
    return composite_transcript_dic   


def checkComTrans(composite_transcript_dic,
                  genome):
    start_codon = "ATG"
    for te_region,transcripts in composite_transcript_dic.items():
        chrom,te_start,te_end,strand = te_region.split(":")
        te_start,te_end = int(te_start),int(te_end)
        for index,transcript in enumerate(transcripts):
            annotation = transcript.annotation
            if annotation == None:
                exon_info = transcript.exon
                exon_coordinate_dic = {y:i for i,x in enumerate(exon_info) for y in x}
                exon_coordinate = list(exon_coordinate_dic.keys())
                st_condon_info = transcript.start_codon
                st_condon_info = list(st_condon_info.keys())
                st_condon_info = [y for x in st_condon_info for y in x]
                st_condon_info = list(set(st_condon_info).difference(set(exon_coordinate)))
                exon_st_list = exon_coordinate + st_condon_info
                exon_st_list = sorted(exon_st_list)
                min_st_index = exon_st_list.index(min(st_condon_info))
                max_st_index = exon_st_list.index(max(st_condon_info))
                if (min_st_index > 0) and (max_st_index + 1 < len(exon_st_list)):
                    left_exon_loc = exon_coordinate_dic[exon_st_list[min_st_index-1]]
                    right_exon_loc = exon_coordinate_dic[exon_st_list[max_st_index+1]]
                else:
                    left_exon_loc = -1
                    right_exon_loc = -2 
                if left_exon_loc != right_exon_loc:
                    te_seqs = genome.sequence({'chr': chrom,
                                               'start': te_start,
                                               'stop': te_end,
                                               'strand': strand})
                    te_seqs = te_seqs.upper()
                    start_index = te_seqs.find(start_codon)
                    if start_index >= 0:
                        if strand == "+":
                            novel_st_codon_st = te_start + start_index
                            novel_st_codon = {(novel_st_codon_st,novel_st_codon_st+2):0}
                        elif strand == "-":
                            novel_st_codon_end = te_end - start_index
                            novel_st_codon = {(novel_st_codon_end-2,novel_st_codon_end):0}
                        composite_transcript_dic[te_region][index].start_codon = novel_st_codon
                    else:
                        composite_transcript_dic[te_region][index].annotation = "non_coding"
    return


def determinCoding(region_transcript_dic,
                   genome):
    stop_codons = ['TAA', 'TAG', 'TGA']
    region_transcript_annotation_dic = dict()
    protein_regions = set()
    for region,transcripts in region_transcript_dic.items():
        transcripts = [transcript.to_dict() for transcript in transcripts]
        # delete same interval
        if len(transcripts) > 1:
            transcript_df = pd.DataFrame(transcripts).drop_duplicates(subset = ["gene_id",
                                                                                "gene_name",
                                                                                "strand",
                                                                                "exon",
                                                                                "five_UTR",
                                                                                "annotation"],keep="first")
            transcripts = list(transcript_df.to_dict(orient="index").values())
        for transcript in transcripts:
            annotation = transcript["annotation"]
            if annotation == "non_coding":
                transcript["CDS"] = None
                transcript["start_codon"] = None
                transcript["stop_codon"] = None
                region_transcript_annotation_dic.setdefault(region,[]).append(transcript)
            elif annotation == "protein_coding":
                region_transcript_annotation_dic.setdefault(region,[]).append(transcript)
                protein_regions.add(region)
            else:
                chrom,gene_id = transcript["gene_id"].split(":")
                strand = transcript["strand"]
                all_exons = list(transcript["exon"])
                start_codon_frame = transcript["start_codon"]
                start_codon = tuple(start_codon_frame.keys())
                start_codon_first = None
                if strand == "+":
                    start_codon = sorted(start_codon)
                    if len(start_codon) == 2:
                        start_codon_first = start_codon[0]
                        start_codon = start_codon[1]
                    elif len(start_codon) == 1:
                        start_codon = start_codon[0]
                    check_start_codon_border = start_codon[1] in np.array(all_exons).flatten()
                elif strand == "-":
                    start_codon = sorted(start_codon,reverse=True)
                    if len(start_codon) == 2:
                        start_codon_first = start_codon[0]
                        start_codon = start_codon[1]
                    elif len(start_codon) == 1:
                        start_codon = start_codon[0]
                    check_start_codon_border = start_codon[0] in np.array(all_exons).flatten()
                all_exons.append(start_codon)
                all_exon_coordinates = np.array(list(set([y for x in all_exons for y in x])))
                cds_regions = transcript["CDS"]
                if strand == "+":
                    all_exon_coordinates = list(all_exon_coordinates[all_exon_coordinates >= start_codon[0]])
                    if not check_start_codon_border:
                        all_exon_coordinates.remove(start_codon[1])
                    all_exon_coordinates = sorted(all_exon_coordinates)
                    potential_cds_regions = [(all_exon_coordinates[x],
                                              all_exon_coordinates[x+1]) for x in range(0,len(all_exon_coordinates),2)]
                    if start_codon_first:
                        potential_cds_regions.append(start_codon_first)
                    potential_cds_regions = sorted(potential_cds_regions)
                    if len(potential_cds_regions) == 1:
                        potential_cds_regions = potential_cds_regions[0]
                        potential_cds_region_seq = genome.sequence(
                                            {
                                                'chr': chrom,
                                                'start': potential_cds_regions[0],
                                                'stop': potential_cds_regions[1],
                                                'strand': strand
                                            }
                                        )
                        stop_codon_position_index = None
                        for i in range(0,len(potential_cds_region_seq),3):
                            if potential_cds_region_seq[i:i+3].upper() in stop_codons:
                                stop_codon_position_index = i
                                break
                        if stop_codon_position_index is not None:
                            stop_codon_position = (potential_cds_regions[0]+stop_codon_position_index,
                                                   potential_cds_regions[0]+stop_codon_position_index+2)
                            stop_codon_position_frame = {stop_codon_position:0}
                            novel_cds = (potential_cds_regions[0],
                                         potential_cds_regions[0]+stop_codon_position_index-1)
                            transcript["annotation"] = "protein_coding"
                            transcript["CDS"] = {novel_cds:0}
                            transcript["stop_codon"] = stop_codon_position_frame
                            protein_regions.add(region)
                        else:
                            transcript["annotation"] = "non_coding"
                            transcript["CDS"] = None
                        region_transcript_annotation_dic.setdefault(region,[]).append(transcript)
                    else:
                        exist_cds_region = potential_cds_regions[:-1]
                        assert len(set(exist_cds_region).difference(set(cds_regions.keys()))) == 0
                        last_exon = potential_cds_regions[-1]
                        penultimate_exon = potential_cds_regions[-2]
                        penultimate_exon_frame = cds_regions[penultimate_exon]
                        penultimate_CDS_seq = genome.sequence(
                            {
                                'chr': chrom,
                                'start': penultimate_exon[0] + penultimate_exon_frame,
                                'stop': penultimate_exon[1],
                                'strand': strand
                            }
                        )
                        bases_left = len(penultimate_CDS_seq) % 3
                        if bases_left == 0:
                            possible_last_exon_frame = 0
                            previous_bases = ''
                        elif bases_left > 0:
                            possible_last_exon_frame = 3 - bases_left
                            previous_bases = penultimate_CDS_seq[-bases_left:]
                        last_exon_CDS_seq = genome.sequence(
                            {
                                'chr': chrom,
                                'start': last_exon[0],
                                'stop': last_exon[1],
                                'strand': strand
                            }
                        )
                        up_tail_merge_last_exon_CDS_seq = previous_bases + last_exon_CDS_seq
                        stop_codon_position_index = None
                        for i in range(0,len(up_tail_merge_last_exon_CDS_seq),3):
                            if up_tail_merge_last_exon_CDS_seq[i:i+3].upper() in stop_codons:
                                stop_codon_position_index = i
                                break
                        if stop_codon_position_index is not None:
                            protein_regions.add(region)
                            if (stop_codon_position_index == 0) and (bases_left > 0):
                                if bases_left == 1:
                                    stop_codon_position_first = (penultimate_exon[-1],penultimate_exon[-1])
                                    stop_codon_position_second = (last_exon[0],last_exon[0]+1)
                                    stop_codon_position_frame = {stop_codon_position_first:0,
                                                                 stop_codon_position_second:2}
                                    novel_cds_frame = {x:cds_regions[x] for x in exist_cds_region[:-1]}
                                    novel_cds_frame.update({(penultimate_exon[0],penultimate_exon[-1]-1):penultimate_exon_frame})
                                    transcript["annotation"] = "protein_coding"
                                    transcript["CDS"] = novel_cds_frame
                                    transcript["stop_codon"] = stop_codon_position_frame
                                elif bases_left == 2:
                                    stop_codon_position_first = (penultimate_exon[-1]-1,penultimate_exon[-1])
                                    stop_codon_position_second = (last_exon[0],last_exon[0])
                                    stop_codon_position_frame = {stop_codon_position_first:0,
                                                                 stop_codon_position_second:1}
                                    novel_cds_frame = {x:cds_regions[x] for x in exist_cds_region[:-1]}
                                    novel_cds_frame.update({(penultimate_exon[0],penultimate_exon[-1]-2):penultimate_exon_frame})
                                    transcript["annotation"] = "protein_coding"
                                    transcript["CDS"] = novel_cds_frame
                                    transcript["stop_codon"] = stop_codon_position_frame                                    
                            elif (stop_codon_position_index >= bases_left):                
                                stop_codon_position = (last_exon[0]-bases_left+stop_codon_position_index,
                                                       last_exon[0]-bases_left+stop_codon_position_index+2)
                                stop_codon_position_frame = {stop_codon_position:0}
                                novel_cds_frame = {x:cds_regions[x] for x in exist_cds_region}                        
                                last_cds = (last_exon[0],last_exon[0]-bases_left+stop_codon_position_index-1)
                                if stop_codon_position_index > 0:
                                    novel_cds_frame.update({last_cds:possible_last_exon_frame})
                                transcript["annotation"] = "protein_coding"
                                transcript["CDS"] = novel_cds_frame
                                transcript["stop_codon"] = stop_codon_position_frame
                        else:
                            transcript["annotation"] = "non_coding"
                            transcript["CDS"] = None
                        region_transcript_annotation_dic.setdefault(region,[]).append(transcript)
                elif strand == "-":                
                    all_exon_coordinates = list(all_exon_coordinates[all_exon_coordinates <= start_codon[1]])
                    if not check_start_codon_border:
                        all_exon_coordinates.remove(start_codon[0])
                    all_exon_coordinates = sorted(all_exon_coordinates,reverse=True)
                    potential_cds_regions = [(all_exon_coordinates[x+1],all_exon_coordinates[x]) for x in range(0,len(all_exon_coordinates),2)]
                    if start_codon_first:
                        potential_cds_regions.append(start_codon_first)
                    potential_cds_regions = sorted(potential_cds_regions,reverse=True)
                    if len(potential_cds_regions) == 1:
                        potential_cds_regions = potential_cds_regions[0]
                        potential_cds_region_seq = genome.sequence(
                                            {
                                                'chr': chrom,
                                                'start': potential_cds_regions[0],
                                                'stop': potential_cds_regions[1],
                                                'strand': strand
                                            }
                                        )
                        stop_codon_position_index = None
                        for i in range(0,len(potential_cds_region_seq),3):
                            if potential_cds_region_seq[i:i+3].upper() in stop_codons:
                                stop_codon_position_index = i
                                break
                        if stop_codon_position_index is not None:
                            stop_codon_position = (potential_cds_regions[1]-stop_codon_position_index-2,
                                                   potential_cds_regions[1]-stop_codon_position_index)
                            stop_codon_position_frame = {stop_codon_position:0}
                            novel_cds = (potential_cds_regions[1]-stop_codon_position_index+1,
                                         potential_cds_regions[1])
                            transcript["annotation"] = "protein_coding"
                            transcript["CDS"] = {novel_cds:0}
                            transcript["stop_codon"] = stop_codon_position_frame
                            protein_regions.add(region)
                        else:
                            transcript["annotation"] = "non_coding"
                            transcript["CDS"] = None
                        region_transcript_annotation_dic.setdefault(region,[]).append(transcript)
                    else:                    
                        exist_cds_region = potential_cds_regions[:-1]
                        assert len(set(exist_cds_region).difference(set(cds_regions.keys()))) == 0
                        last_exon = potential_cds_regions[-1]
                        penultimate_exon = potential_cds_regions[-2]
                        penultimate_exon_frame = cds_regions[penultimate_exon]
                        penultimate_CDS_seq = genome.sequence(
                            {
                                'chr': chrom,
                                'start': penultimate_exon[0],
                                'stop': penultimate_exon[1] - penultimate_exon_frame,
                                'strand': strand
                            }
                        )
                        bases_left = len(penultimate_CDS_seq) % 3
                        if bases_left == 0:
                            possible_last_exon_frame = 0
                            previous_bases = ''
                        elif bases_left > 0:
                            possible_last_exon_frame = 3 - bases_left
                            previous_bases = penultimate_CDS_seq[-bases_left:]
                        last_exon_CDS_seq = genome.sequence(
                            {
                                'chr': chrom,
                                'start': last_exon[0],
                                'stop': last_exon[1],
                                'strand': strand
                            }
                        )
                        up_tail_merge_last_exon_CDS_seq = previous_bases + last_exon_CDS_seq
                        stop_codon_position_index = None
                        for i in range(0,len(up_tail_merge_last_exon_CDS_seq),3):
                            if up_tail_merge_last_exon_CDS_seq[i:i+3].upper() in stop_codons:
                                stop_codon_position_index = i
                                break
                        if stop_codon_position_index is not None:
                            protein_regions.add(region)
                            if (stop_codon_position_index == 0) and (bases_left > 0):
                                if bases_left == 1:
                                    stop_codon_position_first = (penultimate_exon[0],penultimate_exon[0])
                                    stop_codon_position_second = (last_exon[-1]-1,last_exon[-1])
                                    stop_codon_position_frame = {stop_codon_position_first:0,
                                                                 stop_codon_position_second:2}
                                    novel_cds_frame = {x:cds_regions[x] for x in exist_cds_region[:-1]}
                                    novel_cds_frame.update({(penultimate_exon[0]+1,penultimate_exon[1]):penultimate_exon_frame})
                                    transcript["annotation"] = "protein_coding"
                                    transcript["CDS"] = novel_cds_frame
                                    transcript["stop_codon"] = stop_codon_position_frame
                                elif bases_left == 2:
                                    stop_codon_position_first = (penultimate_exon[0],penultimate_exon[0]+1)
                                    stop_codon_position_second = (last_exon[-1],last_exon[-1])
                                    stop_codon_position_frame = {stop_codon_position_first:0,
                                                                 stop_codon_position_second:1}
                                    novel_cds_frame = {x:cds_regions[x] for x in exist_cds_region[:-1]}
                                    novel_cds_frame.update({(penultimate_exon[0]+2,penultimate_exon[1]):penultimate_exon_frame})
                                    transcript["annotation"] = "protein_coding"
                                    transcript["CDS"] = novel_cds_frame
                                    transcript["stop_codon"] = stop_codon_position_frame                                                                                  
                            elif (stop_codon_position_index >= bases_left):                
                                stop_codon_position = (last_exon[1]+bases_left-stop_codon_position_index-2,
                                                       last_exon[1]+bases_left-stop_codon_position_index)
                                stop_codon_position_frame = {stop_codon_position:0}
                                novel_cds_frame = {x:cds_regions[x] for x in exist_cds_region}                        
                                last_cds = (last_exon[1]+bases_left-stop_codon_position_index+1,last_exon[1])    
                                if stop_codon_position_index > 0:
                                    novel_cds_frame.update({last_cds:possible_last_exon_frame})
                                transcript["annotation"] = "protein_coding"
                                transcript["CDS"] = novel_cds_frame
                                transcript["stop_codon"] = stop_codon_position_frame
                        else:
                            transcript["annotation"] = "non_coding"
                            transcript["CDS"] = None
                        region_transcript_annotation_dic.setdefault(region,[]).append(transcript)
    return region_transcript_annotation_dic


def skipTranscript(rt_te_skip):
    skip_transcript_dic = dict()
    for x in rt_te_skip.itertuples():
        region = x.region1
        chrom,te_start,te_end,strand = region.split(":")
        te_start,te_end = int(te_start),int(te_end)
        gene_id = x.gene_id
        possible_up_site = x.up_site
        up_site = int(possible_up_site.split(":")[1])
        exon_info = gtf_iv_df.loc[gtf_iv_df["gene_id"] == gene_id,]
        gene_name = list(set(exon_info["gene_name"]))
        IPA_region = x.region2
        IPUI = x.IPUI
        IPA_ID = x.IPA_ID
        assert len(gene_name) == 1
        gene_name = gene_name[0]
        if strand == "+":
            possible_transcripts = exon_info.loc[exon_info["end"] == up_site,]
            possible_transcripts = possible_transcripts.loc[possible_transcripts["type"] == "exon",]
            transcript_up_exon_dic = {x.transcript_id:(x.start,x.end) for x in possible_transcripts.itertuples()}
            filtered_transcripts = list(transcript_up_exon_dic.keys())
            
            if len(filtered_transcripts) == 1:
                filtered_transcript = filtered_transcripts[0]
            elif len(filtered_transcripts) > 1:
                filtered_transcript_df = exon_info.loc[gtf_iv_df["transcript_id"].isin(filtered_transcripts),]
                filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["type"]=="exon",]
                filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["start"]<=te_start,]
                filtered_transcript_df["len"] = filtered_transcript_df["end"] - filtered_transcript_df["start"] + 1
                filtered_transcript_len_df = filtered_transcript_df.loc[:,["transcript_id","len"]].groupby(by=["transcript_id"]).sum()
                filtered_transcript = filtered_transcript_len_df["len"].idxmax()
            # up exon
            up_exon = transcript_up_exon_dic[filtered_transcript]
            # exon of transcript
            transcript_exon_info = exon_info.loc[(exon_info["type"] == "exon") & (exon_info["transcript_id"] == filtered_transcript)]
            transcript_exon_info = transcript_exon_info.sort_values(by = ["chrom","start"])
            exon_coordinate = list(zip(transcript_exon_info["start"],transcript_exon_info["end"]))
            exon_coordinate = sorted(exon_coordinate,key = lambda x:x[0])
            # location
            up_exon_index = exon_coordinate.index(up_exon)
            novel_transcript_info = transcript_exon_info.loc[transcript_exon_info["end"] <= up_site,]
            novel_transcript_exon = list(zip(novel_transcript_info["start"],
                                             novel_transcript_info["end"]))
            novel_transcript_exon.append((te_start,te_end))
            # info of transcript
            transcript_info = exon_info.loc[(exon_info["transcript_id"] == filtered_transcript)]
            # start codon
            start_codon = transcript_info.loc[transcript_info["type"] == "start_codon"]
            start_codon_frame = dict(zip(list(zip(start_codon["start"],start_codon["end"])),start_codon["frame"]))
            start_codon = tuple(zip(start_codon["start"],start_codon["end"]))
            start_codon_one_dimen = list(set(np.array(start_codon).flatten()))
            # stop_codon 
            stop_codon = transcript_info.loc[transcript_info["type"] == "stop_codon"]
            stop_codon_frame = dict(zip(list(zip(stop_codon["start"],stop_codon["end"])),stop_codon["frame"]))
            stop_codon = tuple(zip(stop_codon["start"],stop_codon["end"]))
            stop_codon_one_dimen = list(set(np.array(stop_codon).flatten()))
            # 5'UTR
            five_UTR = transcript_info.loc[transcript_info["type"] == "5UTR"]
            five_UTR = list(zip(five_UTR["start"],five_UTR["end"]))
            # cds region
            cds_region_df = transcript_info.loc[transcript_info["type"] == "CDS"]
            cds_region = dict(zip(list(zip(cds_region_df["start"],cds_region_df["end"])),cds_region_df["frame"]))
            # annotation
            annotation = None
            if min(start_codon_one_dimen) > te_end:
                annotation = "non_coding"
            elif max(stop_codon_one_dimen) < te_end:
                annotation = "protein_coding"                            
            novel_transcript = novelTranscript(gene_id = gene_id,
                                               gene_name = gene_name,
                                               transcript_id = filtered_transcript,
                                               exon = novel_transcript_exon,
                                               start_codon = start_codon_frame,
                                               stop_codon = stop_codon_frame,
                                               five_UTR = five_UTR,
                                               strand = strand,
                                               annotation = annotation,
                                               cds = cds_region,
                                               IPA_region = IPA_region,
                                               IPUI = IPUI,
                                               IPA_ID = IPA_ID)
            skip_transcript_dic.setdefault(region,[]).append(novel_transcript)
        elif strand == "-":
            possible_transcripts = exon_info.loc[exon_info["start"] == up_site,]
            possible_transcripts = possible_transcripts.loc[possible_transcripts["type"] == "exon",]
            transcript_up_exon_dic = {x.transcript_id:(x.start,x.end) for x in possible_transcripts.itertuples()}
            filtered_transcripts = list(transcript_up_exon_dic.keys())
            
            if len(filtered_transcripts) == 1:
                filtered_transcript = filtered_transcripts[0]
            elif len(filtered_transcripts) > 1:
                filtered_transcript_df = exon_info.loc[gtf_iv_df["transcript_id"].isin(filtered_transcripts),]
                filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["type"]=="exon",]
                filtered_transcript_df = filtered_transcript_df.loc[filtered_transcript_df["end"]>=te_end,]
                filtered_transcript_df["len"] = filtered_transcript_df["end"] - filtered_transcript_df["start"] + 1
                filtered_transcript_len_df = filtered_transcript_df.loc[:,["transcript_id","len"]].groupby(by=["transcript_id"]).sum()
                filtered_transcript = filtered_transcript_len_df["len"].idxmax() 
            
            # up exon
            up_exon = transcript_up_exon_dic[filtered_transcript]
            # exon of transcript
            transcript_exon_info = exon_info.loc[(exon_info["type"] == "exon") & (exon_info["transcript_id"] == filtered_transcript)]
            transcript_exon_info = transcript_exon_info.sort_values(by = ["chrom","start"],
                                                                    ascending = False)
            exon_coordinate = list(zip(transcript_exon_info["start"],transcript_exon_info["end"]))
            exon_coordinate = sorted(exon_coordinate,key = lambda x:x[0],reverse = True)
            # location
            up_exon_index = exon_coordinate.index(up_exon)
            novel_transcript_info = transcript_exon_info.loc[transcript_exon_info["start"] >= up_site,]
            novel_transcript_exon = list(zip(novel_transcript_info["start"],
                                             novel_transcript_info["end"]))
            novel_transcript_exon.append((te_start,te_end))
            # info of transcript
            transcript_info = exon_info.loc[(exon_info["transcript_id"] == filtered_transcript)]
            # start codon
            start_codon = transcript_info.loc[transcript_info["type"] == "start_codon"]
            start_codon_frame = dict(zip(list(zip(start_codon["start"],start_codon["end"])),start_codon["frame"]))
            start_codon = tuple(zip(start_codon["start"],start_codon["end"]))
            start_codon_one_dimen = list(set(np.array(start_codon).flatten()))
            # stop_codon 
            stop_codon = transcript_info.loc[transcript_info["type"] == "stop_codon"]
            stop_codon_frame = dict(zip(list(zip(stop_codon["start"],stop_codon["end"])),stop_codon["frame"]))
            stop_codon = tuple(zip(stop_codon["start"],stop_codon["end"]))
            stop_codon_one_dimen = list(set(np.array(stop_codon).flatten()))
            # 5'UTR
            five_UTR = transcript_info.loc[transcript_info["type"] == "5UTR"]
            five_UTR = list(zip(five_UTR["start"],five_UTR["end"]))
            # cds region
            cds_region_df = transcript_info.loc[transcript_info["type"] == "CDS"]
            cds_region = dict(zip(list(zip(cds_region_df["start"],cds_region_df["end"])),cds_region_df["frame"]))
            # annotation
            annotation = None
            if max(start_codon_one_dimen) < te_start:
                annotation = "non_coding"
            elif max(stop_codon_one_dimen) >= up_site:
                annotation = "protein_coding" 
            novel_transcript = novelTranscript(gene_id = gene_id,
                                               gene_name = gene_name,
                                               transcript_id = filtered_transcript,
                                               exon = novel_transcript_exon,
                                               start_codon = start_codon_frame,
                                               stop_codon = stop_codon_frame,                                                        
                                               five_UTR = five_UTR,
                                               strand = strand,
                                               annotation = annotation,
                                               cds = cds_region,
                                               IPA_region = IPA_region,
                                               IPUI = IPUI,
                                               IPA_ID = IPA_ID)
            skip_transcript_dic.setdefault(region,[]).append(novel_transcript)    
    return skip_transcript_dic


def transGTFormat(region,
                  transcript,
                  te_type,
                  transcript_source = "IPAfinder_annotated"):
    gtf_info = []
    IPA_info = []
    chrom,gene_id = transcript["gene_id"].split(":")
    gene_id = 'gene_id "' + gene_id + '";'
    gene_name = transcript["gene_name"]
    gene_name = 'gene_name "' + gene_name + '";'
    all_exons = transcript["exon"]
    all_exon_coordinates = np.array(all_exons).flatten()
    cds_frame = transcript["CDS"]      
    strand = transcript["strand"]
    annotation = transcript["annotation"]
    template_transcript = transcript["transcript_id"]
    IPA_region = transcript["IPA_region"]
    IPUI = transcript["IPUI"]
    IPA_ID = transcript["IPA_ID"]
    if annotation == "protein_coding":
        start_codon_frame = transcript["start_codon"]
        start_codons = tuple(start_codon_frame.keys())
        stop_codon_frame = transcript["stop_codon"]
        stop_codons = tuple(stop_codon_frame.keys())
        cds_regions = list(cds_frame.keys())
        novel_transcript_id = "".join(["IPA_XM_",te_type[0],IPA_ID]) 
        novel_transcript_id2 = "|".join(["novel_XM",template_transcript,te_type,"_".join(region.split(":")[1:3])]) 
    elif annotation == "non_coding":
        novel_transcript_id = "".join(["IPA_XR_",te_type[0],IPA_ID]) 
        novel_transcript_id2 = "|".join(["novel_XR",template_transcript,te_type,"_".join(region.split(":")[1:3])])
    IPA_info.append([transcript["gene_name"],IPA_region,te_type,IPUI,novel_transcript_id,novel_transcript_id2])
    exon_id = 'exon_id "' + novel_transcript_id + '.%d";'
    novel_transcript_id = 'transcript_id "' + novel_transcript_id + '";'
    novel_transcript_id2 = 'transcript_id_long "' + novel_transcript_id2 + '";'
    # transcript info
    transcript_attr = "".join([gene_id,novel_transcript_id,novel_transcript_id2,gene_name])
    transcript_info = [chrom,transcript_source,"transcript",min(all_exon_coordinates),
                       max(all_exon_coordinates),".",strand,".",transcript_attr]
    gtf_info.append(transcript_info)        
    # exon and cds
    exon_number = 'exon_number "%d";'        
    exon_attr = "".join([gene_id,novel_transcript_id,exon_number,exon_id,gene_name])
    if strand == "+":
        all_exons = sorted(all_exons)
    elif strand == "-":
        all_exons = sorted(all_exons,reverse=True)
    for i,exon_st_end in enumerate(all_exons):
        i += 1
        exon_start = exon_st_end[0]
        exon_end = exon_st_end[1]
        exon_info = [chrom,transcript_source,"exon",exon_start,exon_end,
                     ".",strand,".",exon_attr % (i,i)]
        gtf_info.append(exon_info)
        if annotation == "protein_coding":
            check_cds = [x for x in cds_regions if (exon_start <= x[0]) and (exon_end >= x[1])]
            if check_cds:
                if len(check_cds) != 1:
                    print(region,te_type)
                    print(exon_st_end)
                    print(check_cds)
                assert len(check_cds) == 1
                check_cds = check_cds[0]
                cds_info = [chrom,transcript_source,"CDS",check_cds[0],check_cds[1],
                            ".",strand,cds_frame[check_cds],exon_attr % (i,i)]
                gtf_info.append(cds_info)
            for start_codon in start_codons:
                check_start = (exon_start <= start_codon[0]) and (exon_end >= start_codon[1])
                if check_start:
                    start_info = [chrom,transcript_source,"start_codon",start_codon[0],start_codon[1],
                                  ".",strand,start_codon_frame[start_codon],exon_attr % (i,i)]
                    gtf_info.append(start_info)
            for stop_codon in stop_codons:
                check_stop = (exon_start <= stop_codon[0]) and (exon_end >= stop_codon[1])
                if check_stop:
                    stop_info = [chrom,transcript_source,"stop_codon",stop_codon[0],stop_codon[1],
                                 ".",strand,stop_codon_frame[stop_codon],exon_attr % (i,i)]
                    gtf_info.append(stop_info)
    return gtf_info,IPA_info


def novelTranscriptToGtf(composite_annotation,
                         skip_annotation):
    novel_gtf_infos = []
    novel_IPA_infos = []
    # composite 
    for region,transcripts in composite_annotation.items():
        for transcript in transcripts:
            gtf_info,IPA_info = transGTFormat(region,
                                     transcript,
                                     "Composite")
            novel_gtf_infos.extend(gtf_info)
            novel_IPA_infos.extend(IPA_info)
    # skip
    for region,transcripts in skip_annotation.items():
        for transcript in transcripts:
            gtf_info,IPA_info = transGTFormat(region,
                                     transcript,
                                     "Skip")
            novel_gtf_infos.extend(gtf_info)
            novel_IPA_infos.extend(IPA_info)
    return novel_gtf_infos,novel_IPA_infos


def getAnnotateGeneInfos(gtf_file_path,novel_terminal_genes):
    annotated_gtf_infos = []
    with open(gtf_file_path) as hd:
        for x in hd:
            if "NM_" in x:
                chr_gene_id = ":".join(re.search('(chr.*?)\t.*?gene_id "(.*?)"',x).groups())
                if chr_gene_id in novel_terminal_genes:
                    x = x.strip()
                    annotated_gtf_infos.append(x.split("\t"))
    return annotated_gtf_infos


def typeTrans(iv_type):
    type_code = {
        "transcript":1,
        "exon":2,
        "3UTR":3,
        "CDS":4,
        "5UTR":5,
        "start_codon":6,
        "stop_codon":7
    }
    return type_code[iv_type]


def create_parser(name):
    p = argparse.ArgumentParser(prog=name,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                description='Build transcript from IPA.')
    # input
    g = p.add_argument_group('input')
    g.add_argument(
        '--finder_pa_file',
        type = str,
        help='Predicted terminal exon.')
    g.add_argument(
        '--annotated_gtf',
        type = str,
        help='Annotated GTF file.')
    g.add_argument(
        '--fa_path',
        type = str,
        help='Fasta file.')
    # output
    g = p.add_argument_group('output')
    g.add_argument(
        '--save_file',
        type = str,
        help='Save file.')
    return p

# run
# parse args
args = sys.argv
parser = create_parser(args[0])
args = parser.parse_args(args[1:])
# variable
finder_pa_file = args.finder_pa_file
RefSeq_gtf_file_path = args.annotated_gtf
fasta_path = args.fa_path
save_gtf = args.save_file
assert save_gtf.lower().endswith(".gtf")
save_ipa_info = save_gtf[:-4] + ".IPA_info.txt"
# Build transcript
# upstream site
total_info_df = extractUpInfo(finder_pa_file,RefSeq_gtf_file_path)
novel_terminal_genes,rt_te_skip,rt_te_composite = splitTerminal(total_info_df)
# extrat info of gtf
gtf_iv_df = ivFromGtf(RefSeq_gtf_file_path,novel_terminal_genes)
# read fasta
genome = Fasta(fasta_path)
# composite
composite_transcript_dic = compositeTranscript(rt_te_composite)
checkComTrans(composite_transcript_dic,
              genome)
composite_transcript_annotation_dic = determinCoding(composite_transcript_dic,
                                                     genome)
# skip
skip_transcript_dic = skipTranscript(rt_te_skip)
skip_transcript_annotation_dic = determinCoding(skip_transcript_dic,
                                                genome)
# novel gtf
novel_gtf_infos,novel_IPA_infos = novelTranscriptToGtf(composite_transcript_annotation_dic,
                                                       skip_transcript_annotation_dic)
# annotated gtf
annotated_gtf_infos = getAnnotateGeneInfos(RefSeq_gtf_file_path,novel_terminal_genes)
# save　GTF
annotated_novel_merge_infos = annotated_gtf_infos + novel_gtf_infos
annotated_novel_df = pd.DataFrame(sorted(annotated_novel_merge_infos,
                                         key = lambda x:[x[0],
                                                         re.search('gene_id "(.*?)"',x[8]).groups()[0],
                                                         re.search('transcript_id "(.*?)"',x[8]).groups()[0],
                                                         x[3],
                                                         typeTrans(x[2])]))
annotated_novel_df.to_csv(save_gtf,sep="\t",index=False,header=False,quoting=3)
log.info("Save " + save_gtf)

# save IPA info
novel_IPA_info_df = pd.DataFrame(novel_IPA_infos)
novel_IPA_info_df.columns = ["SYMBOL","Terminal_exon","IPAtype","IPUI","IPA_trans","IPA_trans_long"]
novel_IPA_info_df.to_csv(save_ipa_info,sep="\t",index=False,header=True)
log.info("Save " + save_ipa_info)


