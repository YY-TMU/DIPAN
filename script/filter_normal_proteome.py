import sys
import logging
import argparse
import pandas as pd

# set log
logging.basicConfig(format='%(asctime)s: %(message)s',
                    datefmt = '%Y-%m-%d %H:%M')
log = logging.getLogger()
log.setLevel(logging.INFO)

def create_parser(name):
    p = argparse.ArgumentParser(prog=name,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                description='Filter peptides in normal proteome.')
    # input
    g = p.add_argument_group('input')
    g.add_argument(
        '--normal_proteome_fa',
        type = str,
        help='Fasta file from normal proteome.')
    g.add_argument(
        '--peptide_file',
        type = str,
        help='Results of netMHC.')
    # output
    g = p.add_argument_group('output')
    g.add_argument(
        '--save_file',
        type = str,
        help='Output file.')
    return p


# run
# parse args
args = sys.argv
parser = create_parser(args[0])
args = parser.parse_args(args[1:])
# variable
normal_proteome_pt = args.normal_proteome_fa
mhc_filter_result = args.peptide_file
save_file = args.save_file
# extract normal sequences
normal_aa_seq_dic = {}
with open(normal_proteome_pt,"r") as fa:
    for x in fa.readlines():
        if x.startswith(">"):
            name = x.strip()
        else:
            normal_aa_seq_dic.setdefault(name,[]).append(x.strip())
normal_aa_seq_dic = {k:"".join(v).upper() for k,v in normal_aa_seq_dic.items()}
normal_aa_seqs = list(normal_aa_seq_dic.values())
merge_normal_aa_seqs = "#".join(normal_aa_seqs)
# read bind sequences
mhc_result_df = pd.read_table(mhc_filter_result,delim_whitespace=True)
bind_sequences = set(mhc_result_df["Peptide"])
log.info("Number of total peptides : %d" % len(bind_sequences))
# filter peptides in normal
same_sequences_with_normal = set([bind_sequence for bind_sequence in bind_sequences if bind_sequence in merge_normal_aa_seqs])
log.info("Number of same peptides with normal : %d" % len(same_sequences_with_normal))
unique_sequences = bind_sequences.difference(same_sequences_with_normal)
unique_mhc_result_df = mhc_result_df.loc[mhc_result_df["Peptide"].isin(unique_sequences),]
unique_mhc_result_df.to_csv(save_file,sep="\t",index=False,header=True)
log.info("Save %s" % save_file)

