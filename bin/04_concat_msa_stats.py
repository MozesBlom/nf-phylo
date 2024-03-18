#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import argparse
	import pandas as pd
	from Bio.Seq import *
	from Bio import SeqIO
	from Bio import AlignIO
	from Bio.Align import MultipleSeqAlignment
	from Bio.SeqRecord import SeqRecord
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--msa_path", help="Multiple sequence alignment file --> path", action = "store", required=True)
parser.add_argument("-p", "--prefix", help="String that will be the basename of the msa stats file; i.e. [prefix]_missing_data_per_indiv.tsv", default='genome_msa', action = "store", required=False)
parser.add_argument("-o", "--out_dir", help="Directory to store the MSA, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
msa_fn = args.msa_path
prefix_fn = args.prefix
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Function
#########################

def indiv_stats(msa, df_stats):
	## Loop over genome, gather stats on missing data by indiv and append to pandas df
	msa_len = msa.get_alignment_length()
	for record in msa:
		indiv = record.id
		indiv_n = record.seq.count('N') + record.seq.count('n') + record.seq.count('-')
		indiv_rel = (float(indiv_n)/float(msa_len)) * 100
		indiv_df = pd.DataFrame([[indiv, indiv_n, indiv_rel]], columns=['Indiv','Missing_data_VAL','Missing_data_PERCENTAGE'])
		df_stats = pd.concat([df_stats, indiv_df], ignore_index=True)
	df_stats = df_stats.round({'Missing_data_PERCENTAGE': 2})
	return df_stats

#########################
## Automated Analysis
#########################

## Create a pandas dataframe with a summary of all the read stats
df = pd.DataFrame(columns=['Indiv','Missing_data_VAL','Missing_data_PERCENTAGE'])
msa = AlignIO.read(msa_fn, 'fasta')
df_stats = indiv_stats(msa, df)
tsv_fn = os.path.join(output_path, (prefix_fn + '_missing_data_per_indiv.tsv'))
df_stats.to_csv(tsv_fn, sep='\t', index=False)

