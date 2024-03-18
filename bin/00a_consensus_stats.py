#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import argparse
	import pandas as pd
	from Bio import SeqIO
except ImportError:
	sys.exit("One of the required modules can't be found...")


#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cons_fn_list", help="List of consensus files --> paths", action = "store", nargs='+', default=[], required=True)
parser.add_argument("-i", "--individual", help="Individual to which consensus sequences belong", required=True)
parser.add_argument("-o", "--out_dir", help="Directory to store the MSA, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
consensus_list = args.cons_fn_list
indiv = args.individual
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Automated Analysis
#########################

## Create a pandas dataframe with a summary of all the read stats
df = pd.DataFrame(columns=['Chromo','Length', 'Missing_data', 'Missing_data_PERCENTAGE'])

## Loop over consensus seqs and estimate missing data
for consensus in consensus_list:
	chromo = os.path.basename(consensus).split('_')[1]
	record = SeqIO.read(consensus, "fasta")
	cons_len = len(record.seq)
	cons_n = record.seq.count('N') + record.seq.count('n') + record.seq.count('-')
	cons_n_rel = (float(cons_n)/float(cons_len)) * 100
	cons_df = pd.DataFrame([[chromo, cons_len, cons_n, cons_n_rel]], columns=['Chromo','Length', 'Missing_data', 'Missing_data_PERCENTAGE'])
	df = pd.concat([df, cons_df], ignore_index=True)

## Calculate total missing data
total_len = df['Length'].sum()
total_n = df['Missing_data'].sum()
total_n_rel = (float(total_n)/float(total_len)) * 100
total_df = pd.DataFrame([['Total', total_len, total_n, total_n_rel]], columns=['Chromo','Length', 'Missing_data', 'Missing_data_PERCENTAGE'])
df = pd.concat([df, total_df], ignore_index=True)
## Round off by two decimals
df = df.round({'Missing_data_PERCENTAGE': 2})

## Write out file
tsv_fn = os.path.join(output_path, (indiv + '_missing_data.tsv'))
df.to_csv(tsv_fn, sep='\t', index=False)



