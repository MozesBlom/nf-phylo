#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import argparse
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
parser.add_argument("-c", "--missing_data_cut_off", help="Float between 0 and 1, upper ratio of missing data for a given alignment column to be retained", default=0.5, action = "store", required=False)
parser.add_argument("-p", "--prefix_fn", help="Prefix to be used for the consensus files", action = "store", required=True)
parser.add_argument("-o", "--out_dir", help="Directory to store the MSA, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
chromo_aln_raw_fn = args.msa_path
prefix = args.prefix_fn
col_n_ratio = args.missing_data_cut_off
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Function
#########################


def filt_aln(msa, missing_data_cut_off):
	# Empty msa object (note it's an alignment because indel variation should have been removed so each cons should be the same length)
	chromo_msa_filt = MultipleSeqAlignment([])
	# Now we want to create a filtered concatenated genome where we filter out columns with too much missing data
	counter = 0
	while counter < msa.get_alignment_length():
		col = msa[:, counter]
		col_miss = col.count('N') + col.count('n') + col.count('-') 
		if (float(col_miss)/float(len(col))) < float(missing_data_cut_off):
			if (len(chromo_msa_filt) > 0):
				chromo_msa_filt = chromo_msa_filt + msa[:, counter:(counter + 1)]
			else:
				chromo_msa_filt = msa[:, counter:(counter + 1)]
		else:
			pass
		counter += 1	
	return chromo_msa_filt


#########################
## Automated Analysis
#########################

chromo_aln_raw = AlignIO.read(chromo_aln_raw_fn, 'fasta')
chromo_aln_filt = filt_aln(chromo_aln_raw, col_n_ratio)
chromo_aln_filt_fn = os.path.join(output_path, (prefix + '_filt.fa'))
AlignIO.write(chromo_aln_filt, chromo_aln_filt_fn, 'fasta')
