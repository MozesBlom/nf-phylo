#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
	import os
	import argparse
	from natsort import natsorted
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

parser.add_argument("-m", "--msa_fn_list", help="List of MSA files --> paths", action = "store", nargs='+', default=[], required=True)
parser.add_argument("-p", "--prefix", help="String that will be the basename of the concatenated alignment file; i.e. prefix.fa", default='genome_msa', action = "store", required=False)
parser.add_argument("-o", "--out_dir", help="Directory to store the MSA, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
msa_list = args.msa_fn_list
prefix_fn = args.prefix
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Function
#########################

def concat_aln(msa_path):
	# Empty msa object (note it's an alignment because indel variation should have been removed so each cons should be the same length)
	genome_msa = []
	# Iterate over chromo MSAs and append to msa object
	for msa_fn in msa_path:
		aln = AlignIO.read(msa_fn, 'fasta')
		aln.sort()
		# First let's just concat the raw msa's
		if (len(genome_msa) > 0):
			genome_msa = genome_msa + aln
		else:
			genome_msa = aln
	return genome_msa

#########################
## Automated Analysis
#########################

## First 'natural sort' the alignment list. This is crucial to stitch the subset filtered alignments back together in the right order. Less important for whole genome concat.
## However, one issue is that it will sort by path and not by filename. We therefore first need to sort by filename and then match the paths
msa_fns = []
for msa in msa_list:
	msa_fns.append(os.path.basename(msa))
## Natural sort filenames
msa_fns_sorted = natsorted(msa_fns)


## Now natural sort path based on filenames
sorted_msa_list = []
for msa_fn in msa_fns_sorted:
	for msa_path in msa_list:
		if msa_path.endswith(msa_fn):
			sorted_msa_list.append(msa_path)
		else:
			pass


## Create an output file with the order of MSA paths, just as a sanity check
with open(os.path.join(output_path, (prefix_fn + '_concat_order.txt')), 'w') as msa_concat_order_fn:
	for msa in sorted_msa_list:
		msa_concat_order_fn.write(msa + '\n')


## Create concatenated alignment
genome_aln = concat_aln(sorted_msa_list)
genome_aln_fn = os.path.join(output_path, (prefix_fn + '.fa'))
AlignIO.write(genome_aln, genome_aln_fn, 'fasta')

