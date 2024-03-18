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

parser.add_argument("-m", "--msa_fn", help="MSA file --> path", action = "store", required=True)
parser.add_argument("-s", "--subset_len", help="Integer or value string with length of window", default=10000, type=int, action = "store", required=False)
parser.add_argument("-j", "--jump_len", help="Integer or value string with length BETWEEN windows", default=90000, type=int, action = "store", required=False)
parser.add_argument("-c", "--cut_off_ratio_col", help="Float between 0 and 1, upper ratio of missing data for a given alignment column to be considered as sufficient", default=0.5, action = "store", required=False)
parser.add_argument("-w", "--cut_off_ratio_ncol_per_window", help="Float between 0 and 1, minimum ratio of columns per window", default=0.5, action = "store", required=False)
parser.add_argument("-r", "--cut_off_ratio_row", help="Float between 0 and 1, upper ratio of missing data for a given consensus sequence to be considered as sufficient", default=0.2, action = "store", required=False)
parser.add_argument("-i", "--cut_off_ratio_nrow_per_window", help="Float between 0 and 1, minimum ratio of individuals per window", default=0.5, action = "store", required=False)
parser.add_argument("-p", "--prefix", help="String that will be the basename of the output file; i.e. chromo_0_1000000.fa", default='chromo', action = "store", required=False)
parser.add_argument("-o", "--out_dir", help="Directory to store the MSA, if unspecified assumes current work dir", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
msa = args.msa_fn
window_size = args.subset_len
jump_size = args.jump_len
n_col_max = args.cut_off_ratio_col
n_col_window_min = args.cut_off_ratio_ncol_per_window
n_row_max = args.cut_off_ratio_row
n_row_window_min = args.cut_off_ratio_nrow_per_window
prefix_fn = args.prefix
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()


#########################
## Function
#########################

# Read in alignment, evaluate, if window fulfils criteria create alignment
def eval_window(msa, start, end, window_len, coll_cut_off, ratio_min_window_size, row_cut_off, ratio_min_number_of_indivs, prefix):
	window = msa[:, start:end]
	## Create an empty alignment object where we will append alignment columns with sufficient data
	msa_tmp = []
	for record in window:
		msa_tmp.append(SeqRecord(Seq(''), id=record.id, annotations={"molecule_type": "DNA"}))
	window_col_subset = MultipleSeqAlignment(msa_tmp)
	col_counter = 0
	while (col_counter < window_len):
		col_aln = window[:, col_counter]
		col_n = col_aln.count('N') + col_aln.count('n') + col_aln.count('-')
		if (col_n <= coll_cut_off):
			# The number of missing genotypes for this column is lower than our cut-off: Append to empty MSA
			window_col_subset = window_col_subset + window[:, col_counter:(col_counter+1)]
		else:
			pass
		col_counter += 1
	## Columns with too much missing data have now been removed.
	## If the window is too short, terminate eval for this window. Otherwise continue:
	if float(window_col_subset.get_alignment_length()/window_len) < float(ratio_min_window_size):
		return False
	else:
		## Window is still long enough so let's eval missing data per individual
		msa_subset = []
		msa_all = []
		for indiv in window_col_subset:
			row_n = indiv.seq.count('N') + indiv.seq.count('n') + indiv.seq.count('-')
			if (row_n <= row_cut_off):
				## Sequence sufficient data, append to both the subset/filtered as complete msa
				msa_subset.append(indiv)
				msa_all.append(indiv)
			else:
				## Sequence doesn't have enough data. Still append to complete msa (needed for concordance factor annotation)
				msa_all.append(indiv)
		window_col_indiv_subset = MultipleSeqAlignment(msa_subset)
		window_col_indiv_all = MultipleSeqAlignment(msa_all)
		## If the number of indivs for this window is too low, terminate eval. Otherwise generate alignment
		if float(len(window_col_indiv_subset)/len(window)) < float(ratio_min_number_of_indivs):
			return False
		else:
			## Window passes all criteria, let's create a new sequence alignment and return True.
			## Note we create two alignments:
			## - 1) All individuals present (even if complete Ns)
			## - 2) Only individuals that surpass the filtering criteria
			window_indiv_subset_fn = prefix + '_' + str(start) + '_' + str(end) + '_filt_indivs_aln.fa'
			window_indiv_all_fn = prefix + '_' + str(start) + '_' + str(end) + '_all_indivs_aln.fa'
			AlignIO.write(window_col_indiv_subset, window_indiv_subset_fn, 'fasta')
			AlignIO.write(window_col_indiv_all, window_indiv_all_fn, 'fasta')
			return True


#########################
## Automated Analysis
#########################

## Read in alignment only once
msa_original = AlignIO.read(msa, 'fasta')
msa_len = msa_original.get_alignment_length()
# Max number of genotypes missing for an individual
max_row = int(float(window_size) * float(n_row_max))
# Max number of alignment cols missing per window
indiv_num = len(msa_original)
max_col = int(float(indiv_num) * float(n_col_max))
window_start = 0
window_end = int(window_size)

## Now iterate over chromosome, evaluate and filter window. If window suffices then jump, otherwise move to the adjacent region and try again
while (int(window_end) <= int(msa_len)):
	if eval_window(msa_original, window_start, window_end, window_size, max_col, n_col_window_min, max_row, n_row_window_min, prefix_fn):
		# Window suffices and an alignment has been generated
		window_start = window_end + int(jump_size)
		window_end = window_start + int(window_size)
	else:
		# Window does not pass one or both of the filtering critera!
		# Move to adjacent area and recalculate
		window_start = window_end
		window_end = window_start + int(window_size)





