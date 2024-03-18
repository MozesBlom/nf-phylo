#!/usr/bin/env nextflow

/*
 * Specify that code written in DSL1
 */
nextflow.enable.dsl=1

log.info """\
		From VCF to Phylogeny - NF phylo    
		===================================
		indivfile    : ${params.indivs_file}
		chromos      : ${params.chromos_file}
		reference    : ${params.ref_file}
		vcfdir       : ${params.inputdir}
		outdir       : ${params.outputdir}
		ASTRAL       : ${params.astral_fn}

		Window sizes : ${params.window_sizes}
		Sex chromos  : ${params.sex_chromos}
		"""
		.stripIndent()

/*
 * Create two lists:
 * - All individuals to include
 * - All chromosomes to include
 */

indivs = file(params.indivs_file).readLines()
chromos = file(params.chromos_file).readLines()


process call_consensus {
    
/*
 * Call consensus sequences for each individual and each chromosome.
 *
 * Directed into two channels that are ultimately used for:
 * ch1: process(chromo_raw_msa) -- Create raw chromosome alignments
 * ch2: process(calc_cons_N) -- Calculate missing data per consensus
 *
 */

	tag "Create consensus sequence for $chromo of $indiv"
	publishDir "${params.outputdir}/00.consensus", mode:'copy'

	input:
	val(indiv) from indivs
	each chromo from chromos

	output:
	tuple val(indiv), \
	val(chromo), \
	path("${indiv}_${chromo}_cons.fa") into consensus_ch1, consensus_ch2

	when:
	params.call_consensus == true

	script:
	def chromo_var_fn = file("${params.inputdir}/${indiv}/${chromo}_vars_filt_indels.vcf.gz")
	def chromo_mask_fn = file("${params.inputdir}/${indiv}/${chromo}_mask_cov_het.vcf")
	"""

	samtools faidx ${params.ref_file} ${chromo} | bcftools consensus ${chromo_var_fn} -m ${chromo_mask_fn} -o ${indiv}_${chromo}_cons.fa

	sed -i 's/${chromo}/${indiv}/g' ${indiv}_${chromo}_cons.fa

	"""
}


/*
 * To create MSA per chromosome/scaffold, a channel is needed where all consensus sequences are grouped by chromo/scaffold.
 * Consensus calling done internally or as external input?
 *
 * ch1: For downstream analyses: Creating RAW multiple sequence alignments (chromo_raw_msa)
 *
 */

if (params.call_consensus == true) {
	consensus_ch1
		.groupTuple(by: 1)
		.map { indiv, chromo, cons_fn -> tuple( chromo, indiv.sort{it}, cons_fn.sort{it} ) }
		.set { consensus_by_chromo_ch1 }
} else {
	channel
		.fromPath(params.consensus_tsv_fn)
		.splitCsv(header:true)
		.map { row -> tuple(row.Individual, row.Scaffold, file(row.Consensus_fn)) }
		.groupTuple(by: 1)
		.map { indiv, chromo, cons_fn -> tuple( chromo, indiv.sort{it}, cons_fn.sort{it} ) }
		.set { consensus_by_chromo_ch1 }
}


/*
 * To estimate the amount of missing data, a channel is needed where all consensus sequences are grouped by individual.
 * Consensus calling done internally or as external input?
 *
 * ch1: Calculating missing data per individual
 *
 */

if (params.call_consensus == true) {
	consensus_ch2
		.groupTuple(by: 0)
		.map{ it ->

			def indiv = it[0]
			def cons_fn = it[2].flatten()

			[indiv, cons_fn]
		}
		.set { consensus_by_indiv_ch1 }
} else {
	channel
		.fromPath(params.consensus_tsv_fn)
		.splitCsv(header:true, sep:'\t')
		.map { row -> tuple(row.Individual, row.Scaffold, file(row.Consensus_fn)) }
		.groupTuple(by: 0)
		.map{ it ->

			def indiv = it[0]
			def cons_fn = it[2].flatten()

			[indiv, cons_fn]
		}
		.set { consensus_by_indiv_ch1 }
}


process calc_cons_N {

	/*
	 * Estimate the amount of missing data per individual
	 *
	 * Two stats file:
	 * 1) Amount of missing data for each individual and each chromosome
	 * 2) Overall amount of missing data for each individual (genome wide)
	 */

	tag "Estimate missing data per individual"
	publishDir "${params.outputdir}/00.consensus/00.stats/", mode:'copy'

	input:
	tuple val(indiv), path(cons_fn_list) from consensus_by_indiv_ch1

	output:
	file("${indiv}_missing_data.tsv") into consensus_by_indiv_stats_ch1

	script:
	"""

	00a_consensus_stats.py \
	-c ${cons_fn_list} \
	-i ${indiv}

	"""
}


/*
 * Wait until all summary stats have been calculated
 */

consensus_by_indiv_stats_ch1
	.collect ()
	.set { consensus_by_indiv_stats_ch2 }


process summarise_cons_N {

	/*
	 * Create summary overview of missing data
	 */

	tag "Summarise missing data"
	publishDir "${params.outputdir}/00.consensus/00.stats/", mode:'copy', overwrite:'false'

	input:
	path(cons_stats_fn_list) from consensus_by_indiv_stats_ch2

	output:
	file("all_indivs_missing_data.tsv")
	file("all_indivs_missing_data.png")

	script:
	"""
	unset DISPLAY
	
	00b_consensus_stats.py \
	-c ${cons_stats_fn_list}

	"""
}


process chromo_raw_msa {

	/*
	 * Create multiple sequence alignment for each chromosome
	 * Directed into five channels
	 * 1) For process (subset_MSA_to_windows) -- Subsetting each chromosome into windows
	 * 2) For process (chromo_subset_msa) -- Filtering complete chromos in chunks
	 * 3) For process (concat_raw) -- Concatenation of unfiltered autosomes (first channel processed to exclude sex chromos)
	 * 4) For process (infer_chromo_raw_phy) -- Infer ML phylogeny for each chromosome
	 */

	tag "Raw multiple sequence alignment per chromosome -- raw"
	publishDir "${params.outputdir}/01.alignments/complete/00.chromosomes/00.raw", mode:'copy'

	input:
	tuple val(chromo), \
	val(indivs), \
	path(cons_fn) from consensus_by_chromo_ch1

	output:
	tuple val(chromo), path("${chromo}_raw_msa.fa") into msa_by_chromo_raw_ch1, \
	msa_by_chromo_raw_ch2, \
	msa_by_chromo_raw_ch3, \
	msa_by_chromo_raw_ch4

	script:
	"""

	cat ${cons_fn} > ${chromo}_raw_msa.fa

	"""
}


/*
 * Consensus sequences have now been called and evaluated, and raw sequence alignments have been created for each chromosome.
 *
 * We will now subset the genome, for a range of subset sizes, to facilitate window based analyses.
 * The subsetting is done before the concatenated analyses because for some datasets we don't need entire genomes.
 * A subset of the genome is fine as well for concatenation ML trees etc.
 *
 */


process subset_MSA_to_windows {

	/*
	 * Subset MSA for each chromosome in windows of size X and minimally spaced Y bp. from eachother
	 * Each candidate window is subsequently inspected for missing data: a) Per column and b) per individual
	 * a) Per column: params.missing_data_per_col_threshold, params.missing_cols_per_window_threshold
	 * b) Per individual: params.missing_data_per_row_threshold, params.missing_indivs_per_window_threshold
	 *
	 * If a given window surpasses the filtering criteria specified above. It will generate two alignments for the same window:
	 * "*all_indivs_aln.fa" = This includes all individuals, even if some individuals may not fulfill (missing_data_per_row_threshold)
	 * "*filt_indivs_aln.fa" = This includes only those individuals that fulfill the (missing_data_per_row_threshold) threshold
	 *
	 * The reason that both are generated is because IQtree will not infer a gene tree for a window where a certain individual only has N's
	 * However, for concatenation (and concordance factor annotation) each subset alignment must have the same number of individuals.
	 * If the dataset is of really high quality and the criteria above set loose then both alignments should be effectively the same.
	 *
	 * Directed into three channels:
	 * ch1: process(window_msa_to_phy) -- ML window trees for each window
	 * ch2: process(concat_autosomal_window_MSA_by_subset) -- To create a concatenated MSA of all autosomal windows (for each subset size)
	 * ch3: process(concat_chromo_window_MSA_by_subset) -- To create a concatenated MSA for each chromosome (for each subset size)
	 *
	 * NOTE: THIS SUBSETTING IS ALWAYS DONE BECAUSE IT IS NEEDED FOR THE SUMMARY-COALESCENT ST INFERENCE.
	 * In example, the params.infer_ML_on_subsets refers to whether ML trees are inferred for complete chromosomes/genomes or subsets.
	 */

	label 'SC_SM_HT'
    
	tag "Create window alignments for each chromo"
	publishDir "${params.outputdir}/01.alignments/subset/02.windows/${subset}/${chromo}", mode:'copy'

	input:
	tuple val(chromo), path(chromo_msa) from msa_by_chromo_raw_ch1
	each subset from params.window_sizes

	output:
	tuple val(chromo), val(subset), file("*filt_indivs_aln.fa") optional true into window_msa_by_subset_ch1
	tuple val(chromo), val(subset), file("*all_indivs_aln.fa") optional true into window_msa_by_subset_ch2, window_msa_by_subset_ch3

	script:
	"""

	01_window_subset_MSA.py \
	-m ${chromo_msa} \
	-s ${subset} \
	-j ${params.window_jump} \
	-c ${params.missing_data_per_col_threshold} \
	-w ${params.missing_cols_per_window_threshold} \
	-r ${params.missing_data_per_row_threshold} \
	-i ${params.missing_indivs_per_window_threshold} \
	-p ${chromo}

	"""
}


/* 
 * For each subset size, concatenate subset window alignments for each chromosome and all autosomes.
 *
 * First differentiate between sex and autosomes, then group by subset and direct into channel with a list of
 * all window alignments (per subset).
 *
 * NOTE, using GroupTuple results in double brackets for the MSA list (i.e.: [[aln1, aln2]]) and therefore would
 * impact the .size estimate if not taken care off. After GroupTuple, I therefore included msa_list.flatten in the
 * .map operator. This removes the double brackets
 *
 * IMPORTANT: The concatenation of the windows is not necessarily in the original consecutive order (i.e. by chromo position)!!
 */

window_msa_by_subset_ch2
	.branch { 
		auto : it[0] !in params.sex_chromos 
		sex : it[0] in params.sex_chromos 
	}
	.set { window_msa_by_subset_chromo_type_ch1 }

// Here all autosomes will be combined so it doesn't really matter what chromo a given alignment is from
// We can just use GroupTuple on the subset only and don't have to include the chromo information
window_msa_by_subset_chromo_type_ch1.auto
	.groupTuple(by: 1)
	.map{ it ->

		def chromo = it[0]
		def subset = it[1]
		def msa_list = it[2].flatten()

		[subset, msa_list, msa_list.size()]
	}
	.set { auto_window_msa_by_subset_ch1 }


// This channel is created so alignments that belong to the same chromosome can be concatenated.
// This is across all subsets and GroupTuple should therefore be employed on both subset and chromosome
window_msa_by_subset_ch3
	.groupTuple(by: [0,1])
	.map{ it ->

		def chromo = it[0]
		def subset = it[1]
		def msa_list = it[2].flatten()

		[chromo, subset, msa_list, msa_list.size()]
	}
	.into { window_msa_by_subset_ch4; window_msa_by_subset_ch5  }

// Included for troubleshooting
window_msa_by_subset_ch5
	.view()


process concat_autosomal_window_MSA_by_subset {

	/*
	 * Concatenate the windows for all autosomes. Process only executed if more than 1 window available: window_msa_number > 1
	 *
	 * We now have a MSA for each window and a list of ALL autosomal windows for a given subset
	 * Let's concatenate all windows from autosomes, which is needed to calculate CF using IQtree.
	 * NOTE: THE CONCATENATION DOES NOT NECESSARILY HAPPEN IN THE SEQUENTIAL ORDER OF CHROMOSOME POSITION!
	 * Mostly because it's computationally intensive and not needed for downstream analyses
	 *
	 * Directed into three channels:
	 * - ch1: Process (infer_auto_concat_subset_phy) -- To infer a ML tree based on concatenation of all autosomal windows (OPTIONAL, params.infer_ML_on_subsets = true)
	 * - ch2: Process (windowCF_sCF_auto_subset_CONCAT) -- To calculate sCF/wCF on a ML tree based on concatenation of all autosomal windows (OPTIONAL, params.infer_ML_on_subsets = true)
	 * - ch3: Process (windowCF_sCF_auto_subset_SUMCOAL) -- To calculate sCF/wCF on a summary-coalescent tree based on different sized autosomal windows
	 */

	label 'SC_IM_ST'

	tag "Concatenate all autosomal windows - per subset"
	publishDir "${params.outputdir}/01.alignments/subset/01.autosomes/${subset}", mode:'copy', overwrite:'false'

	input:
	tuple val(subset), path(window_msa_list), val(window_msa_number) from auto_window_msa_by_subset_ch1

	output:
	tuple val(subset), path("${subset}bp_autosomal_windows_concat_msa.fa") into subset_window_concat_autosome_msa_ch1, \
	subset_window_concat_autosome_msa_ch2, \
	subset_window_concat_autosome_msa_ch3

	when:
	window_msa_number > 1

	script:
	"""

	seqkit concat \
	--threads 2 \
	--seq-type dna \
	--out-file ${subset}bp_autosomal_windows_concat_msa.fa \
	${window_msa_list}

	"""
}


process concat_chromo_window_MSA_by_subset {

	/*
	 * Concatenate the windows per chromosome and subset. Process only executed if more than 1 window available: window_msa_number > 1
	 *
	 * Directed into two channels:
	 * - ch1: process (infer_chromo_concat_subset_phy). Infer a phylogeny for each chromosomes (based on concatenated subsets of that chromosome)
	 * - ch2: process (windowCF_sCF_sexchromos). CF annotation of the sex chromosomes
	 */

	label 'SC_IM_ST'

	tag "Concatenate all windows by chromosome - per subset"
	publishDir "${params.outputdir}/01.alignments/subset/00.chromosomes/${subset}", mode:'copy', overwrite:'false'

	input:
	tuple val(chromo), val(subset), path(window_msa_list), val(window_msa_number) from window_msa_by_subset_ch4

	output:
	tuple val(chromo), val(subset), path("${chromo}_filt_msa.fa") into subset_window_concat_chromo_msa_ch1, \
	subset_window_concat_chromo_msa_ch2

	when:
	window_msa_number > 1

	script:
	"""
	
	seqkit concat \
	--threads 2 \
	--seq-type dna \
	--out-file ${chromo}_filt_msa.fa \
	${window_msa_list}

	"""
}


/* 
 * All alignments have been generated for the analyses based on subsets of the genome.
 *
 * Do we also want to generate alignments for complete autosomes?
 * This can be specified by setting: params.infer_ML_on_subsets = false
 * If so, do they also need to be filtered in a similar way as the windows? (params.filter_complete_chromos = true)
 * Filtering for small windows is essential because poor datasets may result in non-informative windows
 */


/*
 * Filtering the MSA based on missing data can be painstakingly slow using biopython.
 * Therefore, if need be, first subset the MSA for each chromo, filter and then later on concat them together
 * NOTE: These intermediate files do not need to be stored, so there's no publishdir call for each process.
 */


process chromo_subset_msa {

	/*
	 * Create subsets of length x for each chromo (only if whole chromosomes are needed and they need to be filtered)
	 * Requires a BED file (params.chromo_info) where the name, start and end position of each scaffold are listed
	 *
	 * Directed into one channel
	 * 1) Process (chromo_filt_subset_msa). To filter the MSA for each subset.
	 */

	tag "Subset chromosome alignment for filtering"

	input:
	tuple val(chromo), \
	path(raw_msa_fn) from msa_by_chromo_raw_ch2

	output:
	tuple val(chromo), file("*_msa.fa") into subset_msa_list_by_chromo_ch1

	when:
	params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	"""
	#!/usr/bin/env python3

	try:
		import os
		import subprocess
		import pandas as pd
	except ImportError:
		sys.exit("One of the required modules can't be found...")

	df = pd.read_csv("${params.chromo_info}",delimiter='\t')
	msa_len = df.loc[df['chromo'] == '${chromo}', 'end'].iloc[0]

	start = 1
	end = int(${params.filt_subset_size}) * 1000000

	while end < int(msa_len):
		subprocess.call("seqkit subseq -r %s:%s -o ${chromo}_%s_%s_msa.fa ${raw_msa_fn}" % (str(start), str(end), str(start), str(end)), shell=True)
		start = start + (int(${params.filt_subset_size}) * 1000000)
		end = end + (int(${params.filt_subset_size}) * 1000000)

	subprocess.call("seqkit subseq -r %s:%s -o ${chromo}_%s_%s_msa.fa ${raw_msa_fn}" % (str(start), str(msa_len), str(start), str(msa_len)), shell=True)

	"""
}


/*
 * subset_msa_list_by_chromo_ch1 is of the format: [chromo, [chr1_0_100_msa.fa, chr1_100_200_msa.fa, ...]]
 * to filter each alignment in parallel, the channel needs to be transposed so it emits [chromo, msa.fa, msa_basename]:
 * [chromo, chr1_0_100_msa.fa, chr1_0_100_msa]
 * [chromo, chr1_100_200_msa.fa, chr1_100_200_msa ]
 * ...
 *
 * NOTE: Even IF the channel will not be used, we still create an empty ELSE channel because otherwise nextflow spits out an error
 * 		 in any of the processes that call for this channel (even if the When statement isn't activated)
 */

if (params.infer_ML_on_subsets == false && params.filter_complete_chromos == true) {
	subset_msa_list_by_chromo_ch1
		.transpose()
		.map{ it ->

			def chromo = it[0]
			def subset_msa_path = it[1]
			def basename = file(subset_msa_path).getSimpleName()

			[chromo, subset_msa_path, basename]
		}
		.set { subset_msa_by_chromo_ch1 }
} else {
	subset_msa_by_chromo_ch1 = Channel.empty()
}

process chromo_filt_subset_msa {

	/*
	 * Filter each MSA column based on a specific criteria
	 * E.g. if params.missing_data_threshold == 0.5. Then 50% of the indivs in MSA need a GT call
	 * for a given column. Otherwise the column is removed
	 *
	 * Directed into one channel
	 * 1) Process (stitch_filt_subset_msa_by_chromo). To stitch the filtered subsets back together to chromosomes.
	 *
	 * NOTE: Here filtering is only done by alignment column. Mostly because filtering is computationally really intensive
	 * and the filtering is done based on an entire chromosome so the amount of missing data per individual is likely 
	 * much less of an issue (otherwise individual shouldn't be included!). 
	 */

	label 'SC_SM_IT'

	tag "Filtering of each chromosome subset alignment"

	input:
	tuple val(chromo), \
	path(raw_msa_fn), \
	val(prefix) from subset_msa_by_chromo_ch1

	output:
	tuple val(chromo), \
	path("${prefix}_filt.fa") optional true into filt_subset_msa_by_chromo_ch1

	when:
	params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	"""

	02_chromo_filt_msa.py \
	-m ${raw_msa_fn} \
	-c ${params.missing_data_per_col_threshold} \
	-p ${prefix}

	"""
}


/*
 * filt_subset_msa_by_chromo_ch1 needs to be grouped, so all subsets for that chromo can be stitched back together
 *
 * NOTE: it's also important that concatenation happens in the right order!
 */

if (params.infer_ML_on_subsets == false && params.filter_complete_chromos == true) {
	filt_subset_msa_by_chromo_ch1
		.groupTuple(by: 0)
		.map{ it ->

			def chromo = it[0]
			def msa_list = it[1].flatten()

			[chromo, msa_list]
		}
		.set { filt_subset_msa_list_by_chromo_ch1 }
} else {
	filt_subset_msa_list_by_chromo_ch1 = Channel.empty()
}


/*
 * Here the order of concatenation is really important and a custom concatenation script is used that will use a natural sort
 * to ensure the order of concatenation is taken care off (this is much slower that seqkit though...)!
 */

process stitch_filt_subset_msa_by_chromo {

	/*
	 * Create multiple sequence alignment for each chromosome
	 * Directed into two channels
	 * 1) For process (concat_filt) -- Concatenation of all filtered autosomes
	 * 2) For process (infer_chromo_filt_phy) -- Infer ML phylogeny for each chromosome
	 */


	label 'SC_SM_IT'

	tag "Stitch the subsets for each chromosome back together"
	publishDir "${params.outputdir}/01.alignments/complete/00.chromosomes/01.filtered", mode:'copy'

	input:
	tuple val(chromo), path(subset_msa_list) from filt_subset_msa_list_by_chromo_ch1

	output:
	tuple val(chromo), path("${chromo}_filt_msa.fa") into msa_by_chromo_filt_ch1, msa_by_chromo_filt_ch2
	file("${chromo}_filt_msa_concat_order.txt")

	when:
	params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	"""

	03_concat_msa.py \
	-m ${subset_msa_list} \
	-p ${chromo}_filt_msa

	"""
}

/*
 * We now have filtered and raw MSA's for each chromosome. 
 *
 * We want to concatenate all autosomes together and keep the sex chromosome alignments separate.
 * The sex chromosome alignments will be used later for sCF annotation of the sex-chromosomes.
 *
 * NOTE: ML inference on the raw autosome alignment is always done. Therefore only (params.infer_ML_on_subsets == false)
 */

if (params.infer_ML_on_subsets == false) {
	msa_by_chromo_raw_ch3
		.branch { 
			auto : it[0] !in params.sex_chromos 
			sex : it[0] in params.sex_chromos 
		}
		.set { msa_by_chromo_raw_chromo_type_ch1 }
} else {
}

/*
 * Note that to concatenate, all previous processes (i.e. MSA generation) need to be finalised
 * We therefore use the collect operator which should wait until all input processes have finalised
 * collect then emits all chromosome alignments as a single list.
 */
if (params.infer_ML_on_subsets == false) {
	msa_by_chromo_raw_chromo_type_ch1.auto
		.map { chromo, chromo_msa_fn -> chromo_msa_fn }
		.collect ()
		.set { msa_auto_raw_ch }
} else {
	msa_auto_raw_ch = Channel.empty()
}

/*
 * Repeat for the filtered alignments, if filtering was done
 */

if (params.infer_ML_on_subsets == false && params.filter_complete_chromos == true) {
	msa_by_chromo_filt_ch1
		.branch { 
			auto : it[0] !in params.sex_chromos 
			sex : it[0] in params.sex_chromos 
		}
		.set { msa_by_chromo_filt_chromo_type_ch1 }
} else {
}

if (params.infer_ML_on_subsets == false && params.filter_complete_chromos == true) {
	msa_by_chromo_filt_chromo_type_ch1.auto
		.map { chromo, chromo_msa_fn -> chromo_msa_fn }
		.collect ()
		.set { msa_auto_filt_ch }
} else {
	msa_auto_filt_ch = Channel.empty()
}


/*
 * For the full autosome concatenation the order is less of importance and we can therefore just take all of the entries
 * in the channel. We can then use Seqkit for a rapid concatenation
 */

 process concat_raw {

	/*
	 * We now have multiple sequence alignments for each chromosome
	 * Let's concatenate those in a whole genome alignment.
	 *
	 * Directed into three channels
	 * 1) process(infer_auto_concat_raw_phy). To infer ML phylogeny based on concatenation of complete autosomes
	 * 2) process(windowCF_sCF_auto_complete_raw_CONCAT). Calculate sCF/wCF on concatenated ML phylogeny -- raw autosomes
	 */

	label 'SC_HM_IT'

	tag "Concatenate all autosomes - raw"
	publishDir "${params.outputdir}/01.alignments/complete/01.autosomes/00.raw", mode:'copy', overwrite:'false'

	input:
	path(msa_fn) from msa_auto_raw_ch

	output:
	path("autosomes_concat_raw_msa.fa") into concat_msa_auto_raw_ch1, \
	concat_msa_auto_raw_ch2

	when:
	params.infer_ML_on_subsets == false

	script:
	"""

	seqkit concat \
	--threads 2 \
	--seq-type dna \
	--out-file autosomes_concat_raw_msa.fa \
	${msa_fn}

	"""
}


process concat_filt {

	/*
	 * We now have filtered multiple sequence alignments for each chromosome, sex chromosomes excluded
	 * Let's concatenate those in a whole genome alignment.
	 *
	 * Directed into three channels
	 * 1) process(infer_auto_concat_filt_phy). To infer a concatenated ML phylogeny 
	 * 2) process(windowCF_sCF_auto_complete_filt_CONCAT). Calculate sCF/wCF on concatenated ML phylogeny -- filtered autosomes
	 * 3) 
	 */

	label 'SC_HM_IT'

	tag "Concatenate all autosomes - filt"
	publishDir "${params.outputdir}/01.alignments/complete/01.autosomes/01.filtered", mode:'copy', overwrite:'false'

	input:
	path(msa_fn) from msa_auto_filt_ch

	output:
	path("autosomes_concat_filt_msa.fa") into concat_msa_auto_filt_ch1, \
	concat_msa_auto_filt_ch2

	when:
	params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	"""

	seqkit concat \
	--threads 2 \
	--seq-type dna \
	--out-file autosomes_concat_filt_msa.fa \
	${msa_fn}

	"""
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 		TIME FOR TREE BUILDING!! 												   *
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 																				   *
 * 		- Maximum likelihood phylogeny for each window in each subset
 * 		- Maximum likelihood phylogeny based on all autosomes combined 			   *
 * 		- Maximum likelihood phylogeny for each chromosome 						
 * 		- Summary-coalescent phylogeny based on window trees from all autosomes    *
 *
 * 																				   *
 ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		MAXIMUM-LIKELIHOOD PHYLOGENY FOR EACH WINDOW 	   *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		NOTE, the 01_window_subset_MSA.py script emits a 
 *		single list of MSA per chromo/subset and to 
 *		parallelize the process, a single item for EACH 
 *		MSA is needed. In example, with the following info:
 * 		[chromo, subset, aln_path, basename of alignment]
 * 
 */


window_msa_by_subset_ch1
	.transpose()
	.map{ it ->

		def chromo = it[0]
		def subset = it[1]
		def aln_path = it[2]
		def basename = file(aln_path).getSimpleName()

		[chromo, subset, aln_path, basename]
	}
	.set { window_msa_by_subset_annotated_ch1 }


process window_msa_to_phy {

	/*
	 * Take each window MSA and infer phylogeny using IQtree2
	 *
	 * Directed into one channel:
	 * - ch1: process (gather_subset_autosomal_window_trees) -- Gather all autosomal trees for ASTRAL 
	 */
    
	label 'SC_SM_IT'

	tag "Create window phylogenies for all windows (across all subsets)"
	publishDir "${params.outputdir}/02.phylogenies/subset/02.windows/${subset}/${chromo}", mode:'copy'

	input:
	tuple val(chromo), val(subset), path(msa_fn), val(basename) from window_msa_by_subset_annotated_ch1

	output:
	tuple val(chromo), val(subset), file("${basename}.treefile") into window_phy_by_subset_ch1

	script:
	"""

	iqtree \
	-s ${msa_fn}  \
	-m MFP \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix ${basename}

	"""
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		MAXIMUM-LIKELIHOOD PHYLOGENY BASED ON ALL AUTOSOMES COMBINED 	   *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		ML phylogeny based on a concatenated alignment of all autosomes
 *		
 *		This can either be on complete autosomes or on a concatenated
 *		dataset of autosomal subsets
 * 		
 * 		params.infer_ML_on_subsets = true/false
 *
 *		IF params.infer_ML_on_subsets == false :
 *		The ML phylogeny will be based on complete autosomes.
 *		Is there also need for a ML phylogeny based on filtered alignment?
 *		
 *		params.filter_complete_chromos = true/false
 *
 *		IF params.filter_complete_chromos == true :
 *		ML phylogenies will be inferred for both the raw and filtered alignment
 *		Otherwise only on the raw alignments
 */


/* 
 * 
 *		ML PHYLOGENY BASED ON CONCATENATION OF WINDOW SUBSETS FROM ALL AUTOSOMES
 *
 */

process infer_auto_concat_subset_phy {

	/*
	 * Take the concat window MSA for each subset and infer phylogeny using IQtree2
	 *
	 * Directed into one channel:
	 * - ch1: process (windowCF_sCF_auto_subset_CONCAT_ST) -- sCF/wCF annotation
	 */
    
	label 'IC_HM_HT'

	tag "Infer phylogeny of concatenated autosomal windows for each subset length"
	publishDir "${params.outputdir}/02.phylogenies/subset/01.autosomes/", mode:'copy', overwrite:'false'

	input:
	tuple val(subset), path(msa_fn) from subset_window_concat_autosome_msa_ch1

	output:
	tuple val(subset), file("${subset}bp_autosomal_windows_concat.treefile") into phy_auto_concat_subset_ch1

	when:
	params.infer_ML_on_subsets == true 

	script:
	"""

	iqtree \
	-s ${msa_fn}  \
	-m ${params.phy_model_concat} \
	-B 1000 \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix ${subset}bp_autosomal_windows_concat

	"""
}


/* 
 * 
 *		ML PHYLOGENY BASED ON CONCATENATION OF COMPLETE AUTOSOMES
 *
 *		Raw (and filtered?)
 *
 */


process infer_auto_concat_raw_phy {

	/*
	 * ML phylogeny for the concatenated autosomes; RAW (always calculated as long as infer_ML_on_subsets == false )
	 *
	 * Directed into two channels
	 * 1) process (windowCF_sCF_auto_complete_raw_CONCAT). Calculate wCF/sCF on ML phylogeny
	 *
	 * NOTE: This can take a lot of RAM and take a while if the genome is very large
	 */
    
	label 'IC_HM_HT'

	tag "Infer phylogeny for concatenated autosomes -- raw"
	publishDir "${params.outputdir}/02.phylogenies/complete/01.autosomes/00.raw/", mode:'copy', overwrite:'false'

	input:
	path(msa) from concat_msa_auto_raw_ch1

	output:
	file("autosomes_concat_raw.treefile") into phy_auto_concat_raw_ch1

	when:
	params.infer_ML_on_subsets == false

	script:
	"""

	iqtree \
	-s ${msa}  \
	-m ${params.phy_model_concat} \
	-B 1000 \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix autosomes_concat_raw

	"""
}


process infer_auto_concat_filt_phy {

	/*
	 * ML phylogeny for the concatenated autosomes; filt
	 *
	 * Directed into two channels
	 * 1) process (windowCF_sCF_auto_complete_filt_CONCAT). Calculate wCF/sCF on ML phylogeny
	 *
	 * NOTE: This can take a lot of RAM and take a while if the genome is very large
	 */
    
	label 'IC_HM_HT'

	tag "Infer phylogeny for concatenated autosomes -- filt"
	publishDir "${params.outputdir}/02.phylogenies/complete/01.autosomes/01.filtered/", mode:'copy', overwrite:'false'

	input:
	path(msa) from concat_msa_auto_filt_ch1

	output:
	file("autosomes_concat_filt.treefile") into phy_auto_concat_filt_ch1

	when:
	params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	"""

	iqtree \
	-s ${msa}  \
	-m ${params.phy_model_concat} \
	-B 1000 \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix autosomes_concat_filt

	"""
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		MAXIMUM-LIKELIHOOD PHYLOGENY FOR EACH CHROMOSOME 	 			   *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		ML phylogeny for each chromosome
 *
 *		Only done if params.infer_chromo_phy = true
 *		
 *		This can either be for complete chromosomes or on a concatenated
 *		dataset that includes all subset windows for a given chromosome
 * 		
 * 		params.infer_ML_on_subsets = true/false
 *
 *		IF params.infer_ML_on_subsets == false :
 *		The ML phylogeny will be based on complete chromosomes.
 *		Is there also need for a ML phylogeny based on filtered alignment?
 *		
 *		params.filter_complete_chromos = true/false
 *
 *		IF params.filter_complete_chromos == true :
 *		ML phylogenies will be inferred for both the raw and filtered alignment
 */



/* 
 * 
 *		ML PHYLOGENY FOR EACH CHROMOSOME, BASED ON CONCATENATION OF SUBSETS FOR THAT CHROMOSOME
 *
 *
 */


process infer_chromo_concat_subset_phy {

	/*
	 * ML phylogeny for each chromosome based on a concatenated subset of the chromosome
	 *
	 * Directed into one channel
	 * 1) process (windowCF_sCF_sexchr_concat_subsets) -- CF annotation on sex chromosomes
	 *
	 */

	label 'IC_IM_HT'
    
	tag "Infer phylogeny for each chromosome, based on concatenated subsets of each chromo"
	publishDir "${params.outputdir}/02.phylogenies/subset/00.chromosomes/${subset}", mode:'copy'

	input:
	tuple val(chromo), val(subset), path(chromo_msa) from subset_window_concat_chromo_msa_ch1

	output:
	tuple val(chromo), val(subset), path("${chromo}_filt.treefile") into phy_by_chromo_concat_subset_ch1

	when:
	params.infer_chromo_phy == true && params.infer_ML_on_subsets == true 

	script:
	"""

	iqtree \
	-s ${chromo_msa}  \
	-m ${params.phy_model_concat} \
	-B 1000 \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix ${chromo}_filt

	"""
}


/* 
 * 
 *		ML PHYLOGENY FOR COMPLETE CHROMOSOMES
 *
 *		Raw (and filtered?)
 *
 */


process infer_chromo_raw_phy {

	/*
	 * ML phylogeny for each chromosome; RAW
	 *
	 * Directed into two channels
	 * 1) process (wCF_sCF_sex_chromo) --> CF annotation on sex chromosomes (RAW)
	 *
	 * NOTE: This can take a lot of RAM if the chromos are very large
	 */

	label 'IC_HM_HT'
    
	tag "Infer phylogeny for all chromosomes -- RAW"
	publishDir "${params.outputdir}/02.phylogenies/complete/00.chromosomes/00.raw/", mode:'copy'

	input:
	tuple val(chromo), path(chromo_msa) from msa_by_chromo_raw_ch4

	output:
	tuple val(chromo), path("${chromo}_raw.treefile") into phy_by_chromo_raw_ch1

	when:
	params.infer_chromo_phy == true && params.infer_ML_on_subsets == false

	script:
	"""

	iqtree \
	-s ${chromo_msa}  \
	-m ${params.phy_model_concat} \
	-B 1000 \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix ${chromo}_raw

	"""
}


process infer_chromo_filt_phy {

	/*
	 * ML phylogeny for each chromosome; FILTERED
	 *
	 * Directed into one channel
	 * 1) process (wCF_sCF_sex_chromo) --> CF annotation on sex chromosomes (FILTERED)
	 *
	 * NOTE: This can take a lot of RAM if the chromos are very large
	 */
    
	label 'IC_HM_HT'

	tag "Infer phylogeny for all chromosomes -- filt"
	publishDir "${params.outputdir}/02.phylogenies/complete/00.chromosomes/01.filtered/", mode:'copy'

	input:
	tuple val(chromo), path(chromo_msa) from msa_by_chromo_filt_ch2

	output:
	tuple val(chromo), path("${chromo}_filt.treefile") into phy_by_chromo_filt_ch1

	when:
	params.infer_chromo_phy == true && params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	"""

	iqtree \
	-s ${chromo_msa}  \
	-m ${params.phy_model_concat} \
	-B 1000 \
	-nt AUTO \
	-ntmax ${task.cpus} \
	--redo \
	--prefix ${chromo}_filt

	"""
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		SUMMARY-COALESCENT PHYLOGENIES	                   *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		Summary-coalescent phylogenies based on all autosomal
 *		window phylogenies for each subset size.
 *		
 *		For each subset, first all phylogenies need to be 
 * 		gathered in a single tree file.
 * 
 */


/* 
 * 		Exclude sex chromosome windows phylogenies 
 */

window_phy_by_subset_ch1
	.branch { 
		auto : it[0] !in params.sex_chromos 
		sex : it[0] in params.sex_chromos 
	}
	.set { window_phy_by_subset_chromo_type_ch1 }


window_phy_by_subset_chromo_type_ch1.auto
	.groupTuple(by: 1)
	.map{ it ->

		def chromo = it[0]
		def subset = it[1]
		def phy_list = it[2].flatten()

		[subset, phy_list]
	}
	.set { auto_window_phy_by_subset_ch1 }


// Check if sex chromos listed
// NOTE: There can be more than 1 sex chromo and each of them we will treat independently. Therefore GroupBy both chromo and subset
if (params.sex_chromos != []) {
	window_phy_by_subset_chromo_type_ch1.sex
		.groupTuple(by: [0,1])
		.map{ it ->

			def chromo = it[0]
			def subset = it[1]
			def phy_list = it[2].flatten()

			[chromo, subset, phy_list]
		}
		.set { sex_window_phy_by_subset_ch1 }
} else {
	sex_window_phy_by_subset_ch1 = Channel.empty()
}


process gather_subset_autosomal_window_trees {

	/*
	 * Gather all autosomal trees, per subset size, into a single tree file
	 *
	 * Directed into four channels:
	 * - ch1: Process (infer_subset_ASTRAL_ST_phy) -- Species tree inference in ASTRAL for each subset size
	 * - ch2: Process (windowCF_sCF_auto_subset_CONCAT) -- CF annotation on ML phylogeny based on all concatenated autosome windows - ALL autosomes SUBSET
	 * - ch3: Process (windowCF_sCF_auto_concat_raw) -- CF annotation on ML phylogeny based on complete concatenated autosomes - All autosomes RAW
	 * - ch4: Process (windowCF_sCF_auto_concat_filt) -- CF annotation on ML phylogeny based on complete concatenated autosomes - All autosomes FILT
	 * - ch5: Process (windowCF_sCF_auto_subset_SUMCOAL) -- CF annotation on summary-coalescent phylogeny
	 */
    
	tag "Gather autosomal window trees -- subset "
	publishDir "${params.outputdir}/02.phylogenies/subset/02.windows/${subset}/", mode:'copy'

	input:
	tuple val(subset), path(auto_window_trees_list) from auto_window_phy_by_subset_ch1

	output:
	tuple val(subset), path("${subset}_bp_autosomal_windows.trees") into subset_window_trees_combined_ch1, \
	subset_window_trees_combined_ch2, \
	subset_window_trees_combined_ch3, \
	subset_window_trees_combined_ch4, \
	subset_window_trees_combined_ch5

	script:
	"""

	cat ${auto_window_trees_list} > ${subset}_bp_autosomal_windows.trees

	"""
}


process gather_subset_sex_window_trees {

	/*
	 * Gather all sex chromosome trees, per subset size, into a single tree file (NEEDED FOR DOWNSTREAM wCF/cCF annotation of sex chromos)
	 *
	 * Directed into two channels:
	 * - ch1: Process (windowCF_sCF_sexchr_complete_filt) -- CF annotation on filtered complete sex chromosomes
	 * - ch2: Process (windowCF_sCF_sexchr_complete_raw) -- CF annotation on raw complete sex chromosomes
	 * - ch3: Process (windowCF_sCF_sexchr_concat_subsets) -- CF annotation on concatenated subsets for each sex chromosome
	 */
    
	tag "Gather sex window trees -- subset "
	publishDir "${params.outputdir}/02.phylogenies/subset/02.windows/${subset}/", mode:'copy'

	input:
	tuple val(chromo), val(subset), path(sex_window_trees_list) from sex_window_phy_by_subset_ch1

	output:
	tuple val(chromo), val(subset), path("${chromo}_${subset}_bp_sex_chromo_windows.trees") into subset_sex_window_trees_combined_ch1, \
	subset_sex_window_trees_combined_ch2, \
	subset_sex_window_trees_combined_ch3

	when:
	params.sex_chromos != []

	script:
	"""

	cat ${sex_window_trees_list} > ${chromo}_${subset}_bp_sex_chromo_windows.trees

	"""
}


process infer_subset_ASTRAL_ST_phy {

	/*
	 * Infer species tree using ASTRAL3, for each window subset size, using the inferred window trees
	 *
	 * Directed into one channel:
	 * - ch1: process (windowCF_sCF_auto_subset_SUMCOAL) -- CF annotation on sum coal phylogeny
	 */

	tag "Create window phylogenies for each chromo and subset window length"
	publishDir "${params.outputdir}/02.phylogenies/subset/03.summary_coalescent/${subset}", mode:'copy'

	input:
	tuple val(subset), file(window_trees_fn) from subset_window_trees_combined_ch1

	output:
	tuple val(subset), file("*_ASTRAL.tree") into subset_window_astral_phy_ch1

	script:
	def prefix = window_trees_fn.getBaseName()
	"""

	java -jar ${params.astral_fn} \
	-i ${window_trees_fn} \
	-o ${prefix}_ASTRAL.tree

	"""
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * 		TIME FOR TREE ANNOTATION!! 												   *
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 																				   *
 * 		- sCF, wCF on ML phylogeny for all autosomes combined
 * 		- sCF, wCF on ML phylogeny of sex chromosome
 * 		- sCF, wCF on Summary-coalescent phylogenies						       *
 *
 * 																				   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		sCF, wCF on ML phylogeny for all autosomes	       *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *		WHICH ML phylogeny depends on:
 *
 *		params.infer_ML_on_subsets = true/false
 *
 *		params.filter_complete_chromos = true/false
 *
 */


/* 
 *
 *		ML PHYLOGENY BASED ON CONCATENATION OF WINDOW SUBSETS FROM ALL AUTOSOMES
 *		params.infer_ML_on_subsets = true
 *
 * 		To calculate sCF and wCF the following is needed:
 *		
 *		- Alignment
 *		- Tree files
 * 		- ML phylogeny 
 *
 *		Therefore create a new channel that combines all three:
 *		[ subset, subset_concat_aln.fa, subset_window_trees.tree, subset_concat_window.tree ]
 */


if (params.infer_ML_on_subsets == true) {
	subset_window_concat_autosome_msa_ch2 
		.join(subset_window_trees_combined_ch2)
		.join(phy_auto_concat_subset_ch1)
		.set{ subset_auto_window_cf_input_ch1 }
} else {
	subset_auto_window_cf_input_ch1 = Channel.empty()
}


process windowCF_sCF_auto_subset_CONCAT {

	/*
	 * Calculate and annotate window and site concordance factors on: CONCAT species tree (per window subset size)
	 */
    
	label 'SC_IM_IT'

	tag "Calculate windowCF and sCF, per window subset, on: Concatenated autosome tree"
	publishDir "${params.outputdir}/02.phylogenies/subset/01.autosomes", mode:'copy', overwrite:'false'

	input:
	tuple val(subset), \
	file(subset_concat_window_aln_fn), \
	file(subset_window_trees_fn), \
	file(subset_concat_window_tree_fn) from subset_auto_window_cf_input_ch1

	output:
	file("*_wCF_sCF.cf.tree")
	file("*_wCF_sCF.cf.tree.nex")
	file("*_wCF_sCF.cf.branch")
	file("*_wCF_sCF.cf.stat")

	when:
	params.infer_ML_on_subsets == true

	script:
	def prefix = "${subset}bp_autosomal_windows_concat"
	"""

	iqtree \
	-t ${subset_concat_window_tree_fn} \
	-s ${subset_concat_window_aln_fn} \
	--gcf ${subset_window_trees_fn} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}_wCF_sCF

	"""
}


/* 
 *
 *		ML PHYLOGENY BASED ON CONCATENATION OF COMPLETE CHROMOSOMES
 *		params.infer_ML_on_subsets = true
 *
 *		1. unfiltered
 *		2. filtered
 *
 * 		To calculate sCF and wCF the following is needed:
 *		
 *		- Alignment 
 *		- Tree files (by subset)
 * 		- ML phylogeny 
 *
 *		Therefore create a new channel that combines all three:
 *		[ autosome_concat_aln.fa, subset, subset_window_trees.tree, autosome_concat.tree ]
 *
 *		NOTE: The same alignment and ML phylogeny file are used.
 * 			  However, different wCF are annotated, based on the different subset sizes
 */


if (params.infer_ML_on_subsets == false) {
	concat_msa_auto_raw_ch2
		.combine(subset_window_trees_combined_ch3)
		.combine(phy_auto_concat_raw_ch1)
		.set { wCF_sCF_concat_raw_input_ch1 }
} else {
	wCF_sCF_concat_raw_input_ch1 = Channel.empty()
}

process windowCF_sCF_auto_complete_raw_CONCAT {

	/*
	 * Calculate and annotate window and site concordance factors on: Concatenated autosome tree -- RAW
	 */
    
	label 'SC_IM_IT'

	tag "Calculate windowCF and sCF, per window subset, on: Concatenated autosome tree"
	publishDir "${params.outputdir}/02.phylogenies/complete/01.autosomes/00.raw/", mode:'copy', overwrite:'false'

	input:
	tuple file(autosome_concat_aln_fn), \
	val(subset), \
	file(subset_window_trees_fn), \
	file(autosome_concat_tree_fn) from wCF_sCF_concat_raw_input_ch1

	output:
	file("*_raw.cf.tree")
	file("*_raw.cf.tree.nex")
	file("*_raw.cf.branch")
	file("*_raw.cf.stat")

	when:
	params.infer_ML_on_subsets == false

	script:
	def prefix = "${subset}bp_wCF_complete_sCF_on_autosomes_concat_raw"
	"""

	iqtree \
	-t ${autosome_concat_tree_fn} \
	-s ${autosome_concat_aln_fn} \
	--gcf ${subset_window_trees_fn} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}

	"""
}


if (params.infer_ML_on_subsets == false && params.filter_complete_chromos == true) {
	concat_msa_auto_filt_ch2
		.combine(subset_window_trees_combined_ch4)
		.combine(phy_auto_concat_filt_ch1)
		.set { wCF_sCF_concat_filt_input_ch1 }
} else {
	wCF_sCF_concat_filt_input_ch1 = Channel.empty()
}


process windowCF_sCF_auto_complete_filt_CONCAT {

	/*
	 * Calculate and annotate window and site concordance factors on: Concatenated autosome tree -- FILT
	 */
    
	label 'SC_IM_IT'

	tag "Calculate windowCF and sCF, per window subset, on: Concatenated autosome tree"
	publishDir "${params.outputdir}/02.phylogenies/complete/01.autosomes/01.filt/", mode:'copy', overwrite:'false'

	input:
	tuple file(autosome_concat_aln_fn), \
	val(subset), \
	file(subset_window_trees_fn), \
	file(autosome_concat_tree_fn) from wCF_sCF_concat_filt_input_ch1

	output:
	file("*_filt.cf.tree")
	file("*_filt.cf.tree.nex")
	file("*_filt.cf.branch")
	file("*_filt.cf.stat")

	when:
	params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	def prefix = "${subset}bp_wCF_complete_sCF_on_autosomes_concat_filt"
	"""

	iqtree \
	-t ${autosome_concat_tree_fn} \
	-s ${autosome_concat_aln_fn} \
	--gcf ${subset_window_trees_fn} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}

	"""
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		sCF, wCF on ML phylogenies of sex-chromosome(s)	       *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */


/* 
 *
 *		ML PHYLOGENY FOR SEX CHROMOSOME
 *
 * 		To calculate sCF and wCF the following is needed:
 *		
 *		- Alignment (sex chromo only)
 *		- Tree files (sex chromo only, by subset sizes)
 * 		- ML phylogeny (sex chromo, by subset sizes)
 *
 */

if (params.sex_chromos != [] && params.infer_chromo_phy == true && params.infer_ML_on_subsets == false) {
	// THERE ARE SEX CHROMOS AND PHYLOGENIES HAVE BEEN INFERRED ON COMPLETE CHROMOSOMES
	if (params.filter_complete_chromos == true) {
		/*
		*	Inferred on complete chromosomes, filtered and raw
		* 
		*	Desired structure of output channel:
 		*	[ chromo, chromo_aln.fa, chromo_phy.tree, subset, subset_chromo_trees.tree ]
 		*/

		// First select the sex chromosome phylogenies
		phy_by_chromo_filt_ch1
			.branch { 
				auto : it[0] !in params.sex_chromos 
				sex : it[0] in params.sex_chromos 
			}
			.set { phy_by_filt_chromo_type_ch1 }

		// Join the alignment and phylogenies by chromo
		msa_by_chromo_filt_chromo_type_ch1.sex
			.join(phy_by_filt_chromo_type_ch1.sex)
			.join(subset_sex_window_trees_combined_ch1)
			.set { wCF_sCF_sexchr_filt_input_ch1 }

		/* 
		*	Repeat for the unfiltered alignments and phylogenies
 		*/

		phy_by_chromo_raw_ch1
			.branch { 
				auto : it[0] !in params.sex_chromos 
				sex : it[0] in params.sex_chromos 
			}
			.set { phy_by_raw_chromo_type_ch1 }

		msa_by_chromo_raw_chromo_type_ch1.sex
			.join(phy_by_raw_chromo_type_ch1.sex)
			.join(subset_sex_window_trees_combined_ch2)
			.set { wCF_sCF_sexchr_raw_input_ch1 }

	} else {

		/*
		*	Inferred on complete chromosomes, RAW only
		* 
		*	Desired structure of output channel:
 		*	[ chromo, chromo_aln.fa, chromo_phy.tree, subset, chromo_trees.tree ]
 		*/

		phy_by_chromo_raw_ch1
			.branch { 
				auto : it[0] !in params.sex_chromos 
				sex : it[0] in params.sex_chromos 
			}
			.set { phy_by_raw_chromo_type_ch1 }

		msa_by_chromo_raw_chromo_type_ch1.sex
			.join(phy_by_raw_chromo_type_ch1.sex)
			.join(subset_sex_window_trees_combined_ch2)
			.set { wCF_sCF_sexchr_raw_input_ch1 }

		wCF_sCF_sexchr_filt_input_ch1 = Channel.empty()
		windowCF_sCF_sex_chr_concat_subsets_input_ch1 = Channel.empty()
	}

} else if (params.sex_chromos != [] && params.infer_chromo_phy == true && params.infer_ML_on_subsets == true) {
	// THERE ARE SEX CHROMOS AND PHYLOGENIES, BUT ONLY HAVE BEEN INFERRED ON CONCATENATED SUBSETS OF CHROMOSOMES

	/*
	*	Inferred on complete chromosomes, ALWAYS filtered. However, there will now be MULTIPLE trees for each sexchromo
	*	Each being of different subset size. These therefore need to be combined correctly. (NOTE difference between combine and join)
	* 
	*	Desired structure of output channel:
	*	[ chromo, subset, chromo_aln.fa, chromo_phy.tree, chromo_trees.tree ]
	*/

	subset_window_concat_chromo_msa_ch2
		.branch { 
			auto : it[0] !in params.sex_chromos 
			sex : it[0] in params.sex_chromos 
		}
		.set { msa_by_subset_window_concat_chromo_type_ch1 }

	phy_by_chromo_concat_subset_ch1
		.branch { 
			auto : it[0] !in params.sex_chromos 
			sex : it[0] in params.sex_chromos 
		}
		.set { phy_by_subset_window_concat_chromo_type_ch1 }

	msa_by_subset_window_concat_chromo_type_ch1.sex
		.combine(phy_by_subset_window_concat_chromo_type_ch1.sex, by:[0,1])
		.combine(subset_sex_window_trees_combined_ch3, by:[0,1])
		.set { windowCF_sCF_sex_chr_concat_subsets_input_ch1 }
	
	wCF_sCF_sexchr_filt_input_ch1 = Channel.empty()
	wCF_sCF_sexchr_raw_input_ch1 = Channel.empty()

} else if (params.sex_chromos == [] || params.infer_chromo_phy == false ) {

	windowCF_sCF_sex_chr_concat_subsets_input_ch1 = Channel.empty()
	wCF_sCF_sexchr_filt_input_ch1 = Channel.empty()
	wCF_sCF_sexchr_raw_input_ch1 = Channel.empty()
	
	}



/*
*	With input channels prepared, now annotate the sex chromosomes
*/

process windowCF_sCF_sexchr_complete_raw {

	/*
	 * Calculate and annotate window and site concordance factors on: Each sex chromosome, complete, RAW
	 *
	 * Input channels needed:
	 * wCF_sCF_sexchr_raw_input_ch1 	-- 		[ chromo, chromo_aln.fa, chromo_phy.tree, subset, subset_chromo_trees.tree ]
	 *
	 */
    
	label 'SC_IM_IT'
   
	tag "Calculate windowCF and sCF, per window subset, on: Sex chromo phylo - raw complete"
	publishDir "${params.outputdir}/02.phylogenies/complete/00.chromosomes/01.raw", mode:'copy', overwrite:'false'

	input:
	tuple val(chromo), \
	file(sex_chromo_msa), \
	file(sex_chromo_phy), \
	val(subset), \
	file(subset_sex_chromo_trees) from wCF_sCF_sexchr_raw_input_ch1

	output:
	file("*_wCF_sCF.cf.tree")
	file("*_wCF_sCF.cf.tree.nex")
	file("*_wCF_sCF.cf.branch")
	file("*_wCF_sCF.cf.stat")

	when:
	params.sex_chromos != [] && params.infer_chromo_phy == true && params.infer_ML_on_subsets == false

	script:
	def prefix = "${chromo}_${subset}bp_concat_raw"
	"""

	iqtree \
	-t ${sex_chromo_phy} \
	-s ${sex_chromo_msa} \
	--gcf ${subset_sex_chromo_trees} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}_wCF_sCF

	"""
}


process windowCF_sCF_sexchr_complete_filt {

	/*
	 * Calculate and annotate window and site concordance factors on: Each sex chromosome, complete, FILTERED
	 *
	 * Input channels needed:
	 * wCF_sCF_sexchr_filt_input_ch1 	-- 		[ chromo, chromo_aln.fa, chromo_phy.tree, subset, subset_chromo_trees.tree ]
	 *
	 */
    
	label 'SC_IM_IT'
   
	tag "Calculate windowCF and sCF, per window subset, on: Sex chromo phylo - filtered complete"
	publishDir "${params.outputdir}/02.phylogenies/complete/00.chromosomes/01.filtered", mode:'copy', overwrite:'false'

	input:
	tuple val(chromo), \
	file(sex_chromo_msa), \
	file(sex_chromo_phy), \
	val(subset), \
	file(subset_sex_chromo_trees) from wCF_sCF_sexchr_filt_input_ch1

	output:
	file("*_wCF_sCF.cf.tree")
	file("*_wCF_sCF.cf.tree.nex")
	file("*_wCF_sCF.cf.branch")
	file("*_wCF_sCF.cf.stat")

	when:
	params.sex_chromos != [] && params.infer_chromo_phy == true && params.infer_ML_on_subsets == false && params.filter_complete_chromos == true

	script:
	def prefix = "${chromo}_${subset}bp_concat_filt"
	"""

	iqtree \
	-t ${sex_chromo_phy} \
	-s ${sex_chromo_msa} \
	--gcf ${subset_sex_chromo_trees} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}_wCF_sCF

	"""
}


process windowCF_sCF_sexchr_concat_subsets {

	/*
	 * Calculate and annotate window and site concordance factors on: Each sex chromosome phylogeny (concatenated subsets)
	 *
	 * Input channels needed:
	 * windowCF_sCF_sex_chr_concat_subsets_input_ch1 	-- 	[ chromo, subset, chromo_aln.fa, chromo_phy.tree, chromo_trees.tree ]
	 *
	 */
    
	label 'SC_IM_IT'
   
	tag "Calculate windowCF and sCF, per window subset, on: Sex chromo phylo - concatenated subsets"
	publishDir "${params.outputdir}/02.phylogenies/subset/00.chromosomes/${subset}", mode:'copy', overwrite:'false'

	input:
	tuple val(chromo), \
	val(subset), \
	file(sex_chromo_msa), \
	file(sex_chromo_phy), \
	file(subset_sex_chromo_trees) from windowCF_sCF_sex_chr_concat_subsets_input_ch1

	output:
	file("*_wCF_sCF.cf.tree")
	file("*_wCF_sCF.cf.tree.nex")
	file("*_wCF_sCF.cf.branch")
	file("*_wCF_sCF.cf.stat")

	when:
	params.sex_chromos != [] && params.infer_chromo_phy == true && params.infer_ML_on_subsets == true

	script:
	def prefix = "${chromo}_filt"
	"""

	iqtree \
	-t ${sex_chromo_phy} \
	-s ${sex_chromo_msa} \
	--gcf ${subset_sex_chromo_trees} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}_wCF_sCF

	"""
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * 		sCF, wCF on Summary-coalescent phylogenies	       *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *		Summary-coalescent phylogeny has been inferred
 *		using different sized subsets.
 *
 * 		To calculate chromoCF and wCF the following is needed:
 *		
 *		- Alignment (Concatenation of all subset window alignments)
 *		- Tree files (All window trees combined, by subset size)
 * 		- Summary-coalescent phylogeny (by subset sizes)
 *
 */


/*
*	Create a single channel for wCF/sCF
* 
*	Desired structure of output channel:
*	[ subset, subset_concat_window_aln.fa, subset_window_trees.tree, subset_summary_coalescent_phylogeny.tree ]
*/

subset_window_concat_autosome_msa_ch3
	.join(subset_window_trees_combined_ch5)
	.join(subset_window_astral_phy_ch1)
	.set { windowCF_sCF_auto_subset_SUMCOAL_input_ch1 }


process windowCF_sCF_auto_subset_SUMCOAL {

	/*
	 * Calculate and annotate window and site concordance factors on: Summary coalescent species tree (per window subset size)
	 *
	 * Input channels needed:
	 * windowCF_sCF_auto_subset_SUMCOAL_input_ch1 	-- 	[ subset, subset_concat_window_aln.fa, subset_window_trees.tree, subset_summary_coalescent_phylogeny.tree ]
	 */
    
	label 'SC_IM_IT'

	tag "Calculate windowCF and sCF on: Summary-Coalescent tree -- per window subset length"
	publishDir "${params.outputdir}/02.phylogenies/subset/03.summary_coalescent/${subset}/", mode:'copy', overwrite:'false'

	input:
	tuple val(subset), \
	file(subset_concat_window_aln_fn), \
	file(subset_window_trees_fn), \
	file(sumcoal_tree_fn) from windowCF_sCF_auto_subset_SUMCOAL_input_ch1

	output:
	file("*_wCF_sCF.cf.tree")
	file("*_wCF_sCF.cf.tree.nex")
	file("*_wCF_sCF.cf.branch")
	file("*_wCF_sCF.cf.stat")

	script:
	def prefix = sumcoal_tree_fn.getBaseName()
	"""

	iqtree \
	-t ${sumcoal_tree_fn} \
	-s ${subset_concat_window_aln_fn} \
	--gcf ${subset_window_trees_fn} \
	--scf 100 \
	-T ${task.cpus} \
	--prefix ${prefix}_wCF_sCF

	"""
}

