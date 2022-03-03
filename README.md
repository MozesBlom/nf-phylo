**A reproducible pipeline for phylogenomic inference using whole-genome data**.
**v1.0**
[TOC]

## Introduction

**nf-phylo** is a bioinformatics pipeline that takes a reference genome, variant calls and mask files, and conducts a thorough phylogenetic analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io) (DSL1), a workflow tool to run tasks across multiple computing infrastructures in a very portable manner. It currently creates a custom Conda environment from scratch and therefore does not require any pre-compiled software packages (besides Nextflow itself). Future development may include containerized versions as well to further enhance reproducibility. If possible, e.g. when running on a HPC cluster, the pipeline will process alignments, infer trees, etc. in parallel and all batch submission jobs are handled internally through the Nextflow workflow.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04) 
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10)
3. Download the pipeline and adjust the config file. Here you need to specify the location of the input data, filtering settings etc.

    ```bash
    ./nf-phylo/nextflow.config
    ```
4. Create a nextflow config profile that matches your cluster set-up ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) and start running your own analysis!

    ```bash
    nextflow run ./nf-phylo/phylo.nf -profile mfn
    ```

6. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```

If a run stalls for a given reason, inspect the error message, adjust and you can rerun using the 'resume' flag.

    ```bash
    nextflow run ./nf-phylo/phylo.nf -profile mfn -resume
    ```


**N.B.** Please make sure that:
* The scripts in `./bin/xxx.py` can be executed by any user. You can check the user permissions via:

    ```bash
    ls -l ./bin/
    ```
And if need be changed via:

    ```bash
    cd ./bin/
    chmod 755 *
    ```

* There is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored). On the MfN cluster, running on `/home/` will easily lead to problems and all pipelines will need to be executed from a `/data4/` project folder.

* The ASTRAL executable is [`downloaded`](https://github.com/smirarab/ASTRAL) and can be succesfully started. In example:
    ```bash
    java -jar astral.5.7.8.jar -h
    ```

## Pipeline Summary

Depending on genome size, computational resources and time available, NF-phylo can be adjusted to find the most optimal balance between inference needed and resources available. To be more precise, phylogenetic analyses can be based on a subset of the genome or on a complete datamatrix of the entire genome. Window based analyses are the default and concatenated analyses can be based on entire chromosomes or on concatenated window subsets that are sampled at a user specified density (eg. 10kb ever 100 kb.). Complete analyses can be repeated with different window lengths to evaluate the importance of window size. The workflow follows the following structure:

* If need be, call consensus sequences for each chromosome/scaffold and each individual (`samtools`, `bcftools`)
* Estimate the amount of missing data for each individual (`bin/00[a/b]_consensus_stats.py`)
* Merge individual consensus sequences into multiple sequence alignments per chromosome/scaffold
* Subset each MSA into windows of size x (length x specified in config) and at distance y (density of window sampling specified in config) (`bin/01_window_subset_MSA.py`)
* Infer a phylogeny for each window using maximum-likelihood (`IQtree2`). Windows are filtered based on user-specified settings.
* Create a concatenated alignment of all autosomes (One of two options below)
    * a) Subset: Concatenation of window alignments of size x (length x specified in config) and at distance y (density of window sampling specified in config) (`Seqkit`)
    * b) Full chromosomes: Concatenation of complete chromosomes.
        * b1) unfiltered 
        * b2) filtered  (OPTIONAL; `bin/02_chromo_filt_msa.py`)
            * b1.1) filter entire chromosomes in one go
            * b1.2) filter chromosome in subsets that are subsequently stitched back together (`bin/03_concat_msa.py`)
* Create an alignment for each chromosome (One of two options below)
    * a) Subset: Concatenation of window alignments of size x (length x specified in config) and at distance y (density of window sampling specified in config) (`Seqkit`)
    * b) Complete chromosomes.
        * b1) unfiltered 
        * b2) filtered  (OPTIONAL; `bin/02_chromo_filt_msa.py`)
            * b1.1) filter entire chromosome in one go
            * b1.2) filter chromosome in subsets that are subsequently stitched back together (`bin/03_concat_msa.py`)
* Infer a phylogeny for a concatenated alignment of all autosomes using maximum-likelihood (`IQtree2`). (One of two options below)
    * a) Subset
    * b) Complete chromosomes.
* Infer a phylogeny for each chromosome using maximum-likelihood (`IQtree2`). [OPTIONAL]
    * a) Subset
    * b) Complete chromosomes
        * b1) unfiltered 
        * b2) filtered [OPTIONAL]
* Infer a summary-coalescent phylogeny based on all autosomal windows (`ASTRAL3`). If windows of different lengths have been generated a SUMCOAL species tree will be inferred for each window length.
* Concordance factor annotation on all phylogenies inferred. (`IQtree2`)


### Input
The overarching aim of the pipeline is to generate species trees/phylogenies and quantify the phylogenetic uncertainty for each inferred bipartition by taking into account (ideally) genome scale data. The pipeline can start with consensus calling (each individual, each scaffold/chromosome) or with consensus sequences themselves. If starting with consensus calling, the pipeline requires a `variant call set` and a `mask set` for each chromosome and stored in a single folder per individual:

variant call set = `/inputdir/[indiv]/[chromo]_vars_filt_indels.vcf.gz`
mask set = `/inputdir/[indiv]/[chromo]_mask_cov_het.vcf`

Alternatively, it is also possible to start with consensus sequences. A single consensus for each individual and chromosome which need to be specified in a tab seperated input file with format (and headers):

| Individual  | Scaffold |        Consensus_fn       |
| :---------: | :-------:| :-----------------------: |
| IndivA      | Chr1     |  /path/to/indivA_chr1.fa  |
| IndivA      | Chr2     |  /path/to/indivA_chr2.fa  |
| IndivB      | Chr1     |  /path/to/indivB_chr1.fa  |

Run specific options can  be specified in the `./nextflow.config` script. All other values for programs are set at default. Boolean values can be specified as `true` or `false` and path names need to be absolute. If preferred, the same options listed in the config file can also be directly modified by using a parameter flag when initiating the nextflow run. Example given:

    ```bash
    nextflow run nf-phylo/phylo.nf -profile mfn --indivs_file /path/to/indivs.txt --chromos_file /path/to/chromos.txt --outdir /path/to/results --consensus_calling true
    ```

**NOTE**:
- NF-phylo only does a 'pseudo-alignment', meaning that it assumes that the consensus sequences for each chromosome/scaffold is the same length across all individuals. Indel variation should therefore have been removed from the variant call set or not incorporated into the existing sequences.
- Partially motivated for the required above, variant and mask sets need to be called on the SAME reference assembly across all individuals.
- Sex chromosomes should be treated differently due to their distint inheritance history. Such chromosomes should either be excluded from the (chromo) input list or assigned in the `sex_chromos` parameter. The difference is that such chromosomes are either ignored (the former option) or sex chromosome phylogenies are independently inferred (i.e. not included in the concatenation analyses)
- NF-phylo runs all inferences on the 'raw' consensus sequences and eventual alignments. However it also includes an option to filter the resulting alignments based on the relative ratio of missing data per data column (i.e. missing individuals per alignment site). Whereas this filtering is optional for the chromosome and concatenation alignments, it's done by default for the window based alignments. The window based alignments should be much smaller and the likelihood of alignments with a large proportion of missing data is therefore much higher. The risk of phylogenies based on large amounts of missing data is therefore much higher for the window based analyses.


### Output
The pipeline stores both alignments, phylogenies and log/summary statistics. Trees can be viewed with commonly used tree viewers such as [Figtree](http://tree.bio.ed.ac.uk/software/figtree/) and where appropriate support or concordance factors are annotated. Please see original software manuals for further details on output files.


