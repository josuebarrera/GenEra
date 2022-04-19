# Introduction
GenEra, which was developed at the Max Planck Institute for Biology Tübingen, is an easy-to-use, low-dependency and highly customizable command-line tool that estimates the age of the earliest common ancestor of protein-coding genes though genomic phylostratigraphy (Domazet-Lošo et al., 2007). GenEra takes advantage of DIAMOND’s speed and sensitivity to search for homolog genes throughout the entire NR database, and combines these results with the NCBI taxonomy to assign an origination date for each gene and gene family in a query species. GenEra can also incorporate protein data from external sources to enrich the analysis, it can search for proteins within nucleotide data (i.e., genome/transcriptome assemblies) using MMseqs2 to improve the classification of orphan genes, it automatically collapses phylostrata that lack genomic data for its inclusion in the analysis, and it calculates a taxonomic representativeness score to assess the reliability of assigning a gene to a specific phylostratum. Additionally, it can calculate the probability of a gene origination age to appear younger than it really is due to homology detection failure.

# Dependencies

GenEra requires the following software dependencies:

-	DIAMOND (https://github.com/bbuchfink/diamond)
-	NCBItax2lin (https://github.com/zyxue/ncbitax2lin)
-	MCL (https://github.com/micans/mcl)
-	MMseqs2 (https://github.com/soedinglab/MMseqs2) (optional)
-	abSENSE (https://github.com/caraweisman/abSENSE) (optional)
-	NumPy (https://numpy.org/) and SciPy (https://scipy.org/) (optional)

Additionally, GenEra requires a locally installed NR database for DIAMOND, as well as internet connection and access to the taxonomy dump from the NCBI.

# Installation

For an easy conda installation, copy and paste this in your terminal:

```console
git clone 
conda create -n genEra python=3.7
conda activate genEra
conda install -c bioconda diamond
conda install -c conda-forge -c bioconda mmseqs2
pip install -U ncbitax2lin
conda install -c bioconda mcl
conda install -c anaconda scipy
CONDABIN=$(which ncbitax2lin | sed 's/ncbitax2lin//g') && mv genEra ${CONDABIN} && mv Erassignation.sh ${CONDABIN}
```

Otherwise, you can install the dependencies independently and then include both genEra and Erassignation.sh to your PATH.
