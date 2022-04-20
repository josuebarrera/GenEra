![GenEra](https://github.com/josuebarrera/GenEra/blob/main/logo.png)

Introduction
============

GenEra is an easy-to-use, low-dependency and highly customizable command-line tool that estimates the age of the earliest common ancestor of protein-coding genes though genomic phylostratigraphy (Domazet-Lošo et al., 2007). GenEra takes advantage of [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND")’s speed and sensitivity to search for homolog genes throughout the entire NR database, and combines these results with the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy "NCBI Taxonomy") to assign an origination date for each gene and gene family in a query species. GenEra can also incorporate protein data from external sources to enrich the analysis, it can search for proteins within nucleotide data (i.e., genome/transcriptome assemblies) using [MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") to improve the classification of orphan genes, it automatically collapses phylostrata that lack genomic data for its inclusion in the analysis, and it calculates a taxonomic representativeness score to assess the reliability of assigning a gene to a specific phylostratum. Additionally, it can calculate the probability of a gene origination age to appear younger than it really is due to homology detection failure using [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE").

Contents
========

-   [Dependencies](#dependencies)
-   [Installation](#installation)
-   [Setting up the databases](#setting-up-the-databases)
-   [Quick start for the impatient](#quick-start-for-the-impatient)
-   [Arguments and input files](#arguments-and-input-files)
-   [Fine-tunning and other useful arguments](#fine-tunning-and-other-useful-arguments)
-   [Output files](#output-files)

Dependencies
============

GenEra requires the following software dependencies:

-	[DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND")
-	[NCBItax2lin](https://github.com/zyxue/ncbitax2lin "NCBItax2lin")
-	[MCL](https://github.com/micans/mcl "MCL")
-	[MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") (optional for protein-against-nucleotide sequence search)
-	[abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") (optional to calculate homology detection failure probabilities)
-	[NumPy](https://numpy.org/ "NumPy") and [SciPy](https://scipy.org/ "SciPy") (needed to run abSENSE in step 4)

Additionally, GenEra requires a locally installed NR database for DIAMOND, as well as internet connection and access to the taxonomy dump from the NCBI.

Installation
============

For an easy conda installation, copy and paste this in your terminal:

```console
git clone https://github.com/josuebarrera/GenEra.git
cd GenEra
conda create -n genEra python=3.7
conda activate genEra
conda install -c bioconda diamond
conda install -c bioconda mcl
pip install -U ncbitax2lin
conda install -c conda-forge -c bioconda mmseqs2
conda install -c anaconda scipy
CONDABIN=$(which ncbitax2lin | sed 's/ncbitax2lin//g') && mv genEra ${CONDABIN} && mv Erassignation.sh ${CONDABIN}
```

Otherwise, you can install the dependencies independently and then include both genEra and Erassignation.sh to your PATH.

Setting up the databases
========================

First, download the nr database (warning: this is a huge FASTA file):
```console
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr{.gz,.gz.md5} && md5sum -c *.md5
gunzip nr.gz
```
NOTE: Alternatively, you can download any other database whose sequence IDs can be traced back to the NCBI Taxonomy.

Then download “prot.accession2taxid” from the NCBI webpage:
```console
wget ftp://ftp.ncbi.nih.gov:21/pub/taxonomy/accession2taxid/prot.accession2taxid.gz && gunzip prot.accession2taxid.gz
```
Then download the taxonomy dump from the NCBI:
```console
wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir -p taxdump && tar zxf taxdump.tar.gz -C ./taxdump
```
Finally, create a local nr database (another huge file):
```console
diamond makedb \
 --in nr \
 --db nr \
 --taxonmap prot.accession2taxid \
 --taxonnodes taxdump/nodes.dmp \
 --taxonnames taxdump/names.dmp \
 --memory-limit 100
```
You can eliminate “prot.accession2taxid”, but keep the taxdump, as GenEra will use it later on.

Quick start for the impatient
=============================

The basic usage of GenEra is:
```console
genEra -q [query_sequences.fasta] -t [query_taxid] -b [path/to/nr] -d [path/to/taxdump]
```
However, several optional arguments can be included to get the best results possible:
```console
genEra -q [query_sequences.fasta] -t [query_taxid] -b [path/to/nr] -d [path/to/taxdump] \ 
-a [protein_list.tsv] -f [nucleotide_list.tsv] -s [evolutionary_distances.tsv] -i [true] -n [many threads]
  ```

Arguments and input files
=========================

### GenEra always requires these two inputs:

-  __-q__ &nbsp;&nbsp;&nbsp; A standard FASTA file containing all the protein sequences of your species of interest. Make sure all your sequence headers are unique.

-  __-t__ &nbsp;&nbsp;&nbsp; The NCBI taxonomy ID of your species of interest (it can be easily found at https://www.ncbi.nlm.nih.gov/taxonomy).

### Additionally, GenEra requires one of the following two files:

-  __-b__ &nbsp;&nbsp;&nbsp; A locally installed NR database for DIAMOND (or any other database whose sequences can be automatically traced back to NCBI taxids).

-  __-p__ &nbsp;&nbsp;&nbsp; A pre-generated BLAST/DIAMOND/MMseqs2 table. This file is automatically generated by GenEra in the first step of the pipeline, and can be given to GenEra in case the user decides to re-run the pipeline (e.g., fine-tunning of the parameters). The first column contains the query sequence IDs (qseqid), the second coulmn contains the tarjet sequence IDs (sseqid), the third column contains the e-value of the matching alignment (evalue), the fourth column contains the bitscore of the alignment (bitscore) and the fifth column contains the taxonomy ID of the tarjet sequence (staxids). 

### And one of the following three files:

-  __-d__ &nbsp;&nbsp;&nbsp; The location of the NCBI taxonomy files, colloquially known as the “taxdump”. Make sure the taxdump directory contains the files “nodes.dmp” and “names.dmp”. This files are used by NCBItax2lin to create a “ncbi_lineages” file.

-  __-r__ &nbsp;&nbsp;&nbsp; The location of the uncompressed “ncbi_lineages” file generated by NCBItax2lin. Once the user runs GenEra for the first time (or runs NCBItax2lin independently), this file can be used to save some time during step 2 of the pipeline. This file is a comma-delimited table with each row representing a specific lineage, and each column representing the taxonomic hierarchies to which that lineage belongs. This file will be automatically modified by genEra to rearrange the columns in hierarchical order.

-  __-c__ &nbsp;&nbsp;&nbsp; Custom "ncbi_lineages" file that is already tailored for the query species. GenEra modifies the raw “ncbi_lineages” file so that the taxid of the query species appears in the first column, and the taxonomic hierarchies of the query species are rearranged from the species level all the way back to the last universal common ancestor (termed “cellular organisms” by the NCBI taxonomy). The file is also modified by GenEra to collapse the phylostrata (i.e., the taxonomic levels) that lack the necessary genomic data to be useful in the analysis. Once this file is generated, the user can re-use this file if they want to run GenEra on the same species with different parameters. By default, GenEra will search for the correct hierarchical order of the query organism’s phylostrata directly from the NCBI webpage. If the user machine is unable to access the NCBI webpage, GenEra will attempt to infer the correct order of the phylostrata directly from the “ncbi_lineages” file. However, there are some taxonomic hierarchies that the NCBI labels as “clades” (e.g., Embryophyta in the plant lineage), which GenEra cannot automatically assign to their correct taxonomic level when working offline. This option is also useful if the user wants to modify the taxonomic hierarchies that were originally assigned by the NCBI (e.g., an outdated taxonomy that does not reflect the phylogenetic relationships of the query species). This custom table can be easily generated by printing the desired columns of the “ncbi_lineages” file with awk:
```console
awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' ncbi_lineages_[date].csv | awk -F "," '{ print $1","$8","$7","$6"," ... }' > custom_table.csv
```
### The user can also incorporate four other inputs, which are optional but can potentially improve the final age assignation:

-  __-a__ &nbsp;&nbsp;&nbsp; Tab-delimited table with additional proteins to be included in the analysis. This option is particularly useful if the user wants to include proteins from species that are absent from the nr database, such as newly annotated genomes or transcriptomes. It can also help fill the gaps of phylostrata that would otherwise be collapsed due to lack of available genomes in the database. Importantly, each protein in this custom dataset should have unique identifiers, otherwise GenEra will not work properly during the taxonomic assignation of the Diamond hits. The table format consists of the location of one protein FASTA file for each additional species in the first column and the NCBI taxonomy ID of that species in the second column:

		   /path/to/species_1.fasta	taxid_1
		   /path/to/species_2.fasta	taxid_2
		   /path/to/species_3.fasta	taxid_3

-  __-f__ &nbsp;&nbsp;&nbsp; Tab-delimited table with additional nucleotide sequences (e.g., non-annotated genome assemblies) to be searched against your query proteins with MMseqs2 in a "tblastn" fashion. This option is extremely useful to improve the detection of orphan genes, since errors in the genome annotations can always lead to the spurious detection of taxonomically-restricted genes (Basile et al., 2019; Weisman et al., 2022). Importantly, each contig/scaffold/chromosome in this dataset should have a unique identifier (e.g., avoid having multiple “>chr1” sequences throughout your genome assemblies). The table format consists of the location of one nucleotide FASTA file for each genome/transcriptome assembly in the first column and the NCBI taxonomy ID of that species in the second column:

		   /path/to/assembly_1.fasta	taxid_1
		   /path/to/assembly_2.fasta	taxid_2
		   /path/to/assembly_3.fasta	taxid_3

-  __-s__ &nbsp;&nbsp;&nbsp; Tab-delimited table with pairwise evolutionary distances, calculated as substitutions/site, between several species in a phylogeny and the query species. This file is necessary to calculate the homology detection failure probabilities of the query genes in step 4 using abSENSE (Weisman et al., 2020). Make sure that the species from which you have evolutionary distances are also represented in the databases that were used for the phylostratigraphic analysis. Extracting the pairwise distances between any pair of species can be easily achieved with the ‘ape’ package in R (Paradis et al., 2004), if you have phylogenetic tree in NEWICK format:
```console
library(ape)
tree<-read.tree(file = "phylogenetic_tree.nwk")
distances<-cophenetic.phylo(x = tree)
write.table(distances,"pairwise_distances.tsv", row.names = TRUE, sep = "\t")
```
&nbsp;&nbsp;&nbsp;&nbsp; Afterwards, the distance table should contain the species taxid in the first column and the distance in substitutions per site in the second column (NOTE that the query species is also included in the table, with an evolutionary distance of 0):

		   query_sp_taxid	0
		   species_1_taxid	distance_1
		   species_2_taxid	distance_2

-  __-i__ &nbsp;&nbsp;&nbsp; Boolean argument that, when true, generates an additional output file with the best sequence hit responsible for the oldest phylostrata assignation for each of the query genes. Established by default as "false" to prevent GenEra from taking a lot of computation time.

Fine-tunning and other useful arguments
=======================================

### These arguments can be modified to suit the user's needs, but keeping the default parameters is usually fine:

-  __-n__ &nbsp;&nbsp;&nbsp; Number of threads to run GenEra. This is the only parameter we HIGHLY SUGGEST to modify, in accordance to the user's needs and resources. Running GenEra is computationally expensive, so using a small amount of threads will result in long running times. By default, GenEra uses 20 threads to run, but we suggest to use as many threads as possible.
-  __-l__ &nbsp;&nbsp;&nbsp; Taxonomic representativeness threshold below which a gene will be flagged as putative genome contamination or the product of a horizontal gene transfer event. The threshold is established as 30% by default, but it can be freely modified by the user. Please refer to GenEra's paper (in prep) for a detailed explanation of taxonomic representativeness. 
-  __-e__ &nbsp;&nbsp;&nbsp; E-value threshold for DIAMOND and MMseqs2, in order to consider a sequence match as a homolog. This value is established as 1e-5 by default, but the user can modify it to a more relaxed or stringent threshold.
-  __-o__ &nbsp;&nbsp;&nbsp; Additional options to feed DIAMOND, based on the user's preferences. This can be useful if, for example, the user wants to filter the homology results based on identity thresholds or query coverage thresholds, instead of the e-value. Users should input the additional commands in quotes, using the original arguments from DIAMOND (for example: -o "--id 30").
-  __-m__ &nbsp;&nbsp;&nbsp; Minimum percentage of matches between your query sequences and another species to consider it useful for the gene age assignation. This parameter only determines whether a phylostratum will be collapsed or not during step 2 of the pipeline. The parameter is implemented so that GenEra collapses the phylostrata where all the representative species lack WGS data in the databases. By default, this threshold is empirically established as 10% of sequence matches with respect to the number of genes in the query species. For example, if the query species has 26,000 genes, and all the representative species for a given phylostratum contain less that 2,600 matches to the query proteins, then the phylostratum is "collapsed" or removed from the analysis.
-  __-x__ &nbsp;&nbsp;&nbsp; Alternative path where the user would like to store the temporary files, as well as the DIAMOND/MMseqs2 results. GenEra generates HUGE temporary files (usually within the range of hundreds of Gigabytes in size), so users might have a hard time storing all these files within the working directory. Thus, the users can redirect the temporary files to a location where there is enough storage space for GenEra to run correctly. By default, the files are stored stored in a temporary directory within the working directory (tmp\_[TAXID]\_[RANDOMNUM]/), which is automatically created by GenEra.
-  __-y__ &nbsp;&nbsp;&nbsp; Modify the sensitivity parameter in DIAMOND. By default, GenEra runs DIAMOND in ultra-sensitive mode to retrieve the highest ammount of homologs in a reasonable amount of time. However, super-sensitive mode can achive similar results but much faster, or the user might want to sacrifice sensitivity in exchange for faster results (please refer to Figure SX of GenEra's paper to make an informed decision while modifying this parameter). 

Output files
============

### The main output files of GenEra are the following:

-  __[TAXID]\_phylostrata\_assignation.tsv__ &nbsp;&nbsp;&nbsp; Tab-delimited table that contains the phylostratigraphic assignation for every gene in the query species, the phylostratum rank that ranges from 1 in the oldest phylostratum (_i.e._, conserved genes throughout all cellular organisms) to the Nth youngest phylostratum (_i.e._, putative orphans at species-level), and the taxonomic representativeness score for each phylostratigraphic assignation. This table can be pared as the input to perform evolutionary transcriptomics through [myTAI](https://github.com/drostlab/myTAI "myTAI").
-  __[TAXID]\_phylostrata\_gene\_count.txt__ &nbsp;&nbsp;&nbsp; Summary file with the number of genes in the query species that could be assigned to each phylostratum.
-  __[TAXID]\_founder\_events.tsv__ &nbsp;&nbsp;&nbsp; Tab-delimited table that contains the oldest phylostratigraphic assignation for each gene family (as defined by MCL clustering), with its respective phylostratum rank and with the number of genes that are contained within each gene family. These could be regarded as the putative gene-family founder events (_i.e._, the point in time where gene families are expected to have originared _de novo_ from a single ancestral gene).
-  __[TAXID]\_founder\_summary.tsv__ &nbsp;&nbsp;&nbsp; Summary file with the number of putative gene-family founder events per phylostratum.
-  __[TAXID]\_abSENSE\_results/Detection\_failure\_probabilities (Optional)__ &nbsp;&nbsp;&nbsp; Main output file generated by [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") when the user specified a table with pairwise evolutionary distances using __-s__. The file contains the probability of a homolog to be undetected due to fast substitution rates, rather than a _de novo_ emergence event at that point in time. The user should look at the probability values in organisms belonging to phylostrata where the gene is no longer detected. Values close to 0 support the inference of a founder event at that phylostratum, while values close to 1 suggest that the phylostratigraphic assignation of that gene may be older but cannot be traced back due to fast substitution rates. Genes labeled with "not_enough_data" are usually putative orphans with no other homologs to be able to compute detection failure probabilities. The age of these genes can also be verified using complementary approaches such as [TOGA](https://github.com/hillerlab/TOGA "TOGA") or [fagin](https://github.com/arendsee/fagin "fagin").

### Other output files that are relevant:

-   __[TAXID]\_ambiguous\_phylostrata.tsv__ &nbsp;&nbsp;&nbsp; Tab-delimited table with genes that ranked low in taxonomic representativeness (by default, below 30%), which were flagged as potential contaminants in the genome or putative horizontal gene transfer events. The table gives a list of possible phylostrata to which these genes could be assigned.
-   __[TAXID]\_deepest_homolog.tsv (Optional)__ &nbsp;&nbsp;&nbsp; Additional output file that is generated when __-i__ is established as true. The file contains the best sequence hit (as defined by the bitscore value) responsible for the oldest phylostrata assignation for each of the query genes. This file is useful to identify erroneous phylostratigraphic assignments due to false positive matches, and to manually evaluate genes with a low taxonomic representativeness.
-   __[TAXID]\_Diamond\_results.bout__ &nbsp;&nbsp;&nbsp; Homology table generated by [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND") (and [MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") when using __-f__) with all the traceable homologs for each query protein. This is an intermediate file generated in the first step of the pipeline, which can be used with the __-p__ argument, in case the user desires to resume GenEra from step 2 onwards (thus, saving a considerable amount of time). This file is usually HUGE, so it is stored by default as a temporary file within a directory made by GenEra (tmp\_[TAXID]\_[RANDOMNUM]/), but it can also be redirected to any specified location using the __-x__ argument.
-   __[TAXID]\_ncbi\_lineages\_yyyy-mm-dd.csv__ &nbsp;&nbsp;&nbsp; A modified version of the lineage table generated by [NCBItax2lin](https://github.com/zyxue/ncbitax2lin "NCBItax2lin"), with the phylostrata ordered in accordance to the query species, and without the phylostrata that lack genomic data for a reliable age assignation. This is an intermediate file generated in the second step of the pipeline, which can be used with the __-c__ argument, in case the user desires to resume GenEra while skipping step 2 (thus, saving some time).
-   __[TAXID]\_abSENSE\_results/[everything\_else] (Optional)__ &nbsp;&nbsp;&nbsp; Additional files generated by [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") when the user specified a table with pairwise evolutionary distances using __-s__. These files contain the bitscore predictions that were used used to calculate the homology detection failure probabilities, as well as other parameters and general information about the analysis. Please refer to the [abSENSE README](https://github.com/caraweisman/abSENSE/blob/master/README.md "abSENSE README") for more detailed information.
