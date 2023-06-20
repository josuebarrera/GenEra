

[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![DOI](https://zenodo.org/badge/483209866.svg)](https://zenodo.org/badge/latestdoi/483209866)
[![Paper link](https://img.shields.io/badge/Published_in-Genome_Biology-badge?labelColor=%23CBE059&color=%231D3050)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02895-z)


![GenEra](https://github.com/josuebarrera/GenEra/blob/main/logo.png)

Introduction
============

GenEra is an easy-to-use and highly customizable command-line tool that estimates gene-family founder events (i.e., the age of the last common ancestor of protein-coding gene families) through the reimplementation of genomic phylostratigraphy (Domazet-Lošo et al., 2007). GenEra takes advantage of [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND")’s speed and sensitivity to search for homolog genes throughout the entire NR database, and combines these results with the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy "NCBI Taxonomy") to assign an origination date for each gene and gene family in a query species. GenEra can also incorporate protein data from external sources to enrich the analysis, it can search for proteins within nucleotide data (i.e., genome/transcriptome assemblies) using [MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") to improve the classification of orphan genes, and it calculates a taxonomic representativeness score to assess the reliability of assigning a gene to a specific age. Additionally, GenEra can calculate homology detection failure probabilities using [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") to help distinguish fast-evolving genes from high-confidence gene-family founder events. 

As of v1.1.0, users can now use [Foldseek](https://github.com/steineggerlab/foldseek "Foldseek") to search protein structural predictions against the [AlphaFold DB](https://alphafold.ebi.ac.uk/ "AlphaFold DB") for fast and sensitive structural alignments. Alternatively, the user can choose to perform a reassessment of gene ages by running [JackHMMER](http://hmmer.org/ "JackHMMER") on top of DIAMOND (be aware, this additional step significantly slows down the analysis). 

Contents
========

-   [Dependencies](#dependencies)
-   [Installation](#installation)
-   [Setting up the databases (if using DIAMOND)](#setting-up-the-databases-if-using-diamond)
-   [Setting up the databases (if using Foldseek)](#setting-up-the-databases-if-using-foldseek)
-   [Quick start for the impatient](#quick-start-for-the-impatient)
-   [Arguments and input files](#arguments-and-input-files)
-   [Fine-tunning and other useful arguments](#fine-tunning-and-other-useful-arguments)
-   [Output files](#output-files)
-   [Citations](#citations)

Dependencies
============

GenEra requires the following software dependencies:

-	[DIAMOND v2.0.0 or higher](https://github.com/bbuchfink/diamond "DIAMOND")
-	[Foldseek v3.915ef7d or higher](https://github.com/steineggerlab/foldseek "Foldseek")
-	[NCBItax2lin](https://github.com/zyxue/ncbitax2lin "NCBItax2lin")
-	[MCL](https://github.com/micans/mcl "MCL")
-	[MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") (optional for protein-against-nucleotide sequence search)
-	[abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") (optional to calculate homology detection failure probabilities)
-	[NumPy](https://numpy.org/ "NumPy") and [SciPy](https://scipy.org/ "SciPy") (needed to run abSENSE in step 4)
-	[R](https://www.r-project.org/ "R") alongside the libraries [optparse](https://cran.r-project.org/web/packages/optparse/index.html "optparse"), [Bio3D](http://thegrantlab.org/bio3d/ "Bio3D"), [Tidyverse](https://www.tidyverse.org/ "Tidyverse"), and [SeqinR](https://cran.r-project.org/web/packages/seqinr/index.html "SeqinR") (optional for JackHMMER reassessment)

Additionally, GenEra requires access to the taxonomy dump from the NCBI and either a locally installed NR database for DIAMOND or a locally installed AlphaFold database for Foldseek.

Installation
============

### Docker installation

A [Docker image](https://hub.docker.com/r/josuebarrera/genera "docker image") can be pulled using the following command:

```console
docker pull josuebarrera/genera
```
GenEra can then be used through Docker by running the following command: 
```console
docker run --rm -v $(pwd):/working-dir -w /working-dir josuebarrera/genera genEra -q [query_sequences.fasta] -t [query_taxid] -b [path/to/nr] -d [path/to/taxdump]
```

### Conda installation

__NOTE: We will release a single-command conda installation soon!!!__

For a manual [conda](https://docs.conda.io/en/latest/ "conda") installation, copy and paste this in your terminal:

```console
git clone https://github.com/josuebarrera/GenEra.git && cd GenEra
chmod +x genEra && chmod +x Erassignment && chmod +x hmmEra
git clone https://github.com/caraweisman/abSENSE.git && mv abSENSE/Run_abSENSE.py . && chmod +x Run_abSENSE.py
conda create -n genEra python=3.7
conda activate genEra
conda install -c bioconda -c conda-forge diamond
conda install -c bioconda mcl
pip install -U ncbitax2lin
conda install -c conda-forge -c bioconda mmseqs2
conda install -c conda-forge -c bioconda foldseek
conda install -c anaconda scipy
conda install -c bioconda r-seqinr
conda install -c bioconda r-optparse
conda install -c bioconda r-bio3d
conda install -c r r-tidyverse
CONDABIN=$(which ncbitax2lin | sed 's/ncbitax2lin//g') && mv genEra ${CONDABIN} && mv Erassignment ${CONDABIN} && mv Run_abSENSE.py ${CONDABIN} && mv hmmEra ${CONDABIN}
```

Once you finished, you can test whether genEra was installed correctly by running:

```console
bash test_installation.sh
```

If all the tests appear as PASSED, then you are ready to use the software!

Setting up the databases (if using DIAMOND)
===========================================

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
wget -N ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
mkdir -p taxdump && tar zxf new_taxdump.tar.gz -C ./taxdump
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

Setting up the databases (if using Foldseek)
============================================

Foldseek uses 3D structure predictions in PDB format as input. So first make sure that you have structural predictions for each protein of your query species. This can be done using tools such as [AlphaFold](https://github.com/deepmind/alphafold "AlphaFold") or [OmegaFold](https://github.com/HeliXonProtein/OmegaFold "OmegaFold"). for example:
```console
omegafold query_sequences.fasta output_directory 
```
The folding prediction of each protein should be stored on independent PDB files within a single directory. Make sure all the PDB files are uncompressed before running the analysis! 

Once you have that, use Foldseek to download the AlphaFold database (warning: this database is huge):
```console
foldseek databases Alphafold/UniProt alphafoldDB tmp 
```
Finally, download the taxonomy dump from the NCBI:
```console
wget -N ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
mkdir -p taxdump && tar zxf new_taxdump.tar.gz -C ./taxdump
```

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
If you prefer to use Foldseek instead of DIAMOND, the basic usage is:
```console
genEra -Q [query_structure_directory] -t [query_taxid] -B [path/to/alphafoldDB] -d [path/to/taxdump]
```
And of course, it can also be optimized to get better results:
```console
genEra -Q [query_structure_directory] -t [query_taxid] -B [path/to/alphafoldDB] -d [path/to/taxdump] \ 
-A [structural_prediction_list.tsv] -s [evolutionary_distances.tsv] -i [true] -n [many threads]
  ```

Arguments and input files
=========================

### GenEra always requires this input:

-  `-t` &nbsp;&nbsp;&nbsp; The NCBI taxonomy ID of your species of interest (it can be easily found at https://www.ncbi.nlm.nih.gov/taxonomy).

### Depending on which alignment method you want to use, you'll have to choose between these two input files:

-  `-q` &nbsp;&nbsp;&nbsp; A standard FASTA file containing all the protein sequences of your species of interest. Make sure all your sequence headers are unique and make sure to avoid regular expressions like `\t` or `\n` in the sequence headers, as these will interfere with the pipeline.

-  `-Q` &nbsp;&nbsp;&nbsp; A directory containing all the protein structural predictions of your species of interest in PDB format. Make sure that all your PDB files have unique names and that they are uncompressed. Also make sure to avoid regular expressions like `\t` or `\n` in the PDB file names, as these will interfere with the pipeline.

### Additionally, GenEra requires one of the following three database files:

-  `-b` &nbsp;&nbsp;&nbsp; A locally installed NR database for DIAMOND (or any other database whose sequences can be automatically traced back to NCBI taxids). This database option is only compatible with FASTA input files (`-q`).

-  `-B` &nbsp;&nbsp;&nbsp; A locally installed AlphaFold database for Foldseek (you can choose between several options that can be automatically downloaded using Foldseek). This database option is only compatible with PDB input files (`-Q`).

-  `-p` &nbsp;&nbsp;&nbsp; A pre-generated DIAMOND/Foldseek table. This file is automatically generated by GenEra in the first step of the pipeline, either by choosing DIAMOND or Foldseek, and can be given to GenEra in case the user decides to re-run the pipeline (e.g., fine-tunning the parameters). The first column contains the query sequence IDs (qseqid), the second coulmn contains the target sequence IDs (sseqid), the third column contains the e-value of the matching alignment (evalue), the fourth column contains the bitscore of the alignment (bitscore) and the fifth column contains the taxonomy ID of the target sequence (staxids). 

### And one of the following three files:

-  `-d` &nbsp;&nbsp;&nbsp; The location of the NCBI taxonomy files, colloquially known as the `taxdump`. Make sure the taxdump directory contains the files `nodes.dmp`, `names.dmp` and `fullnamelineage.dmp`. These files are used by NCBItax2lin to create a `ncbi_lineages` file and for GenEra to organize this file in the correct hierarchical order.

-  `-r` &nbsp;&nbsp;&nbsp; The location of the uncompressed `ncbi_lineages` file generated by NCBItax2lin. Once the user runs GenEra for the first time (or runs NCBItax2lin independently), this file can be used to save some time during step 2 of the pipeline. This file is a comma-delimited table with each row representing a specific lineage, and each column representing the taxonomic hierarchies to which that lineage belongs. This file will be automatically modified by genEra to rearrange the columns in hierarchical order using the `fullnamelineage.dmp` found in the `taxdump`. Alternatively, GenEra will try to download the lineage information from the NCBI Taxonomy webpage. IMPORTANT: Some users have experienced memory errors while using the `-d` argument. If your analysis stops after this is written in the STDOUT `Preparings all lineages into a dataframe to be written to disk ...`, you are likely experiencing the same issue, which is related to the dependency ncbitax2lin. [We uploaded a compressed "ncbi_lineages" file](https://github.com/josuebarrera/GenEra/tree/main/ncbi_lineages_raw_file) that can be used with `-r` after uncompressing it. This should give you the exact same result as if you ran the pipeline from scratch with `-d`.

-  `-c` &nbsp;&nbsp;&nbsp; Custom `ncbi_lineages` file that is already tailored for the query species. GenEra modifies the raw `ncbi_lineages` file so that the taxid of the query species appears in the first column, and the taxonomic hierarchies of the query species are rearranged from the species level all the way back to the last universal common ancestor (termed “cellular organisms” by the NCBI taxonomy). The file is also modified by GenEra to collapse the phylostrata (i.e., the taxonomic levels) that lack the necessary genomic data to be useful in the analysis. Once this file is generated, the user can re-use this file if they want to run GenEra on the same species with different parameters. By default, GenEra will search for the correct hierarchical order of the query organism’s phylostrata directly from the NCBI webpage. If the user machine is unable to access the NCBI webpage, GenEra will attempt to infer the correct order of the phylostrata directly from the `ncbi_lineages` file. However, there are some taxonomic hierarchies that the NCBI labels as “clades” (e.g., Embryophyta in the plant lineage), which GenEra cannot automatically assign to their correct taxonomic level when working offline. This option is also useful if the user wants to modify the taxonomic hierarchies that were originally assigned by the NCBI (e.g., an outdated taxonomy that does not reflect the phylogenetic relationships of the query species). This custom table can be easily generated by printing the desired columns of the `ncbi_lineages` file with awk:
```console
awk -F'"' -v OFS='' '{ for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' ncbi_lineages_[date].csv | awk -F "," '{ print $1","$8","$7","$6"," ... }' > custom_table.csv
```
### The user can also incorporate five other inputs, which are optional but can potentially improve the final age assignment:

-  `-a` &nbsp;&nbsp;&nbsp; Tab-delimited table with additional proteins to be included in the analysis __(only compatible when using FASTA sequences as input)__. This option is particularly useful if the user wants to include proteins from species that are absent from the nr database, such as newly annotated genomes or transcriptomes. It can also help fill the gaps of phylostrata that would otherwise be collapsed due to lack of available genomes in the database. Importantly, each protein in this custom dataset should have unique identifiers, otherwise GenEra will not work properly during the taxonomic assignment of the Diamond hits. The table format consists of the location of one protein FASTA file for each species in the first column and the NCBI taxonomy ID of that species in the second column:

		/path/to/species_1.fasta	taxid_1
		/path/to/species_2.fasta	taxid_2
		/path/to/species_3.fasta	taxid_3

-  `-f` &nbsp;&nbsp;&nbsp; Tab-delimited table with additional nucleotide sequences (e.g., non-annotated genome assemblies) to be searched against your query proteins with MMseqs2 in a "tblastn" fashion __(only compatible when using FASTA sequences as input)__. This option is extremely useful to improve the detection of orphan genes, since errors in the genome annotations can oftentimes lead to the spurious detection of taxonomically-restricted genes ([Basile et al., 2019](https://doi.org/10.1101/185983 "Basile et al., 2019"); [Weisman et al., 2022](https://doi.org/10.1016/j.cub.2022.04.085 "Weisman et al., 2022")). Importantly, each contig/scaffold/chromosome in this dataset should have a unique identifier (e.g., avoid having multiple `>chr1` sequences throughout your genome assemblies). The table format consists of the location of one nucleotide FASTA file for each genome/transcriptome assembly in the first column and the NCBI taxonomy ID of that species in the second column:

		/path/to/assembly_1.fasta	taxid_1
		/path/to/assembly_2.fasta	taxid_2
		/path/to/assembly_3.fasta	taxid_3

-  `-A` &nbsp;&nbsp;&nbsp; Tab-delimited table with additional structural predictions to be included in the analysis __(only compatible when using PDB files as input)__. This option is useful if the user wants to include protein folding predictions from species that are absent from the AlphaFold database, such as newly annotated genomes. It can also help fill the gaps of phylostrata that would otherwise be collapsed due to lack of available genomes in the database. Importantly, each PDB file in this custom dataset should have a unique identifier, both within and between species, otherwise GenEra will not work properly during the taxonomic assignment of the Foldseek hits. The table format consists of the location of one directory containing all the uncompressed PDB files for each species in the first column and the NCBI taxonomy ID of that species in the second column:

		/path/to/species_directory_1	taxid_1
		/path/to/species_directory_2	taxid_2
		/path/to/species_directory_3	taxid_3

-  `-s` &nbsp;&nbsp;&nbsp; Tab-delimited table with pairwise evolutionary distances, calculated as substitutions/site, between several species in a phylogeny and the query species. This file is necessary to calculate the detection failure probabilities of the query genes in step 4 using abSENSE ([Weisman et al., 2020](https://doi.org/10.1371/journal.pbio.3000862 "Weisman et al., 2020")), which allows GenEra to find gene age assignments and gene-family founder events that cannot be explained by homology detection failure (HDF). Please note that HDF-tested gene age assignments can only be calculated whenever the user gives GenEra evolutionary distances of species that are found in the taxonomic level that inmediately precedes the one that is being studied (_i.e._, the next older taxonomic level), meaning that genEra will not be able to evaluate DHF-tested assignments for the oldestaxonomic level that contains species representatives in the file with evolutionary distances, nor for taxonomic levels that lack an outgroup species in the file with evolutionary distances. Also note that HDF-tested gene age assignments cannot be calculated for species-specific genes, reason being that detection failure probabilities cannot be calculated for genes that lack any traceable homologs (please read [Weisman et al., 2020](https://doi.org/10.1371/journal.pbio.3000862 "Weisman et al., 2020") for more details). Finally, please make sure that the species from which you have evolutionary distances are also represented in the databases that are being used for the homology search using either DIAMOND or Foldseek. In order to calculate reliable detection failure probabilities, it is important to have as many closely-related species as posible in the evolutionary distance file. Extracting the pairwise distances between any pair of species can be easily achieved with the [ape package in R](http://ape-package.ird.fr/ "ape package in R") ([Paradis et al., 2019](https://doi.org/10.1093/bioinformatics/bty633 "Paradis et al., 2019")), if you have phylogenetic tree in NEWICK format:
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

-  `-i` &nbsp;&nbsp;&nbsp; Boolean argument that, when `true`, generates an additional output file with the best sequence/structure hit responsible for the oldest gene age assignment for each of the query genes. This file is mainly used to verify whether the oldest age assignment of genes are not driven by false positive matches to the database. Established by default as `false` to save computing time.

Fine-tunning and other useful arguments
=======================================

### These arguments can be modified to suit the user's needs, but keeping the default parameters is usually fine:

-  `-n` &nbsp;&nbsp;&nbsp; Number of threads to run GenEra. This is the only parameter we HIGHLY SUGGEST to modify, in accordance to the user's needs and resources. Running GenEra is computationally expensive, so using a small amount of threads will result in long running times. By default, GenEra uses `20` threads to run, but we suggest to use as many threads as possible.

-  `-o` &nbsp;&nbsp;&nbsp; Path to the directory where you want to store the final output files (GenEra will automatically create this directory if it doesn't exist). Many users like to keep things tidy and separate the output files to an specific location in their computer. By default, GenEra stores the output files within the working directory.

-  `-l` &nbsp;&nbsp;&nbsp; Taxonomic representativeness threshold below which a gene will be flagged as putative genome contamination or the product of a horizontal gene transfer event. The threshold is established as 30% by default, but it can be freely modified by the user. Please refer to [GenEra's preprint](https://doi.org/10.1101/2022.07.07.498977 "GenEra paper") for a detailed explanation of taxonomic representativeness. 

-  `-e` &nbsp;&nbsp;&nbsp; E-value threshold for DIAMOND or Foldseek, in order to consider a sequence match as a homolog. This value is established as `1e-5` by default, but the user can modify it to a more relaxed or stringent threshold.

-  `-u` &nbsp;&nbsp;&nbsp; Additional options to feed DIAMOND or Foldseek, based on the user's preferences. This can be useful if, for example, the user wants to filter the homology results based on identity thresholds or query coverage thresholds, instead of the e-value. Users should input the additional commands in quotes, using the original arguments from DIAMOND (for example: `-u "--id 30"`) or from Foldseek (for example: `-u "-c 0.8 --cov-mode 0"`).

-  `-m` &nbsp;&nbsp;&nbsp; Minimum percentage of matches between your query sequences and another species to consider it useful for the gene age assignment. This parameter only determines whether a taxonomic level will be collapsed or not during step 2 of the pipeline. The parameter is implemented so that GenEra collapses the taxonomic levels where all the representative species lack WGS data in the databases. By default, this threshold is empirically established as 10% of sequence matches with respect to the number of genes in the query species. For example, if the query species has 26,000 genes, and all the representative species for a given phylostratum contain less that 2,600 matches to the query proteins, then the taxonomic level is "collapsed" or removed from the analysis.

-  `-x` &nbsp;&nbsp;&nbsp; Alternative path where the user would like to store the temporary files, as well as the DIAMOND/Foldseek results. GenEra generates HUGE temporary files (usually within the range of dozens to hundreds of Gigabytes in size), so users might have a hard time storing all these files within the working directory. Thus, the users can redirect the temporary files to a location where there is enough storage space for GenEra to run correctly. By default, the files are stored stored in a temporary directory within the working directory (`tmp_[TAXID]_[RANDOMNUM]/`), which is automatically created by GenEra.

-  `-y` &nbsp;&nbsp;&nbsp; Modify the sensitivity parameter in DIAMOND. By default, GenEra runs DIAMOND in sensitive mode to retrieve the highest ammount of homologs in a reasonable amount of time. Users can decide to run DIAMOND in ultra-sensitive mode to achieve slighlty higher sensitivity at the cost of speed, or the users might prefer to sacrifice sensitivity in exchange for faster results by using fast mode (please refer to Figure 3 of [GenEra's paper](https://doi.org/10.1186/s13059-023-02895-z "GenEra paper") to make an informed decision while modifying this parameter). 

-  `-M` &nbsp;&nbsp;&nbsp; Adjust the amount of prefilter made by Foldseek (i.e., the maximum number of hits that are reported in the output file). Unlike DIAMOND, Foldseek lacks the option to print an unlimited amount of sequence hits, which may hamper the hability of the user to trace back very deep evolutionary homologs for the query proteins. We prefer to keep this parameter as high as possible, although this desicion is usually limited by the amout of RAM and storage that Foldseek will use. By default, GenEra limits the maximum number of hits to 10000, but the user can modify this parameter in accordance with their needs and hardware limitations.

-  `-j` &nbsp;&nbsp;&nbsp; When true, GenEra runs an additional search step using JackHMMER against the online Uniprot database for all the genes that could not be traced back to "cellular organisms" using DIAMOND. For this, GenEra uses the R package [Bio3D](http://thegrantlab.org/bio3d/ "Bio3D") to submit individual jobs to the EMBL-EBI server. We decided to use this method because the local version of HMMER does not output taxids (which are crucial for GenEra), while the online version does. __This severely limits the speed of this search step, since it cannot be multithreaded. Additionally, it requires constant access to the internet to work, which might complicate things.__ Nonetheless, we implemented this step for users that prefer to search for protein homologs using HMM-profile methods over pairwise alignment methods.

-  `-k` &nbsp;&nbsp;&nbsp; The user can modify this parameter to decide the starting taxonomic level of the genes that will be re-analyzed when running JackHMMER with `-j`. By default, GenEra runs JackHMMER over every gene that could be traced back from the second oldest taxonomic level onwards, namely the taxonomic rank 2. However, genes that were assigned within the second oldest taxonomic rank (e.g., Eukaryota, Bacteria, Archaea) will usually have a graet amount os sequence hits using JackHMMER, which slows down the process a lot. Therefore, the user can decide to limit the JackHMMER searches to younger taxonomic levels to speed up the process and make the computation times feasible.

Output files
============

### The main output files of GenEra are the following:

-  `[TAXID]_gene_ages.tsv` &nbsp;&nbsp;&nbsp; Tab-delimited table that contains the age assignment for every gene in the query species, the taxonomic rank that ranges from 1 in the oldest taxonomic level (_i.e._, conserved genes throughout all cellular organisms) to the Nth youngest taxonomic level (_i.e._, putative orphans at species-level), and the taxonomic representativeness score for each gene age assignment. This table can be used as the input to perform evolutionary transcriptomics through [myTAI](https://github.com/drostlab/myTAI "myTAI").

-  `[TAXID]_gene_age_summary.tsv` &nbsp;&nbsp;&nbsp; Summary file with the number of genes in the query species that could be assigned to each taxonomic level.

-  `[TAXID]_founder_events.tsv` &nbsp;&nbsp;&nbsp; Tab-delimited table that contains the oldest age assignment for each gene family (as defined by MCL clustering), with its respective taxonomic rank and with the number of genes that are contained within each gene family. These could be regarded as the putative gene-family founder events (_i.e._, the point in time where gene families are expected to have originared).

-  `[TAXID]_founder_summary.tsv` &nbsp;&nbsp;&nbsp; Summary file with the number of putative gene-family founder events per taxonomic level.

-  `[TAXID]_HDF_gene_ages.tsv` __(Optional)__ &nbsp;&nbsp;&nbsp; Tab-delimited table that contains the genes whose age assignment cannot be explained by homology detection failure (HDF). This file is created by GenEra whenever the user specifies a table with pairwise evolutionary distances using `-s`. The genes are selected based on the detection failure probabilities (calculated with [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE")) of the closest outgroup for each given taxonomic level in the analysis. All the genes whose detection failure probabilities are lower than 0.05 in the closest outgroup are deemed as gene age assignments that passed the HDF test.

-  `[TAXID]_HDF_gene_age_summary.tsv` __(Optional)__ &nbsp;&nbsp;&nbsp; Summary file with the number of gene-age assignments that passed the HDF test for each taxonomic level. This file is created by GenEra whenever the user specifies a table with pairwise evolutionary distances using `-s`. GenEra will assign an NA in the gene count of all the taxonomic levels that lacked an outgroup in the file with evolutionary distances. Species-specific genes are also treated as NA, since detection failure probabilities cannot be calculated for single datapoints. Taxonomic levels with a gene count of 0 means that an appropriate outgroup was available in the analysis, but GenEra could not detect any high-confidence gene for that taxonomic level (sometimes due to the lack of enough datapoints to calulate detection failure probabilities). The fouth column contains the NCBI taxonomy ID of the outgroup species that was used to determine whether these genes passed the HDF test or not.

-  `[TAXID]_HDF_founder_events.tsv` __(Optional)__ &nbsp;&nbsp;&nbsp; Tab-delimited table that contains the oldest age assignment of the gene-families that contain at least one gene that passed the HDF test for that taxonomic level (_i.e._, gene-family founder events whose age assignment cannot be explained by HDF). This file is created by GenEra whenever the user specifies a table with pairwise evolutionary distances using `-s`.

-  `[TAXID]_HDF_founder_summary.tsv` __(Optional)__ &nbsp;&nbsp;&nbsp; Summary file with the number of gene-family founder events that passed the HDF test per taxonomic level. This file is created by GenEra whenever the user specifies a table with pairwise evolutionary distances using `-s`. The fouth column contains the NCBI taxonomy ID of the outgroup species that was used to calculate determine whether these gene families passed the HDF test or not.

### Other output files that are relevant:

-   `[TAXID]_ambiguous_phylostrata.tsv` &nbsp;&nbsp;&nbsp; Tab-delimited table with genes that ranked low in taxonomic representativeness (by default, below 30%), which were flagged as potential contaminants in the genome or putative horizontal gene transfer events. The table gives a list of possible taxonomic levels to which these genes could be assigned.

-   `[TAXID]_deepest_homolog.tsv` __(Optional)__ &nbsp;&nbsp;&nbsp; Additional output file that is generated when `-i` is established as `true`. The file contains the best sequence hit (as defined by the bitscore value) responsible for the oldest gene age assignment for each of the query genes. This file is useful to identify erroneous age assignments due to false positive matches, and to manually evaluate genes with a low taxonomic representativeness.

-   `[TAXID]_Diamond_results.bout` &nbsp;&nbsp;&nbsp; Homology table generated by [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND") (and [MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") when using `-f`) with all the traceable homologs for each query protein (this file is only generated when using a FASTA file as input). This is an intermediate file generated in the first step of the pipeline, which can be used with the `-p` argument, in case the user desires to resume GenEra from step 2 onwards (thus, saving a considerable amount of time). This file is usually HUGE, so it is stored by default as a temporary file within a directory made by GenEra (`tmp_[TAXID]_[RANDOMNUM]/`), but it can also be redirected to any specified location using the `-x` argument.

-   `[TAXID]_Foldseek_results.bout` &nbsp;&nbsp;&nbsp; Homology table generated by [Foldseek](https://github.com/steineggerlab/foldseek "Foldseek") with all the traceable homologs for each query protein structure (this file is only generated when using PDB files as input). This is an intermediate file generated in the first step of the pipeline, which can be used with the `-p` argument, in case the user desires to resume GenEra from step 2 onwards (thus, saving a considerable amount of time). This file is usually big, so it is stored by default as a temporary file within a directory made by GenEra (`tmp_[TAXID]_[RANDOMNUM]/`), but it can also be redirected to any specified location using the `-x` argument.

-   `[TAXID]_HMMER_results.bout` __(Optional)__ &nbsp;&nbsp;&nbsp; Homology table generated by [JackHMMER](http://hmmer.org/ "JackHMMER") with all the traceable homologs of all the query proteins that were re-analyzed (this file is only generated when `-j` is changed to `true`). This file is concatenated on to of the DIAMOND results (`[TAXID]_Diamond_results.bout`), which is then used by GenEra to re-calculate gene ages. This file is stored by default as a temporary file within a directory made by GenEra (`tmp_[TAXID]_[RANDOMNUM]/`), but it can also be redirected to any specified location using the `-x` argument.

-   `[TAXID]_ncbi_lineages.csv` &nbsp;&nbsp;&nbsp; A modified version of the lineage table generated by [NCBItax2lin](https://github.com/zyxue/ncbitax2lin "NCBItax2lin"), with the phylostrata ordered in accordance to the query species, and without the phylostrata that lack genomic data for a reliable age assignment. This is an intermediate file generated in the second step of the pipeline, which can be used with the `-c` argument, in case the user desires to resume GenEra while skipping step 2 (thus, saving some time).

-   `[TAXID]_abSENSE_results/` __(Optional)__ &nbsp;&nbsp;&nbsp; Files generated by [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") when the user specified a table with pairwise evolutionary distances using `-s`. These files include a table with all the calculated detection failure probabilities, the bitscore predictions that were used used to calculate these probabilities, as well as other parameters and general information about the analysis. Please refer to the [abSENSE README](https://github.com/caraweisman/abSENSE/blob/master/README.md "abSENSE README") for more detailed information.

Citations
=========

The paper describing the method implemented in GenEra:
```console
Barrera-Redondo, J., Lotharukpong, J.S., Drost, H.G., Coelho, S.M. (2023). Uncovering gene-family founder events during major evolutionary transitions in animals, plants and fungi using GenEra. Genome Biology, 24, 54. https://doi.org/10.1186/s13059-023-02895-z
```
GenEra makes use of several dependencies that should also be cited, if implemented within the pipeline.

For a standard run:
```console
Buchfink, B., Reuter, K., Drost, H.G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature methods, 18(4), 366-368.
Enright, A.J., Van Dongen, S., Ouzounis, C.A. (2002). An efficient algorithm for large-scale detection of protein families. Nucleic acids research, 30(7), 1575-1584.
```
When implementing protein-against-genome alignments: 
```console
Steinegger, M., Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028.
```
If using homology detection failure probabilities:
```console
Weisman, C.M., Murray, A.W., Eddy, S.R. (2020). Many, but not all, lineage-specific genes can be explained by homology detection failure. PLoS biology, 18(11), e3000862.
```
When using structural alignments:
```console
van Kempen, M., Kim, S., Tumescheit, C., Mirdita, M., Söding, J., & Steinegger, M. (2022). Foldseek: fast and accurate protein structure search. bioRxiv.
Varadi, M., Anyango, S., Deshpande, M., Nair, S., Natassia, C., Yordanova, G., ... & Velankar, S. (2022). AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic acids research, 50(D1), D439-D444.
```
Finally, if using HMMER on top of the standard analysis:
```console
Grant, B. J., Rodrigues, A. P., ElSawy, K. M., McCammon, J. A., & Caves, L. S. (2006). Bio3d: an R package for the comparative analysis of protein structures. Bioinformatics, 22(21), 2695-2696.
Finn, R. D., Clements, J., & Eddy, S. R. (2011). HMMER web server: interactive sequence similarity searching. Nucleic acids research, 39(suppl_2), W29-W37.
Charif, D., & Lobry, J. R. (2007). SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. In Structural approaches to sequence evolution (pp. 207-232). Springer, Berlin, Heidelberg.
```
