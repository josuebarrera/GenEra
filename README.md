

[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![DOI](https://zenodo.org/badge/483209866.svg)](https://zenodo.org/badge/latestdoi/483209866)
[![Paper link](https://img.shields.io/badge/Published_in-Genome_Biology-badge?labelColor=%23CBE059&color=%231D3050)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02895-z)
![Visitors](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fjosuebarrera%2FGenEra&label=VISITORS&countColor=%23263759&style=flat)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/genera) [![Downloads](https://anaconda.org/bioconda/genera/badges/downloads.svg)](https://anaconda.org/bioconda/genera)

![GenEra](https://github.com/josuebarrera/GenEra/blob/main/logo.png)

Introduction
============

`GenEra` is an easy-to-use and highly customizable command-line tool that estimates gene-family founder events (i.e., the age of the last common ancestor of protein-coding gene families) through the reimplementation of genomic phylostratigraphy ([Domazet-Lošo et al., 2007](https://www.sciencedirect.com/science/article/pii/S0168952507002995)).
-	`GenEra` takes advantage of [DIAMOND](https://github.com/bbuchfink/diamond "DIAMOND")’s speed and sensitivity to search for homolog genes throughout the entire NR database, and combines these results with the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy "NCBI Taxonomy") to assign an origination date for each gene and gene family in a query species.
-	`GenEra` can also incorporate protein data from external sources to enrich the analysis, it can detect very recently evolved proteins by incorporating different strains os varieties within the same species, it can search for proteins within nucleotide data (i.e., genome/transcriptome assemblies) using [MMseqs2](https://github.com/soedinglab/MMseqs2 "MMseqs2") to improve the classification of orphan genes, and it calculates a taxonomic representativeness score to assess the reliability of assigning a gene to a specific age.
-	Additionally, `GenEra` can calculate homology detection failure probabilities using [abSENSE](https://github.com/caraweisman/abSENSE "abSENSE") to help distinguish fast-evolving genes from high-confidence gene-family founder events. 

As of v1.1.0, users can now use [Foldseek](https://github.com/steineggerlab/foldseek "Foldseek") to search protein structural predictions against the [AlphaFold DB](https://alphafold.ebi.ac.uk/ "AlphaFold DB") for fast and sensitive structural alignments. Alternatively, the user can choose to perform a reassessment of gene ages by running [JackHMMER](http://hmmer.org/ "JackHMMER") on top of DIAMOND (be aware, this additional step significantly slows down the analysis).

### Precomputed gene ages (or 'phylomaps') made using `GenEra` or from previous studies using other tools can be found [here](https://github.com/HajkD/published_phylomaps).

[Documentation](https://github.com/josuebarrera/GenEra/wiki)
============

We recommend users to consult the `GenEra` wiki for details on [installation](https://github.com/josuebarrera/GenEra/wiki/Installing-GenEra) (via [Conda](https://github.com/josuebarrera/GenEra/wiki/Installing-GenEra#conda-installation) or [Docker](https://github.com/josuebarrera/GenEra/wiki/Installing-GenEra#docker-installation)), [database setup](https://github.com/josuebarrera/GenEra/wiki/Setting-up-the-database(s)) and [how to run `GenEra`](https://github.com/josuebarrera/GenEra/wiki/Running-GenEra), as well as the [output files](https://github.com/josuebarrera/GenEra/wiki/GenEra-output).
We also discuss potential [downstream analyses](https://github.com/josuebarrera/GenEra/wiki/Downstream-analyses) that can be performed on the GenEra output.

Please [cite](https://github.com/josuebarrera/GenEra/wiki/Citations) the appropriate tools when using the [dependencies](https://github.com/josuebarrera/GenEra/wiki#dependencies) of `GenEra`. These citations are valuable in furthering bioinformatics research.

The paper describing the method implemented in [GenEra](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02895-z):
```console
Barrera-Redondo, J., Lotharukpong, J.S., Drost, H.G., Coelho, S.M. (2023). Uncovering gene-family founder events during major evolutionary transitions in animals, plants and fungi using GenEra. Genome Biology, 24, 54. https://doi.org/10.1186/s13059-023-02895-z
```

Acknowledgement
=========

We ([Josué Barrera-Redondo](https://github.com/josuebarrera), [Jaruwatana Sodai Lotharukpong](https://github.com/LotharukpongJS) & [Hajk-Georg Drost](https://github.com/HajkD)) would like to thank several individuals for making this project possible.

We gratefully thank Susana M. Coelho, the Max Planck Institute for Biology Tübingen and the Max Planck Society for hosting and facilitating this research.
We thank Caroline M. Weisman for her helpful comments on how to analyze and interpret HDF probabilities of her software [abSENSE](https://github.com/caraweisman/abSENSE).
We thank the Max Planck Computing and Data Facility for access to and support of the HPC infrastructure, as well as the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A532B, 031A533A, 031A533B, 031A534A, 031A535A, 031A537A, 031A537B, 031A537C, 031A537D, 031A538A).

Lastly, we are very grateful to Alice Laigle, Erica Dinatale, Laura Piovani, Michael Borg, Alexandra Dallaire and all the early adopters for their testing and feedback.

Funding
=========

This work was supported by the European Research Council Grant “THETYS” (Grant agreement ID 864038), the Alexander von Humboldt Foundation, the Gordon and Betty Moore Foundation, and the Max Planck Society.
