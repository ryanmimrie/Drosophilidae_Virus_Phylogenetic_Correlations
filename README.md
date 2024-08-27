
<img src=img/logo.png" width="200" alt="Logo" style="display: block; margin: auto;">


# *Drosophilidae* Virus Phylogenetic Correlations

This repository contains data and code used in Imrie et al., (2024) "Evidence of positive correlations in susceptibility to a diverse panel of viruses across *Drosophilidae* hosts."

## Files

### Code
- `MCMCglmms.R`: This file contains code to run all univariate and bivariate phylogenetic MCMCglmms used in this study.

### Data
- `data/data.vials`: This file contains fold-changes in viral load and other information for each biological replicate (vial) used in this study.
- `data/phylogeny_host.nwk`: This file contains the maximum clade credibility *Drosophilidae* host phylogeny used to fit phylogenetic MCMCglmms.
- `data/phylogeny_dicistrovirus.tree`: This file contains the maximum clade credibility dicistrovirus phylogeny used to select DCV and CrPV isolates for inclusion in this study.
- `data/dicistrovirus_NC_concatenated.fas`: This file contains the concatenated alignment of DCV and CrPV isolate non-coding genome sequences. Each non-coding region was aligned separately using MUSCLE V5.1 before concatenating.
- `data/dicistrovirus_ORF_concatenated.fas`: This file contains the concatenated alignment of DCV and CrPV isolate coding genome sequences. Each coding region was translation aligned separately using MUSCLE V5.1 before concatenating.
- `data/dicistrovirus_genome_annotation.gb`: This file contains annotated genomes of DCV and CrPV isolates indicating the boundaries of coding and non-coding regions separated during alignment.
- `data/dicistrovirus_BEAST.xml`: This file contains the BEAST run settings used to infer the dicistrovirus phylogeny.
