# *Drosophilidae* Virus Phylogenetic Correlations

This repository contains data and code used in Imrie et al., (2024) Evidence of positive correlations in susceptibility to a diverse panel of viruses across *Drosophilidae* hosts.

### Files
<ul>
  <li>MCMCglmms.R</li>
  <ul>
    <li>This file contains code to run all univariate and bivariate phylogenetic MCMCglmms used in this study.</li>
  </ul>
<li>data/data.vials</li>
<ul>
  <li>This file contains fold-changes in viral load and other information for each biological replicate (vial) used in this study. </li>
</ul>
  <li>data/phylogeny_host.nwk</li>
<ul>
  <li>This file contains the maximum clade credibility <em>Drosophilidae</em> host phylogeny used to fit phylogenetic MCMCglmms.</li>
</ul>
  <li>data/phylogeny_dicistrovirus.tree</li>
<ul>
  <li>This file contains the maximum clade credibility dicistrovirus phylogeny used to select DCV and CrPV isolates for inclusion in this study.</li>
</ul>
    <li>data/dicistrovirus_NC_concatenated.fas</li>
<ul>
  <li>This file contains the concatenated alignment of DCV and CrPV isolate non-coding genome sequences. Each non-coding region was aligned separately using MUSCLE V5.1 before concatenating.</li>
</ul>
    <li>data/dicistrovirus_ORF_concatenated.fas</li>
<ul>
  <li>This file contains the concatenated alignment of DCV and CrPV isolate coding genome sequences. Each coding region was translation aligned separately using MUSCLE V5.1 before concatenating.</li>
</ul>
<li>data/dicistrovirus_genome_annotation.gb</li>
<ul>
  <li>This file contains annotated genomes of DCV and CrPV isolates indicating the boundaries of coding and non-coding regions separated during alignment.</li>
</ul>
</ul>
<li>data/dicistrovirus_BEAST.xml</li>
<ul>
  <li>This file contains the BEAST run settings used to infer the dicistrovirus phylogeny.</li>
</ul>
</ul>
