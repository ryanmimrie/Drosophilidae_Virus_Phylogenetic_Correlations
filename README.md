# *Drosophilidae* Virus Phylogenetic Correlations

This repository contains data and code used in Imrie et al., (2024) No Evidence of Trade-Offs in Susceptiblity to Different Viruses across *Drosophilidae* hosts.

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
  <li>This file contains the ultrametric <em>Drosophilidae</em> host phylogeny used to fit phylogenetic MCMCglmms.</li>
</ul>
  <li>data/phylogeny_cripavirus.nwk</li>
<ul>
  <li>This file contains the Cripavirus phylogeny used to select DCV and CrPV isolates for inclusion in this study.</li>
</ul>
    <li>data/cripavirus_NC_concatenated.fas</li>
<ul>
  <li>Concatenated alignment of DCV and CrPV isolate non-coding genome sequences. Each non-coding region was aligned separately using MUSCLE V5.1 before concatenating.</li>
</ul>
    <li>data/cripavirus_NC_concatenated.fas</li>
<ul>
  <li>Concatenated alignment of DCV and CrPV isolate coding genome sequences. Each coding region was translation aligned separately using MUSCLE V5.1 before concatenating.</li>
</ul>
<li>data/cripavirus_genome_annotation.gb</li>
<ul>
  <li>Annotated genomes of DCV and CrPV isolates indicating the boundaries of coding and non-coding regions separated during alignment.</li>
</ul>
</ul>
