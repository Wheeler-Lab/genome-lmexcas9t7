[![DOI](https://zenodo.org/badge/544940552.svg)](https://zenodo.org/badge/latestdoi/544940552)
# _L. mexicana_ genome remapping, polishing and de-novo assembly code

## 1_remap_Fiebig_et_al

This folder contains code to take the genome annotations in [Fiebig et al. (2015)](https://doi.org/10.1371/journal.ppat.1005186) and remap these annotations to the newer and updated version 6.0 of the [TriTrypDB](https://tritrypdb.org/) [_L. mexicana_ reference genome](https://tritrypdb.org/common/downloads/release-50/LmexicanaMHOMGT2001U1103/). This remapping of contigs/chromosomes is used in the second step below. To run, please follow these steps:
```bash
cd 1_remap_Fiebig_et_al
sh 1_makeGFF.sh
```
The file `1_remap_Fiebig_et_al/work/genome.gff` contains the genome annotations according to [Fiebig et al. (2015)](https://doi.org/10.1371/journal.ppat.1005186) with the gene model coordinates updated to the newer TriTrypDB reference genome.

## 2_cas9t7_polish

This folder contains code to Pilon-polish the TriTrypDB v50 reference genome with new Illumina sequencing reads from a Cas9/T7-containing laboratory strain. To run, please follow these steps (this assumes the previous remapping step has already been run):
```bash
cd 2_cas9t7_polish
sh 0_getSoftware.sh   # Download and compile (when needed) necessary software components (everything is kept locally, nothing installed globally).
sh 1_polishAndRemap.sh   # Download Illumina sequencing files from sequence read archive and commence Pilon polish run.
```
The files `genome.polish.gff` and `genome.polish.fasta` in `2_cas9t7_polish/work` contain the polished genome annotation and sequences, respectively.

## 3_nanopore_genome_build

This folder contains code for a _de-novo_ _L. mexicana_ Cas9/T7 laboratory strain genome assembly from a Nanopore sequencing run, Pilon-polished with Illumina sequencing reads (the same as in step 2). To run, please follow these steps (independent of the above two steps):
```bash
cd 3_nanopore_genome_build
sh 0_getSoftware.sh   # Download and compile (when needed) necessary software components (everything is kept locally, nothing installed globally).
sh 1_buildGenome.sh   # Download Nanopore and Illumina sequencing files from sequence read archive and commence assembly.
```
The files `assemply.final.fasta` in `3_nanopore_genome_build/work` contains the assembled genome.
