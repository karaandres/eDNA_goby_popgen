# Goby popgen with eDNA
## This repository contains the files related to the data processing and data analysis for eDNA analysis of round goby microsatellites in North America. 

### Overview
The goal of this project is to assess the capacity for eDNA to be used for population genetic analysis of an invasive species by assessing the similarity between tissue-derived and eDNA-derived metrics of population diversity, differentiation, and structure. We collected 3 2-L eDNA samples at each of 15 sites across the Great Lakes and Finger Lakes regions. We collected round gobies at 13 of those sites. We amplified a panel of 35 microsatellite loci (grouped into 7 multiplexes) for each individual and eDNA sample. Libraries were sequenced across 4 runs of MiSeq 2x250 Paired End sequencing (1 for tissues, 3 for eDNA). Microsatellite multiplexes were barcoded per sample and sequenced together, and samples were parsed into 2 loci during bioinformatic processing (see below). 

### Nuclear eDNA vs. mitochondrial eDNA (qPCR results)
- Input: [2L_COI_msat_qPCR_results.csv](2L_COI_msat_qPCR_results.csv)
- Analysis: [qPCR_results_COI_msat.R](scripts/qPCR_results_COI_msat.R)

### Bioinformatic processing steps
- Sequence files available on NCBI PRJNA997584
- Rename eDNA samples with "e_" to distinguish runs from one another
```
for f in *J82DH_*; do mv -i -- "$f" "${f//J82DH_/J82DH_e_1_}"; done
for f in *J8RH5_e_*; do mv -i -- "$f" "${f//J8RH5_e_/J8RH5_e_2_}"; done
for f in *JB3WG_e_*; do mv -i -- "$f" "${f//JB3WG_e_/JB3WG_e_3_}"; done
```

**Trimmomatic**: [trimmomatic_loop.sh](scripts/trimmomatic_loop.sh) Trim reads that match Nextera adapters or are shorter than 35 bp
- Input is raw reads demultiplexed by sample (.fastq.gz), output is trimmed reads (paired.fastq.gz and unpaired.fastq.gz)
- Make script executable, run on all files

```
chmod u+x trimmomatic_loop.sh
./trimmomatic_loop.sh
```

**Denoising**: 
- Run [amplicon_dada2.py](https://bitbucket.org/cornell_bioinformatics/amplicon/src/master/amplicon_dada2.py): This script enables DADA2 processing of sequences from multiple loci at once with an added Smith-Waterman alignment of ASVs to the most common allele at that locus
- The following commands should be installed and in the PATH: [bbmerge.sh](https://sourceforge.net/projects/bbmap/), [cutadapt](https://cutadapt.readthedocs.io/en/stable/) (v3 or above), python package [swalign](https://pypi.org/project/swalign/), R package [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html)
- Generate reference files: [keyfile](keyfile_FR_goby_35loci_7july2021.txt) and [sample file](sample_file_goby_field_28oct2021.txt)
- Input: trimmed reads and reference files

```
amplicon_dada2.py -s sample_file_goby_field_28oct2021.txt -k keyfile_FR_goby_35loci_7july2021.txt -o outputDir -j 24 -r 1 -f 2
```
- Output is sequence reads per allele per locus per site, as well as a fasta file of each allele per locus

### Calculate allele frequencies for tissues and eDNA: [calculate_allele_frequencies.R](scripts/calculate_allele_frequencies.R)
- Inputs: Sequence reads per allele per locus per site (amplicon_dada2.py output)

### Tissue-based and eDNA-based populations genetics: [eDNA_tissue_popgen_analyses.R](scripts/eDNA_tissue_popgen_analyses.R), [structure](structure)
- Inputs: [hap_genotype_matrix.csv](hap_genotype_matrix.csv), [tissue_allele_freqs.csv](tissue_allele_freqs.csv), [edna_allele_freqs.csv](edna_allele_freqs.csv), [goby_strata.csv](goby_strata.csv)

