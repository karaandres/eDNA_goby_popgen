# Goby popgen with eDNA
## This repository contains the files related to the data processing and data analysis for eDNA analysis of round goby microsatellites in North America. 

### Overview
The goal of this project is to assess the capacity for eDNA to be used for population genetic analysis of an invasive species by assessing the similarity between tissue-derived and eDNA-derived metrics of population diversity, differentiation, and structure. We collected 3 2-L eDNA samples at each of 15 sites across the Great Lakes and Finger Lakes regions. We collected round gobies at 13 of those sites. We amplified a panel of 35 microsatellite loci (grouped into 7 multiplexes) for each individual and eDNA sample. Libraries were sequenced across 4 runs of MiSeq 2x250 Paired End sequencing (1 for tissues, 3 for eDNA). Microsatellite multiplexes were barcoded per sample and sequenced together, and samples were parsed into 2 loci during bioinformatic processing (see below). 

### Bioinformatic processing steps
- Create directory on HPC, move all sequence files into directory, clone this repo
- Rename eDNA samples with "e_" to distinguish runs from one another
```
mkdir /workdir/kja68/ && cd /workdir/kja68/
cp /home/kja68/round_goby_field/all_sites_edna/raw_edna_sequences/*.gz ./
cp /home/kja68/round_goby_field/all_sites_edna/raw_tissue_sequences/*.gz ./
git clone https://github.com/karaandres/eDNA_goby_popgen

for f in *J82DH_*; do mv -i -- "$f" "${f//J82DH_/J82DH_e_1_}"; done
for f in *J8RH5_e_*; do mv -i -- "$f" "${f//J8RH5_e_/J8RH5_e_2_}"; done
for f in *JB3WG_e_*; do mv -i -- "$f" "${f//JB3WG_e_/JB3WG_e_3_}"; done
```
**Trimmomatic**: Trim reads that match Nextera adapters or are shorter than 35 bp
- Input is raw reads demultiplexed by sample (.fastq.gz), output is trimmed reads (paired.fastq.gz and unpaired.fastq.gz)
- Make script executable, run on all files
```
cd eDNA_goby_popgen/scripts/
chmod u+x trimmomatic_loop.sh
cd /workdir/kja68/
eDNA_goby_popgen/scripts/trimmomatic_loop.sh
```

**Split_on_Primer**: Split reads into separate files by primer
- Input is paired trimmed reads (paired.fastq.gz), output is separate files for each locus
- Clone script, make executable
```
mkdir split_files/ && cd split_files/
git clone https://github.com/marcomeola/Split_on_Primer.git
chmod u+x Split_on_Primer/src/Split_on_Primer_fixed.py
```
- Copy F and R primer files (format: primer name,sequence) and trimmed sequences into directory 
```
cp /eDNA_goby_popgen/files/keyfile_F_goby_35loci.csv ./
cp /eDNA_goby_popgen/files/keyfile_R_goby_35loci.csv ./
```
- Make directory for split sequence files
```
mkdir sample_files
cp /workdir/kja68/*_paired.fastq.gz ./sample_files/
```
- Unzip files (code only works on fasta or fastq) 
```
gzip -d sample_files/*_paired.fastq.gz
```
- Run script separately for F and R reads; separated files will appear in `split_files/` directory
- Separated files will get an extension on their file name with the locus name, e.g. filename-locusname_F.fastq.gz
- Might take a while for several loci -- use screen so the session stays active
```
screen
for file in sample_files/*R1.fastq.gz_paired.fastq; do Split_on_Primer/src/Split_on_Primer_fixed.py -f "$file" -p keyfile_F_goby_35loci.csv -m 2; done
for file in sample_files/*R2.fastq.gz_paired.fastq; do Split_on_Primer/src/Split_on_Primer_fixed.py -f "$file" -p keyfile_R_goby_35loci.csv -m 2; done
```

**fastq-pair**: Previous script did not provide option to retain matching F and R reads in identical order (required for downstream analyses)
- Input is trimmed files demultiplexed by locus (locusname.fastq.gz), output is matched F and R reads (paired.fq)
- Clone script to match forward and reverse reads and make executable
```
git clone https://github.com/linsalrob/fastq-pair
cd fastq-pair/
mkdir build && cd build
gcc -std=gnu99   ../main.c ../robstr.c ../fastq_pair.c ../is_gzipped.c  -o fastq_pair
```
- Change directory and run script
```
cd ../sample_files/
for f1 in *_F.fastq; do f2a=${f1%%_F.fastq}"_R.fastq"; f2=${f2a/R1.fastq.gz/R2.fastq.gz}; ../fastq-pair/build/fastq_pair -t 1000 $f1 $f2; done
```
- Make directories for each locus (use the same locus names that are appended to your file names) and move split files into the appropriate directory by locus name 
```
mkdir Nmel140 Nmel155 Nmel185 Nmel248 Nmel262 Nmel299 Nmel32 Nmel344 Nmel351 Nmel361 Nmel363 Nmel403 Nmel405 Nmel411 Nmel422 Nmel505 Nmel549 Nmel625 Nmel660 Nmel726 Nmel729 Nmel746 Nmel810 Nmel815 Nmel821 Nmel89 Nmel914 Nmel990 Nmel994 Nmel1103 Nmel1132 Nmel1462 Nmel1486 Nmel1531 Nmel1566
```
- Loops thru the locus directories you just created and grab files by locus name
```
for dir in */; do mv *paired-"${dir%/}"_F.fastq.paired.fq* $dir; done
for dir in */; do mv *paired-"${dir%/}"_R.fastq.paired.fq* $dir; done
```

### Denoise sequences with DADA2 (Callahan et al. 2016): 
- Load packages
```
/programs//R-4.0.5/bin/R
library(ShortRead); packageVersion("ShortRead") # 1.48.0
library(dada2); packageVersion("dada2") # 1.19.2
```
- Run separately for each locus
- Input is F and R sequences demultiplexed by sample and locus and matched in order
- Output is ASV by sample matrix (number of reads per ASV in each sample)
- Some extra steps added because we have lots of samples with no reads (blanks) and DADA2 doesn't like that
```
  path <- paste("/workdir/kja68/split_files/sample_files/", loci[1],sep="")
  setwd(path)
  fnFs <- sort(list.files(path, pattern="F.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="R.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_F.fastq.gz
	sample.names <- sapply(strsplit(basename(fnFs), "-R1"), `[`, 1)

# Inspect read quality profiles -- only for a subset of samples (1:20) and for those with lengths > 0
	fnFs_qual <- NULL
	fnRs_qual <- NULL
	for (i in 1:20){
	srq <- readFastq(fnFs[i])
	seqlen.tab <- table(width(srq))
	if (nrow(seqlen.tab)>0)
	fnFs_qual <- c(fnFs_qual,i)
	srq <- readFastq(fnRs[i])
	seqlen.tab <- table(width(srq))
	if (nrow(seqlen.tab)>0)
	fnRs_qual <- c(fnRs_qual,i)
	}

	pdf("QualityProfile.pdf")
	plotQualityProfile(fnFs[fnFs_qual])
	plotQualityProfile(fnRs[fnRs_qual])
	dev.off()

# Place filtered files in filtered/ subdirectory
	filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names

# Filter and trim
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,290), trimLeft = c(26,26),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE, verbose=TRUE)
	head(out)

# Check how many samples passed filter and remove files that did not pass the filter
	table(file.exists(filtFs))
	table(file.exists(filtRs))
	filtFs <- filtFs[file.exists(filtFs)]
	filtRs <- filtRs[file.exists(filtRs)]

# Learn errors
	errF <- learnErrors(filtFs, multithread=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE)
	pdf("Error_plot.pdf")
	plotErrors(errF, nominalQ=TRUE)
	dev.off()
	dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
	dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
	mergers_20_1 <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=20, maxMismatch=1, verbose=TRUE)

# Construct ASV table 
	seqtab <- makeSequenceTable(mergers_20_1)
	dim(seqtab)

# Inspect distribution of sequence lengths (ASV table)
	table(nchar(getSequences(seqtab)))
	write.csv(table(nchar(getSequences(seqtab))), "Sequence_lengths.csv")

# remove chimeras 
	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	sum(seqtab.nochim)/sum(seqtab)
	dim(seqtab.nochim)

# Save the tables and workspace
	write.csv(t(seqtab), "seqtab.csv")
	write.csv(t(seqtab.nochim), "seqtab.nochim.csv")
	save.image(file='Dada2.RData')

# Save fasta files
	uniquesToFasta(getUniques(seqtab), fout="uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab)))))
	uniquesToFasta(getUniques(seqtab.nochim), fout="uniqueSeqs.nochim.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))
	head(out)
	mergers <- mergers_20_1

# Track reads through the pipeline
	getN <- function(x) sum(getUniques(x))
	track <- cbind(sum(out[,1]), sum(out[,2]), sum(sapply(dadaFs, getN)), sum(sapply(dadaRs, getN)), sum(sapply(mergers, getN)), sum(rowSums(seqtab.nochim)))
	colnames(track) <- c("input", "filtered", "denoisF", "denoisR", "merged", "nonchim")
	track

	pdf("Track_reads.pdf")
	barplot(colSums(track))
	dev.off()
```
