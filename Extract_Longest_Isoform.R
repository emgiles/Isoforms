### Script for parsing .faa file to get the longest isoform ###
### Found at https://support.bioconductor.org/p/9140688/ ####

library(AnnotationHub)
library(AnnotationDbi)
library(intervals)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(devtools)
library(bedr)
library(purrr)
source("asynt.R")
#install.packages("splitstackshape")
library(splitstackshape)
#install.packages("rstatix")
library(rstatix)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("multcomp")
library(multcomp)
#install.packages("Biostrings")
library(Biostrings)

getwd()
setwd("/Users/emily/Dropbox/School/Thesis/Genomics-Ch1")

############ Scurria scurra ################
#read in gtf file with pkg GenomicFeatures
Scurra_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/01-SG_genome/assembly_annotation_v1/Scurria_scurra_annotation_v1_ch10_top10.gtf", format="gtf")
saveDb(Scurra_txdb_parsed, "Scurria_scurra_parsed")
S.scurra.p <- loadDb("Scurria_scurra_parsed")
columns(S.scurra.p)

tx_lens <- transcriptLengths(S.scurra.p, with.cds_len=TRUE)
head(tx_lens)

cds_len_per_gene <- splitAsList(tx_lens$cds_len, tx_lens$gene_id)
head(cds_len_per_gene)

tx_id_per_gene <- splitAsList(tx_lens$tx_id, tx_lens$gene_id)
head(tx_id_per_gene)

which_max <- which.max(cds_len_per_gene)
length_of_longest_cds <- unlist(cds_len_per_gene[as.list(which_max)])
tx_id_of_longest_cds <- unlist(tx_id_per_gene[as.list(which_max)])

## 'length_of_longest_cds' is a named integer vector that maps each gene id to the length
## of the longest CDS for that gene:
head(length_of_longest_cds)

## 'tx_id_of_longest_cds' is a named integer vector that maps each gene id to the transcript id
## of the longest CDS for that gene:
head(tx_id_of_longest_cds)

## 'length_of_longest_cds' and 'tx_id_of_longest_cds' are parallel vectors with identical names:
identical(names(length_of_longest_cds), names(tx_id_of_longest_cds))  # TRUE

## Let's summarize these results in a 3-column data.frame with 1 row per gene:
longest_cds <- data.frame(gene_id=names(length_of_longest_cds),
                          length_of_longest_cds=length_of_longest_cds,
                          tx_id_of_longest_cds=tx_id_of_longest_cds)
head(longest_cds)

## Let's remove rows for which the longest CDS has length 0 (non-coding genes):
longest_cds <- subset(longest_cds, length_of_longest_cds != 0)

###Now let's extract all the CDS fragments grouped by transcript and keep only the longest ones:
cds_by_tx <- cdsBy(S.scurra.p, by="tx")
cds_by_tx <- cds_by_tx[as.character(longest_cds$tx_id_of_longest_cds)]
head(cds_by_tx)
####cds_by_tx is a GRangesList object _parallel_ to longest_cds i.e. it has 1 list element per row in longest_cds

####The names on cds_by_tx are transcript ids. You can replace them with gene ids with:
names(cds_by_tx) <- longest_cds$gene_id

############ Scurria viridula ################
#read in gtf file with pkg GenomicFeatures
Viridula_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/02-VG_genome/assembly_annotation_v1/Scurria_viridula_annotation_v1.gff3", format="gtf")
saveDb(Viridula_txdb_parsed, "Scurria_viridula_parsed")
S.viridula.p <- loadDb("Scurria_viridula_parsed")
columns(S.viridula.p)

tx_lens <- transcriptLengths(S.viridula.p, with.cds_len=TRUE)
head(tx_lens)
colnames(tx_lens)

cds_len_per_gene <- splitAsList(tx_lens$cds_len, tx_lens$gene_id)
head(cds_len_per_gene)

tx_id_per_gene <- splitAsList(tx_lens$tx_id, tx_lens$gene_id)
head(tx_id_per_gene)

which_max <- which.max(cds_len_per_gene)
length_of_longest_cds <- unlist(cds_len_per_gene[as.list(which_max)])
tx_id_of_longest_cds <- unlist(tx_id_per_gene[as.list(which_max)])

## 'length_of_longest_cds' is a named integer vector that maps each gene id to the length
## of the longest CDS for that gene:
head(length_of_longest_cds)

## 'tx_id_of_longest_cds' is a named integer vector that maps each gene id to the transcript id
## of the longest CDS for that gene:
head(tx_id_of_longest_cds)

## 'length_of_longest_cds' and 'tx_id_of_longest_cds' are parallel vectors with identical names:
identical(names(length_of_longest_cds), names(tx_id_of_longest_cds))  # TRUE

## Let's summarize these results in a 3-column data.frame with 1 row per gene:
longest_cds <- data.frame(gene_id=names(length_of_longest_cds),
                          length_of_longest_cds=length_of_longest_cds,
                          tx_id_of_longest_cds=tx_id_of_longest_cds)
head(longest_cds)
nrow(longest_cds)

## Let's remove rows for which the longest CDS has length 0 (non-coding genes):
longest_cds <- subset(longest_cds, length_of_longest_cds != 0)
head(longest_cds)

###Now let's extract all the CDS fragments grouped by transcript and keep only the longest ones:
cds_by_tx <- cdsBy(S.viridula.p, by="tx")
cds_by_tx <- cds_by_tx[as.character(longest_cds$tx_id_of_longest_cds)]
head(cds_by_tx)
####cds_by_tx is a GRangesList object _parallel_ to longest_cds i.e. it has 1 list element per row in longest_cds

####The names on cds_by_tx are transcript ids. You can replace them with gene ids with:
names(cds_by_tx) <- longest_cds$gene_id

# Get dna seq
dna <- readDNAStringSet("/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/02-VG_genome/assembly_annotation_v1/Scurria_viridula_assembly_v1.fasta")
head(dna)

# use lapply function and unlist function and diretcly convert into a DNAStringSet
LargeDNAStringSet <- DNAStringSet(lapply(cds_by_tx, function(x) {unlist(x)})) 




