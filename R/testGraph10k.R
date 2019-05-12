library(seqinr)
source('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/functions.R')
read.fasta('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/data/ecoliGenome.fasta') -> ecoliSeq
View(ecoliSeq)
paste0(ecoliSeq$NC_011750.1)
paste0(ecoliSeq$NC_011750.1,collapse = '') -> ecoliGenome
ecoliGenome <- ecoliGenome[1]
ecoli10k <- substr(ecoliGenome,1,1e4)