library(seqinr)
source('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/functions.R')
read.fasta('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/data/ecoliGenome.fasta') -> ecoliSeq
View(ecoliSeq)
paste0(ecoliSeq$NC_011750.1)
paste0(ecoliSeq$NC_011750.1,collapse = '') -> ecoliGenome
ecoliGenome <- ecoliGenome[1]
ecoli10k <- substr(ecoliGenome,1,1e4)

reads10k <- generateReads( genome = ecoli10k, readLen = 1000)
reads10kWith5PcentError <- flipReads(reads10k,0.05)

kmers10kWith5pcentError <- sapply(1:length(reads10k),function(x) getKmersFromRead(reads10kWith5PcentError[x]))
kmers10k <- sapply(1:length(reads10k),function(x) getKmersFromRead(reads10k[x]))
kmers10k <- unlist(kmers10k)
beepr::beep(8)
kmers10kWith5pcentError <- unlist(kmers10kWith5pcentError) #Unlist kmers
#get Table
kmers10kE05X30T <- table(kmers10kWith5pcentError) 
#get Df
kmers10kE05X30Df <- as.data.frame(kmers10kE05X30T, stringAsFactors = FALSE )
#remove Table
rm(kmers10kE05X30T)
write.table(kmers10kE05X30Df,'/work-zfs/mschatz1/genomescopeLR/output/kmers10kE05X30Df.txt',sep = '\t',quote = F,row.names = F)