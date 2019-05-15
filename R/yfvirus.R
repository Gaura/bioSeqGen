library(seqinr)
library(dplyr)
source('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/R/functions.R')
read.fasta('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/data/yellow_fever_virus.fasta') -> ecoliSeq
print('Loading Genome')
paste0(ecoliSeq$NC_011750.1,collapse = '') -> ecoliGenome
ecoliGenome <- ecoliGenome[1]
ecoli10k <- substr(ecoliGenome,1,1e4)
ecoli10k <- toupper(ecoli10k)
reads10k <- generateReads( genome = ecoli10k, readLen = 1000)
reads10kWith5PcentError <- flipReads(reads10k,0.05)
print('Getting reads and kmers')
kmers10kWith5pcentError <- sapply(1:length(reads10k),function(x) getKmersFromRead(reads10kWith5PcentError[x]))
kmers10k <- sapply(1:length(reads10k),function(x) getKmersFromRead(reads10k[x]))
kmers10k <- unlist(kmers10k)
kmers10kWith5pcentError <- unlist(kmers10kWith5pcentError) #Unlist kmers
#get Table
kmers10kE05X30T <- table(kmers10kWith5pcentError) 
#get Df
kmers10kE05X30Df <- as.data.frame(kmers10kE05X30T, stringAsFactors = FALSE )
#remove Table
rm(kmers10kE05X30T)
library(stringdist)
print('Calculating Distance')
refKmer <- paste0(rep('A',21),collapse = '')
refDist10kE05X30 <- sapply(1:nrow(kmers10kE05X30Df),function(x){stringdist(kmers10kE05X30Df[x,1],refKmer,method = "hamming")})
kmers10kE05X30Df$refDist <- refDist10kE05X30
#remove Inf distance, later add check to see if it exists
if(Inf %in% kmers10kE05X30Df$refDist){
      kmers10kE05X30Df <- kmers10kE05X30Df[-which(kmers10kE05X30Df$refDist == Inf),]}
kmers10kE05X30Df <- kmers10kE05X30Df %>% arrange(refDist)
kmers10kE05X30Df$index <- seq(1,nrow(kmers10kE05X30Df))
kmers10kE05X30Df$connected <- rep('0',nrow(kmers10kE05X30Df))
print('Finding connected kmers')
dists <- sort(unique(kmers10kE05X30Df$refDist))
kmers10kE05X30Df$connected <- rep('0',nrow(kmers10kE05X30Df))
sapply(1:length(dists), function(x){
      km10kE05X30D67 <- kmers10kE05X30Df %>% filter(refDist == dists[x])
      distMatKm10kE05X30D67 <- stringdistmatrix(km10kE05X30D67[,1], method = "hamming")
      distMatKm10kE05X30D67 <- as.matrix(distMatKm10kE05X30D67)
      dist1 <- which(distMatKm10kE05X30D67 == 1)
      if(length(dist1) != 0){
            #print(dist1)
            #later add the previous rows but for here since it is the first 2 clusters don't need to add
            dimMat <- nrow(distMatKm10kE05X30D67)
            rowNo <- dist1 %% dimMat
            colNo <- ceiling(dist1/dimMat)
            rc <- data.frame(r = rowNo, c = colNo)
            rc <- rc[rc$r > rc$c,]
            previous <- which(kmers10kE05X30Df$refDist == dists[x])[1] -1
            rc <- rc + previous
            #print(rc)
            sapply(1:nrow(rc),function(i){
                  #print(x)
                  if (kmers10kE05X30Df$connected[rc$c[i]] == '0') {
                        kmers10kE05X30Df$connected[rc$c[i]] <<- as.character(rc$r[i])
                        #print(kmers10kE05X30Df$connected[rc$c[x]])
                  }
                  else {
                        kmers10kE05X30Df$connected[rc$c[i]] <<- paste0(kmers10kE05X30Df$connected[rc$c[i]],',',as.character(rc$r[i]))
                        #print(kmers10kE05X30Df$connected[rc$c[x]])
                  }
            })
      }
      
      cat(x,'done\n')
})
write.table(kmers10kE05X30Df,'/work-zfs/mschatz1/genomescopeLR/output/YFVkmers10kE05X30Df.txt',sep = '\t',quote = F,row.names = F)
kmers10kE05X30DfCross <- getCrossConnection(kmers10kE05X30Df)
write.table(kmers10kE05X30DfCross,'/work-zfs/mschatz1/genomescopeLR/output/YFVkmers10kE05X30DfCross.txt',sep = '\t',quote = F,row.names = F)