library(seqinr)
source('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/R/functions.R')
read.fasta('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/data/ecoliGenome.fasta') -> ecoliSeq
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
kmers10kE05X30Df <- kmers10kE05X30Df[-which(kmers10kE05X30Df$refDist == Inf),]
kmers10kE05X30Df <- kmers10kE05X30Df %>% arrange(refDist)
kmers10kE05X30Df$index <- seq(1,nrow(kmers10kE05X30Df))
kmers10kE05X30Df$connected <- rep('0',nrow(kmers10kE05X30Df))
print('Finding connected kmers')
dists <- sort(unique(kmers10kE05X30Df$refDist))
sapply(1:(length(dists)-1), function(x){
      km10kE05X30D67 <- kmers10kE05X30Df %>% filter(refDist == dists[x] | refDist == dists[x+1])
      distMatKm10kE05X30D67 <- stringdistmatrix(km10kE05X30D67[,1], method = "hamming")
      dist1 <- which(as.matrix(distMatKm10kE05X30D67) == 1)
      #later add the previous rows but for here since it is the first 2 clusters don't need to add
      dimMat <- nrow(as.matrix(distMatKm10kE05X30D67))
      rowNo <- dist1 %% dimMat
      colNo <- ceiling(dist1/dimMat)
      rc <- data.frame(r = rowNo, c = colNo)
      rc <- rc[rc$r > rc$c,]
      if(x > 1){previous <- which(kmers10kE05X30Df$refDist == dists[x])[1] -1}else{previous <- 0}
      rc <- rc + previous
      sapply(1:nrow(rc),function(x){
            #print(x)
            if(kmers10kE05X30Df$connected[rc$c[x]] == '0'){
                  kmers10kE05X30Df$connected[rc$c[x]] <<- as.character(rc$r[x])
                  #print(kmers10kE05X30Df$connected[rc$c[x]])
            }
            else {
                  kmers10kE05X30Df$connected[rc$c[x]] <<- paste0(kmers10kE05X30Df$connected[rc$c[x]],',',as.character(rc$r[x]))
                  #print(kmers10kE05X30Df$connected[rc$c[x]])
            }
      })
      
      cat(x,'done\n')
})
write.table(kmers10kE05X30Df,'/work-zfs/mschatz1/genomescopeLR/output/kmers10kE05X30Df.txt',sep = '\t',quote = F,row.names = F)