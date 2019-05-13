generateReads <- function(genome,depth = 30, readLen = 5000 ){
      set.seed(123)
      genomeLen <- nchar(genome)
      nreads <- genomeLen*depth/readLen
      startPos <- sample(genomeLen,nreads)
      reads <- NULL
      for(i in 1:nreads){
            if(startPos[i] + (readLen - 1) > genomeLen){
                  reads <- c(reads,substr(genome,startPos[i],genomeLen))
            }
            else{
                  reads <- c(reads,substr(genome,startPos[i],(startPos[i] + readLen - 1)))
            }
            
      }
      return(reads)
}

nucleotides <- c('A','T','C','G')

flipReads <- function(rds,error){
      flippedReads <- NULL
      for(i in seq(1,length(rds))){
            subread <- rds[i]
            subreadLen <- nchar(subread)
            nErrors <- floor(subreadLen * error)
            subread <- strsplit(subread,'')[[1]]
            if (nErrors != 0) { subread <- flipBase(subread,subreadLen,nErrors)}
            flippedReads <- c(flippedReads, paste0(subread,collapse = ''))
      }
      return(flippedReads)
}

flipBase <- function(read,readLen,nErr){
      samp <- sample(readLen,nErr)
      for(i in 1:length(samp)){
            basePos <- which(nucleotides == read[samp[i]])
            read[samp[i]] <- sample(nucleotides[-basePos],1)
      }
      return(read)
}

getKmersFromRead <- function(flippedReads,kmerLen = 21){
      readLen <- nchar(flippedReads)
      readKmers <- sapply(1:(readLen - kmerLen +1),function(i){
            substr(flippedReads,i,i+kmerLen-1)
      })
      return(readKmers)
}

selfClusterDist <- function(kmers10kE05X30Df){
require(dplyr)
kmers10kE05X30Df   <- kmers10kE05X30Df %>% arrange(refDist)
kmers10kE05X30Df$index <- seq(1,nrow(kmers10kE05X30Df))
kmers10kE05X30Df$connected <- rep('0',nrow(kmers10kE05X30Df))
print('Finding connected kmers')
dists <- sort(unique(kmers10kE05X30Df$refDist))
kmers10kE05X30Df$connected <- rep('0',nrow(kmers10kE05X30Df))
sapply(1:7, function(x){
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
return(kmers10kE05X30Df)
}