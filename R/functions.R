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

getKmersFromRead <- function(flippedReads,kmerLen){
      readKmers <- NULL
      for(i in 1:length(flippedReads)){
            if(nchar(flippedReads[i]) > kmerLen){
                  kmerEndPos = nchar(flippedReads[i]) - kmerLen + 1
                  for(j in 1:kmerEndPos){
                        readKmers <- c(readKmers,substr(flippedReads[i],j,j+kmerLen-1))
                  }}
            #freqReadKmers <- table(readKmers)
      }
      beepr::beep(8)
      return(readKmers)
}