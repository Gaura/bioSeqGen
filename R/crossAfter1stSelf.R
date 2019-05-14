#Find cross clusters after finding self clusters
library(dplyr)
library(stringdist)
source('/work-zfs/mschatz1/genomescopeLR/bioSeqGen/R/functions.R')
kmer10kE05X30DfRefA <- read.table('/work-zfs/mschatz1/genomescopeLR/output/kmers10kE05X30Df.txt',header = 1,stringsAsFactors = F, sep = '\t')
#Got self-distance reference, now proceed to get cross
getCrossConnection <- function(kmers10kE05X30Df){
      require(dplyr)
      require(stringdist)
      kmers10kE05X30Df$index <- seq(1:nrow(kmers10kE05X30Df))
      dists <- unique(kmers10kE05X30Df$refDist)
      #add refDist, connected and index columns before inputing into function and df sorted according to refDist
      sapply(1:(length(dists)-1), function(x){
            km10kE05X30D6 <- kmers10kE05X30Df %>% filter(refDist == dists[x])
            km10kE05X30D7 <- kmers10kE05X30Df %>% filter(refDist == dists[x+1])
            distMatKm10kE05X30D67 <- stringdistmatrix(a = km10kE05X30D6[,1],b = km10kE05X30D7[,1], method = "hamming")
            distMatKm10kE05X30D67 <- as.matrix(distMatKm10kE05X30D67)
            dist1 <- which(distMatKm10kE05X30D67 == 1)
            #later add the previous rows but for here since it is the first 2 clusters don't need to add
            if(length(dist1) > 0){dimMatR <- nrow(distMatKm10kE05X30D67)
            dimMatC <- ncol(distMatKm10kE05X30D67)
            rowNo <- dist1 %% dimMatR
            if(sum(rowNo == 0) > 0){rowNo[which(rowNo == 0)] <- dimMatR}
            colNo <- ceiling(dist1/dimMatR)
            rc <- data.frame(r = rowNo, c = colNo)
            previousR <- which(kmers10kE05X30Df$refDist == dists[x])[1] -1
            previousC <- which(kmers10kE05X30Df$refDist == dists[x + 1])[1] -1
            rc$r <- rc$r + previousR
            rc$c <- rc$c + previousC
            #print(rc)
            sapply(1:nrow(rc),function(i){
                  if(kmers10kE05X30Df$connected[rc$c[i]] == '0'){
                        kmers10kE05X30Df$connected[rc$c[i]] <<- as.character(rc$r[i])
                  }else {
                        kmers10kE05X30Df$connected[rc$c[i]] <<- paste0(kmers10kE05X30Df$connected[rc$c[i]],',',as.character(rc$r[i]))
                  }
            })
            }
            cat(dists[x],'done\n')
      })
      return(kmers10kE05X30Df)
}

kmers10kE05X30DfCross <- getCrossConnection(kmer10kE05X30DfRefA)
write.table(kmers10kE05X30Df,'/work-zfs/mschatz1/genomescopeLR/output/kmers10kE05X30DfCross.txt',sep = '\t',quote = F,row.names = F)