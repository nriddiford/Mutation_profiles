require(dplyr)


trinuc.freqs <- function(genome=NA, write=NA){
  if(is.na(genome)){
    cat("No genome specfied, defaulting to 'BSgenome.Dmelanogaster.UCSC.dm6'\n")
    library(BSgenome.Dmelanogaster.UCSC.dm6, quietly = TRUE)
    genome <- BSgenome.Dmelanogaster.UCSC.dm6
  }
  
  params <- new("BSParams", X = Dmelanogaster, FUN = trinucleotideFrequency, exclude = c("M", "_"), simplify = TRUE)
  data<-as.data.frame(bsapply(params))
  data$genome<-as.integer(rowSums(data))
  data$x<-(data$genome*2)
  gen_wide <- data['x']
  #gen_wide<-cbind(tri = rownames(gen_wide), gen_wide)
  #colnames(gen_wide) <- c("tri", "x")
  #rownames(gen_wide) <- NULL
  
  if(is.na(write)){
    return(gen_wide)
  }
  else{
    cat("Writing genome-wide trinucleotide frequencies to 'data/tri.counts.dmel6.rda'\n")
    write.table(gen_wide, "data/tri.counts.dmel6.txt", sep="\t", quote=FALSE)
    save(gen_wide, file="data/tri.counts.dmel6.rda")
  }
}


possible.trinucs <- function(){
  all.tri = c()
    for(i in c("A", "C", "G", "T")){
      for(j in c("C", "T")){
        for(k in c("A", "C", "G", "T")){
          if(j != k){
            for(l in c("A", "C", "G", "T")){
              tmp = paste(i, "[", j, ">", k, "]", l, sep = "")
              all.tri = c(all.tri, tmp)
            }
          }
        }
      }
    }
  all.tri <- all.tri[order(substr(all.tri, 3, 5))]
  return(all.tri)
}