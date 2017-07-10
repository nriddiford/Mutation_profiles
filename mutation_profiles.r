library(ggplot2)
library(dplyr)

filteR <- function(df=x){
  data<-df
  #filter on chroms
  data<-filter(data, chrom != "Y" & chrom != "4")
  #filter out samples
  data<-filter(data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  data<-droplevels(data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(data)
}

clean_theme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=20),
    axis.title = element_text(size=30)
  )
}


global_transitions <- function(stat="identity"){
  if(stat=="identity"){
    cat("Showing relative contribution of tri class to trans")
  }
  
  GW_snps <- read.table("GW.trinucs.txt", header = FALSE)
  colnames(GW_snps)=c("tri", "trans", "freq", "sample")

  p<-ggplot(GW_snps)
  p<-p + geom_bar(aes(x = tri, y = ..count.., group = tri, fill = trans), position="dodge",stat="count")
  p<-p + facet_wrap(~trans, ncol = 2, scale = "free_x" )
  p
}


genome_wide_trinucs <- function(){
  GW_snps <- read.table("GW.trinucs.txt", header = FALSE)
  colnames(GW_snps)=c("tri", "trans", "freq", "sample")
  p<-ggplot(GW_snps, aes(tri,freq, group = sample))
  p<-p + geom_jitter(aes(colour = sample), size = 1, alpha = 0.9)
  p<-p + facet_wrap(~ trans,scale="free_x")
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  p
}


chrom_wide_trinucs <- function(chrom_filt=NA){
  if(is.na(chrom_filt)){
    chrom_filt<-"X"
  }
  by_chrom <- read.table("chroms.trinucs.txt", header = FALSE)
  colnames(by_chrom)=c("chrom", "tri", "trans", "freq", "sample")
  
  chrom_snps <- filter(by_chrom, chrom == chrom_filt)

  p<-ggplot(chrom_snps, aes(tri,freq, group = sample))
  p<-p + geom_jitter(aes(colour = sample), size = 1, alpha = 0.9)
  p<-p + facet_wrap(~ trans,scale="free_x")
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  p<-p + ggtitle( paste( "Chromosome", chrom_filt))
  p
}

genome_wide_snvs <- function(){
  GW_snv_dist <- read.table("GW.snv.dist.txt", header = FALSE)
  colnames(GW_snv_dist)=c("chrom", "bp", "snv", "tri", "trans", "sample")
  
  GW_snv_dist<-filteR(GW_snv_dist)
  
  p<-ggplot(GW_snv_dist)
  p<-p + geom_point(aes(bp/1000000, sample, colour = trans))
  # p<-p + guides(color = FALSE)
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  
  p
}
