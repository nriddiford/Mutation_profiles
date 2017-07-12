library(ggplot2)
library(dplyr)

filteR <- function(df=x){
  data<-df
  #filter on chroms
  #data<-filter(data, chrom != "Y" & chrom != "4")
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
    axis.text = element_text(size=12),
    axis.title = element_text(size=30)
  )
}


ts_tv <- function(){
  GW_snvs <- read.table("GW.snv.dist.txt", header = FALSE)
  colnames(GW_snvs)=c("chrom", "bp", "snv", "tri", "trans", "decomposed_tri", "grouped_trans", "sample")
  
  #GW_snvs<-filteR(GW_snvs)
  
  all_ts<-nrow(filter(GW_snvs, trans == "A>G" | trans == "C>T" | trans == "G>A" | trans == "T>C"))
  all_tv<-nrow(filter(GW_snvs, trans != "A>G" & trans != "C>T" & trans != "G>A" & trans != "T>C"))
  ts_tv<-all_ts/all_tv
  cat("ts/tv =", ts_tv)
}

mutation_spectrum <- function(){
  GW_snvs <- read.table("GW.snv.dist.txt", header = FALSE)
  colnames(GW_snvs)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans")
  
  #GW_snvs<-filteR(GW_snvs)
  cat("Showing global contribution of tri class to mutation load", "\n")
 
  p<-ggplot(GW_snvs)
  p<-p + geom_bar(aes(x = decomposed_tri, y = (..count..)/sum(..count..), group = decomposed_tri, fill = grouped_trans), position="dodge",stat="count")
  p<-p + scale_y_continuous("Relative contribution to mutation load", expand = c(0.0, .0005))
  p<-p + scale_x_discrete("Genomic context", expand = c(.005, .005))
  p<-p + clean_theme() + 
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
          )
  p<-p + labs(fill="Mutation class")
  p<-p + facet_wrap(~grouped_trans, ncol = 3, scale = "free_x" )
  
  mut_spectrum<-paste("mutation_spectrum.pdf")
  cat("Writing file", mut_spectrum, "\n")
  ggsave(paste("plots/", mut_spectrum, sep=""), width = 20, height = 10)
  p

}

samples_plot <- function(){
  GW_snvs <- read.table("GW.snv.dist.txt", header = FALSE)
  colnames(GW_snvs)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans")
  
  #GW_snvs<-filteR(GW_snvs)
  
  p<-ggplot(GW_snvs)
  p<-p + geom_bar(aes(x = grouped_trans, y = (..count..)/sum(..count..), group = sample, fill = sample), position="dodge",stat="count")
  p<-p + scale_y_continuous("Relative contribution to total mutation load", expand = c(0.0, .001))
  p<-p + scale_x_discrete("Mutation class")
  p<-p + clean_theme() + 
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 10)
    )
  p<-p + facet_wrap(~sample, ncol = 4, scale = "free_x" )
  
  samples_mut_spect<-paste("mutation_spectrum_samples.pdf")
  cat("Writing file", samples_mut_spect, "\n")
  ggsave(paste("plots/", samples_mut_spect, sep=""), width = 20, height = 10)
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
  
  GW_snps <- read.table("GW.snv.dist.txt", header = FALSE)
  colnames(GW_snps)=c("chrom", "bp", "snv", "tri", "trans", "sample")
  
  p<-ggplot(GW_snps, aes(tri, y = (..count..)/sum(..count..), group = sample))
  p<-p + geom_jitter(aes(colour = sample), size = 1, alpha = 0.9,stat="count")
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
