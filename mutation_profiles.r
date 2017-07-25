library(ggplot2)
suppressMessages(library(dplyr))


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


get_data <- function(infile = "GW.snv.dist.txt"){
  data<-read.delim(infile, header = F)
  
  colnames(data)=c("sample", "chrom", "pos", "ref", "alt", "tri", "trans", "decomposed_tri", "grouped_trans", "type")
  levels(data$type) <- tolower(levels(data$type))
  #data <- filter(data, type == 'germline')
  
  #filter on chroms
  # data<-filter(data, chrom != "Y" & chrom != "4")
  #filter out samples
  # data<-filter(data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  data<-droplevels(data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(data)
}


stats <- function(){
  data<-get_data()
  cat("sample", "snvs", sep='\t', "\n")
  rank<-sort(table(data$sample), decreasing = TRUE)
  rank<-as.array(rank)
  
  for (i in 1:nrow(rank)){
    cat(names(rank[i]), rank[i], sep='\t', "\n")
  }
  
  all_ts<-nrow(filter(data, trans == "A>G" | trans == "C>T" | trans == "G>A" | trans == "T>C"))
  all_tv<-nrow(filter(data, trans != "A>G" & trans != "C>T" & trans != "G>A" & trans != "T>C"))
  ts_tv<-all_ts/all_tv
  cat("ts/tv =", ts_tv)
}


mutation_spectrum <- function(){
  data<-get_data()
  cat("Showing global contribution of tri class to mutation load", "\n")
 
  p<-ggplot(data)
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


samples_plot <- function(count=NA){
  data<-get_data()

  p<-ggplot(data)
  
  if(is.na(count)){
    p<-p + geom_bar(aes(x = grouped_trans, y = (..count..)/sum(..count..), group = sample, fill = sample), position="dodge",stat="count")
    tag='_freq'
  }
  else{
    p<-p + geom_bar(aes(x = grouped_trans, y = ..count.., group = sample, fill = sample), position="dodge",stat="count")
    tag='_count'
  }
  p<-p + scale_y_continuous("Relative contribution to total mutation load", expand = c(0.0, .001))
  p<-p + scale_x_discrete("Mutation class")
  p<-p + clean_theme() + 
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 10)
          )
  p<-p + facet_wrap(~sample, ncol = 4, scale = "free_x" )
  
  samples_mut_spect<-paste("mutation_spectrum_samples", tag, ".pdf", sep = '')
  cat("Writing file", samples_mut_spect, "\n")
  ggsave(paste("plots/", samples_mut_spect, sep=""), width = 20, height = 10)
  p
}


mutational_signatures <- function(samples=NA, pie=NA){
  suppressMessages(require(BSgenome.Dmelanogaster.UCSC.dm6))
  suppressMessages(require(deconstructSigs))
  
  if(!exists('scaling_factor')){
    source('R/dmel6.trinucs.R')
    cat("calculationg trinucleotide frequencies in genome\n")
    scaling_factor <-trinuc.freqs()
  }
  
  data<-get_data()
  genome <- BSgenome.Dmelanogaster.UCSC.dm6
  
  if(is.na(samples)){
    data$tissue = 'All'
    sigs.input <- mut.to.sigs.input(mut.ref = data, sample.id = "tissue", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
    sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.cosmic, sample.id = 'All',
                              contexts.needed = TRUE,
                              tri.counts.method = scaling_factor
                              )

    cat("Writing to file 'plots/all_signatures.pdf'\n")
    pdf('plots/all_signatures.pdf', width = 20, height = 10)
    plotSignatures(sig_plot)
    dev.off()
    plotSignatures(sig_plot)
    

    if(!is.na(pie)){
      makePie(sig_plot)
    }
  }
  
  else{
  	sigs.input <- mut.to.sigs.input(mut.ref = data, sample.id = "sample", chr = "chrom", pos = "pos", alt = "alt", ref = "ref", bsg = genome)
  	cat("sample", "snv_count", sep="\t", "\n")
      for(s in levels(data$sample)) {
        snv_count<-nrow(filter(data, sample == s))
        
        if(snv_count > 50){
          cat(s, snv_count, sep="\t", "\n")
        
          sig_plot<-whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = s,
                        contexts.needed = TRUE,
                        tri.counts.method = scaling_factor)
          
          outfile<-(paste('plots/', s, '_signatures.pdf', sep = ''))
          cat("Writing to file", outfile, "\n")
          pdf(outfile, width = 20, height = 10)
          plotSignatures(sig_plot)
          dev.off()
          plotSignatures(sig_plot)
          
          if(!is.na(pie)){
            makePie(sig_plot)
          }
        }
      }
  }
}


notch_hits <- function(){
  data<-get_data()
  data<-filter(data, chrom == "X", pos >= 3000000, pos <= 3300000)
  
  if(nrow(data) == 0){
    stop("There are no snvs in Notch. Exiting\n")
  }
  
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans, size = 2))
  p<-p + guides(size = FALSE, sample = FALSE)
  p<-p + clean_theme() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(3,3.3,by=0.05), limits=c(3, 3.301))
  p<-p + annotate("rect", xmin=3.000000, xmax=3.134532, ymin=0, ymax=0.1, alpha=.2, fill="green")
  p<-p + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=0, ymax=0.1, alpha=.2, fill="skyblue")
  p<-p + annotate("rect", xmin=3.176440, xmax=3.300000, ymin=0, ymax=0.1, alpha=.2, fill="red")
  
  p
}


genome_wide_snvs <- function(){
  data<-get_data()
  data<-filter(data, chrom != "Y" & chrom != "4")
  p<-ggplot(data)
  p<-p + geom_point(aes(pos/1000000, sample, colour = trans))
  # p<-p + guides(color = FALSE)
  p<-p + theme(axis.text.x = element_text(angle=45, hjust = 1))
  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33), expand = c(0.01, 0.01))
  
  p
}


