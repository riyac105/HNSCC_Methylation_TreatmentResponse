library(tidyverse)
library(Gviz)

#Example analysis for EGFR

#STEP 1: Read in methylation M values and protein expression data
#Input file was preprocessed and organized using TCGA-HNSC data
#Illumina methylation data was obtained from GDC
#Protein expression data was obtained from UCSC Xena: https://xenabrowser.net/datapages/?dataset=TCGA-HNSC.protein.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443, total samples=354
#Protein expression was log2 transformed for statistical stability
#Our processed input file has been provided as supplementary data
EGFR_correlation_data <- read.csv("EGFR_methylation_and_expression.csv")

#STEP 2: Run Spearman's Rank correlation for each CpG site vs log transformed protein expression data
#Prepare dataframe for results
spearman_results_EGFR <- data.frame(
  Probe = character(),
  Spearman_Coefficient = numeric(),
  p_value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE 
)

for (i in 3:55) { #Choose thresholds corresponding to CpG site methylation data
  for (j in 57:57) { #Choose column corresponding to log transformed expression
    spr <- cor.test(as.numeric(EGFR_correlation_data[, i]), as.numeric(EGFR_correlation_data[, j]), method = "spearman", use = "complete.obs") #Only use rows without N/A values
    if(is.nan(spr$p.value)){ 
      sig <- "FALSE"
    } else if (spr$p.value <= 0.05) { #Assigning statistical signifiance based on significance level (α) 0.05
      sig <- "TRUE"
    } else {
      sig <- "FALSE"
    }
    
    spearman_results_EGFR <- rbind(spearman_results_EGFR, data.frame(
      Probe = names(EGFR_correlation_data)[i],
      Spearman_Coefficient = spr$estimate,
      P_value = spr$p.value,
      Significance = sig,
      stringsAsFactors = FALSE
    ))
  }
}

#STEP 3: Plot resulting coefficient values based on genomic coordinates
#Input file was filtered through the 
#Genome Coordinate obtained from: https://xenabrowser.net/datapages/?dataset=TCGA-HNSC.methylation450.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#File download available from this link: https://api.gdc.cancer.gov/v0/data/021a2330-951d-474f-af24-1acd77e7664f
gene_coord <- read.csv("EGFR_GENCODE.csv")
gene_coord <- gene_coord[,-1] #Removes unnecessary column in dataframe

map_data <- inner_join(spearman_results_EGFR, gene_coord, by="Probe")
map_data <- map_data[,-10]

#The Gviz package was used to generate these plots. The vignette used to guide our script can be found here: https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html
#Create gRanges object for gene mape
gRanges_EGFR <- GRanges( 
  seqnames = map_data$CpG_chrm,
  ranges = IRanges(start = map_data$CpG_beg,
                   end = map_data$CpG_end),
  CpG_ID = map_data$Probe,
  gene_id = map_data$genesUniq
)
atrack_EGFR <- AnnotationTrack(gRanges_EGFR, name = "CpG sites") #Can be included as a parameter in plotTracks to display CpG site locations
gtrack <- GenomeAxisTrack() #shows location
dtrack_EGFR <- DataTrack(data = map_data$Spearman_Coefficient, start = map_data$CpG_beg,
                           end = map_data$CpG_end, chromosome = map_data$CpG_chrm, genome = map_data$genesUniq, 
                           name = "Spearman Coefficient Values")
plotTracks(list(gtrack, dtrack_EGFR), type="histogram")
