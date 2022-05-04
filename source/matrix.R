# Author : Guyllaume Dufresne
# Description : Generate heatmap of the composition of reads
# Last update : 03/05/22

packages = c("ggplot2","viridis","ShortRead", "stringr", "tidyr", "tibble", "ComplexHeatmap","paletteer")

package.check <- lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(x)
    library(x, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
  }
}
)
print("Package loaded")


# Only keeps the first result of igBlast
name_clean_up <- function(string){
  out <- substring(string, first=1 ,
                   last = (str_locate(string = string, pattern = ","))[1] -1)
  if(is.na(out)){
    string
  } else {
    out
  }
}


# Generate the matrix for the heatmap
heatmapper <- function(data, VDJ, group, prop = TRUE){
  IGH.df <- data.frame()
  if(VDJ == "V"){
    for(x in data){
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$v_call,
                                        Echantillon = rep(x, nrow(igblast[[x]]))))
    }
  }
  if(VDJ == "D"){
    for(x in echantillons){
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$d_call,
                                        Echantillon = rep(x, nrow(igblast[[x]]))))
    }
  }
  if(VDJ == "J"){
    for(x in echantillons){
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$j_call,
                                        Echantillon = rep(x, nrow(igblast[[x]]))))
    }
  }
  IGH.df <- data.frame(echantillon = IGH.df$Echantillon, IGH = IGH.df$IGH )
  
  IGH.tb <- as.data.frame.matrix(table(IGH.df))
  
  if(prop){
    IGH.tb <- IGH.tb[1:ncol(IGH.tb)]/rowSums(IGH.tb[,1:ncol(IGH.tb)]) 
  } else {
    IGH.tb <- IGH.tb[1:ncol(IGH.tb)]
  }
  
  as.matrix(IGH.tb)
}



#Gets the names of the reads to use
args <- commandArgs(trailingOnly=TRUE)

# For testing
  args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","10-618-IgG3-1_S10","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")

id_all <-  sort(unique(args))
lookup <- read.csv(file = "./data/info.csv", row.names = 1)


# Create the groups of reads
groupement <- list("G1" = grep(pattern = "IgG1", x = id_all),
                   "G2" = grep(pattern = "IgG2", x = id_all),
                   "G3" = grep(pattern = "IgG3", x = id_all),
                   "GM1" = grep(pattern = "IgGM1", x = id_all),
                   "GM2" = grep(pattern = "IgGM2", x = id_all),
                   "All" = grep(pattern = "", x = id_all))

IGHV_possible <-  c("IGHV1-7","IGHV1-10","IGHV1-14","IGHV1-17","IGHV1-21/33","IGHV1-25","IGHV1-27","IGHV1-30","IGHV1-37","IGHV1-39","IGHV1-20","IGHV1-32")
IGHD_possible <-  c("IGHD1-1","IGHD1-2/4","IGHD1-3","IGHD2-1/2/3/4","IGHD3-1/3/4","IGHD4-1","IGHD5-2","IGHD5-3/4","IGHD6-2","IGHD6-3/4","IGHD7-3","IGHD7-4","IGHD8-2","IGHD9-1/4")
IGHJ_possible <-  c("IGHJ1-4","IGHJ1-6","IGHJ2-4")

for(group in names(groupement)){
  
  
  # Select the group of reads to use
  echantillons <- id_all[groupement[[group]]]
  
  
  # read the files
  chemin_raw_reads <- sapply(echantillons, function(x) paste("./data/rawReads/", x,"_L001_R1_001.fastq.gz", sep =""),
                             simplify = FALSE, USE.NAMES = TRUE)
  
  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
                          simplify = FALSE, USE.NAMES = TRUE)
  
  chemins_igblast <- sapply(echantillons, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv1", sep = ""),
                            simplify = FALSE, USE.NAMES = TRUE)
  
  igblast <- sapply(chemins_igblast, function(x) read.csv(x, sep= "\t", na.strings=c("","NA")),
                    simplify = FALSE, USE.NAMES = TRUE)
  
  print("Files Loaded")
  
  # Clean the names of region and remove unwanted reads
  for(i in (1:length(igblast))){
    igblast[[i]]$v_call <- sapply(igblast[[i]]$v_call, function(x) name_clean_up(x))
    igblast[[i]]$j_call <- sapply(igblast[[i]]$j_call, function(x) name_clean_up(x))
    igblast[[i]]$d_call <- sapply(igblast[[i]]$d_call, function(x) name_clean_up(x))
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$stop_codon == FALSE,]
    igblast[[i]]        <- igblast[[i]][!is.na(igblast[[i]]$stop_codon),]
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$v_frameshift == FALSE,]
    igblast[[i]]        <- igblast[[i]][!is.na(igblast[[i]]$v_frameshift),]
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$productive == TRUE,]
    igblast[[i]]        <- igblast[[i]][!is.na(igblast[[i]]$productive),]
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$complete_vdj  == TRUE,]
    igblast[[i]]        <- igblast[[i]][!is.na(igblast[[i]]$complete_vdj),]
    igblast[[i]]        <- igblast[[i]][rowSums(is.na(igblast[[i]])) < 4,]
    
    # I want to keep empty IGHD
    #  igblast[[i]] <- na.omit(igblast[[i]])
  }
  print("Removed unwanted reads")
  
  to.hm <- cbind(heatmapper(data = echantillons, VDJ= "V", group = group, prop = TRUE),
       heatmapper(data = echantillons, VDJ= "D", group = group, prop = TRUE),
       heatmapper(data = echantillons, VDJ= "J", group = group, prop = TRUE))
  to.hm <- cbind(to.hm, lookup[rownames(to.hm),])
  saveRDS(to.hm, file = paste("out/distribution",group,"_prop.rds",sep=""))
  
  to.hm <- cbind(heatmapper(data = echantillons, VDJ= "V", group = group, prop = FALSE),
                 heatmapper(data = echantillons, VDJ= "D", group = group, prop = FALSE),
                 heatmapper(data = echantillons, VDJ= "J", group = group, prop = FALSE))
  to.hm <- cbind(to.hm, lookup[rownames(to.hm),])
  saveRDS(to.hm, file = paste("out/distribution",group,"_nbr.rds",sep=""))
}
