
# Author : Guyllaume Dufresne
# Description : Generate csv files of length of reads
# Last update : 01/05/22



#### Chargement des librairies ####

packages = c("ShortRead", "stringr", "tidyr", "tibble", "dplyr")

package.check <- lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(x)
    library(x, character.only = TRUE)
  }
}
)
print("Package loaded succesfully!")

#### Functions ####

# Only keeps the first result of igBlast
name_clean_up <- function(string){
  out <- substring(string, first=1 ,last = (str_locate(string = string, pattern = ","))[1] -1)
  if(is.na(out)){
    string
  } else {
    out
  }
}


#### Arguments ####
args <- commandArgs(trailingOnly=TRUE)

# For testing
  #args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")

id_all <-  sort(unique(args))



groupement <- list("G1" = grep(pattern = "IgG1", x = id_all),
                   "G2" = grep(pattern = "IgG2", x = id_all),
                   "G3" = grep(pattern = "IgG3", x = id_all),
                   "GM1" = grep(pattern = "IgGM1", x = id_all),
                   "GM2" = grep(pattern = "IgGM2", x = id_all))



igblast.lst <- list()

for(group in names(groupement)){
  echantillons <- id_all[groupement[[group]]]
  chemins_igblast <- sapply(echantillons, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv2", sep = ""),
                            simplify = FALSE, USE.NAMES = TRUE)
  igblast.lst[[group]] <- sapply(chemins_igblast, function(x) read.csv(x, sep= "\t", na.strings=c("","NA")),
                                 simplify = FALSE, USE.NAMES = TRUE)
}

for(j in names(igblast.lst)){
  for(i in (1:length(igblast.lst[[j]]))){
    print(j)
    print(i)
    igblast.lst[[j]][[i]]$v_call <- sapply(igblast.lst[[j]][[i]]$v_call, function(x) name_clean_up(x))
    igblast.lst[[j]][[i]]$j_call <- sapply(igblast.lst[[j]][[i]]$j_call, function(x) name_clean_up(x))
    igblast.lst[[j]][[i]]$d_call <- sapply(igblast.lst[[j]][[i]]$d_call, function(x) name_clean_up(x))
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][!is.na(igblast.lst[[j]][[i]]$stop_codon),]
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][igblast.lst[[j]][[i]]$stop_codon == FALSE,]
    
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][!is.na(igblast.lst[[j]][[i]]$v_frameshift),]
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][igblast.lst[[j]][[i]]$v_frameshift == FALSE,]
    
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][!is.na(igblast.lst[[j]][[i]]$vj_in_frame),]
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][igblast.lst[[j]][[i]]$vj_in_frame == TRUE,]
    
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][!is.na(igblast.lst[[j]][[i]]$complete),]
    igblast.lst[[j]][[i]] <- igblast.lst[[j]][[i]][igblast.lst[[j]][[i]]$complete_vdj  == TRUE,]
    
    igblast.lst[[j]][[i]]$sequence_alignment <- nchar(igblast.lst[[j]][[i]]$sequence_alignment)
    igblast.lst[[j]][[i]]$cdr3_aa <- nchar(igblast.lst[[j]][[i]]$cdr3_aa)
    
    # I want to keep empty IGHD
    #igblast.lst[[j]][[i]] <- na.omit(igblast.lst[[j]][[i]])
    
  }
}


for(group in names(groupement)){
  print(group)
  for(i in names(igblast.lst[[group]])){
    write.csv(igblast.lst[[group]][[i]][,-1], paste("csv/",i,"_length.csv",sep = ""))
    }
}

write.table(x = c("delete me only to clear the folder"), file = "csv/done.txt")
