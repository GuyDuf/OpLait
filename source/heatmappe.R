packages = c("ggplot2","viridis","ShortRead", "stringr", "tidyr", "tibble", "ComplexHeatmap","paletteer")


library(paletteer) 


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
name_clean_up <- function(string){
  out <- substring(string, first=1 ,
                   last = (str_locate(string = string, pattern = ","))[1] -1)
  if(is.na(out)){
    string
  } else {
    out
  }
}

heatmapper <- function(data, VDJ, group){
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
  IGH.tb <- IGH.tb[1:ncol(IGH.tb)]/rowSums(IGH.tb[,1:ncol(IGH.tb)])
  #heatmap(as.matrix(IGH.tb),main = paste("IGH",VDJ," : ",group, sep = ""))
  as.matrix(IGH.tb)
}

args <- commandArgs(trailingOnly=TRUE)
#args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")
id_all <-  sort(unique(args))

lookup <- read.csv(file = "./info.csv", row.names = 1)

groupement <- list("G1" = grep(pattern = "IgG1", x = id_all),
                   "G2" = grep(pattern = "IgG2", x = id_all),
                   "G3" = grep(pattern = "IgG3", x = id_all),
                   "GM1" = grep(pattern = "IgGM1", x = id_all),
                   "GM2" = grep(pattern = "IgGM2", x = id_all),
                   "ALL" = grep(pattern = "", x = id_all))


#groupement <- list("ALL" = grep(pattern = "", x = id_all))

for(group in names(groupement)){
  
  echantillons <- id_all[groupement[[group]]]
  print(group)
  print(echantillons)
  
  
  IGHV_possible = c("IGHV1-7","IGHV1-10","IGHV1-14","IGHV1-17","IGHV1-21/33","IGHV1-25","IGHV1-27","IGHV1-30","IGHV1-37","IGHV1-39","IGHV1-20","IGHV1-32")
  IGHD_possible = c("IGHD1-1","IGHD1-2/4","IGHD1-3","IGHD1-4","IGHD2-1/2/3/4","IGHD3-1/3/4","IGHD4-1","IGHD5-2","IGHD5-3/4","IGHD6-2","IGHD6-3/4","IGHD7-3","IGHD7-4","IGHD8-2","IGHD9-1/4")
  IGHJ_possible = c("IGHJ1-4","IGHJ1-6","IGHJ2-4")
 
  chemin_raw_reads <- sapply(echantillons, function(x) paste("./data/rawReads/", x,"_L001_R1_001.fastq.gz", sep =""),
                             simplify = FALSE, USE.NAMES = TRUE)
  
  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
                          simplify = FALSE, USE.NAMES = TRUE)
  print(227)
  ## IgBlast
  chemins_igblast <- sapply(echantillons, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv1", sep = ""),
                            simplify = FALSE, USE.NAMES = TRUE)
  igblast <- sapply(chemins_igblast, function(x) read.csv(x, sep= "\t", na.strings=c("","NA")),
                    simplify = FALSE, USE.NAMES = TRUE)
  
  #for(i in (1:length(igblast))){
  #  igblast[[i]] <- igblast[[i]][,(colnames(igblast[[i]]) %in% to.keep)]
  #  nbr_ligne_igblast_raw[[i]] <- nrow(igblast[[i]])
  #}
  #print("igblast 1/2")
  #names(nbr_ligne_igblast_raw) <- names(igblast)
  
  for(i in (1:length(igblast))){
    igblast[[i]]$v_call <- sapply(igblast[[i]]$v_call, function(x) name_clean_up(x))
    igblast[[i]]$j_call <- sapply(igblast[[i]]$j_call, function(x) name_clean_up(x))
    igblast[[i]]$d_call <- sapply(igblast[[i]]$d_call, function(x) name_clean_up(x))
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$stop_codon == FALSE,]
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$v_frameshift == FALSE,]
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$productive == TRUE,]
    igblast[[i]]        <- igblast[[i]][igblast[[i]]$complete_vdj  == TRUE,]
    igblast[[i]]        <- igblast[[i]][rowSums(is.na(igblast[[i]])) < 4,]
  }
  print("igblast 2/2 loaded")
  
  
  col_fun = viridis(100)
  
  to.hm <- cbind(heatmapper(data = echantillons, VDJ= "V", group = group),
                 heatmapper(data = echantillons, VDJ= "D", group = group),
                 heatmapper(data = echantillons, VDJ= "J", group = group))
  to.hm <- cbind(to.hm, lookup[rownames(to.hm),])

  #>   colorTime <- colorRampPalette(c("blue","deeppink","goldenrod2","green","darkred","blue"))(16)[1:13] 
  #>   colorTime <- colorRampPalette(c("blue","deeppink","goldenrod2","green","darkred","blue"))(16)[1:13] 
  #>   colorTime[14] = "#999999"
  #colorRace = 
  #pie(rep(1,5), col=colorRace)
  #>   colorTime[14] = "#999999"
  #> pie(rep(1,14), col=colorTime)
  
  
  colorRace <- c("Holstein" = "darkgoldenrod2","Ayrshire"= "chartreuse2")
  
  colorTime <- colorRampPalette(c("blue","deeppink","goldenrod2","green","darkred","blue"))(16)[1:13] 
  colorTime[14] = "grey85"
  names(colorTime) <- c("234","241","253","278","279","281","291","603","618","637","641","720","9004","NO")
  colorDuplicate = colorTime[c("234","241","253","603","618","NO")]
  colorClass <- as.vector(plasma(5))
  names(colorClass) <- c("IgG1","IgG2","IgG3","IgGM1","IgGM2")
  
  #Heatmap(heatmapper(data = echantillons, VDJ= "J", group = group))
  if(group == "ALL"){
  HA <-  HeatmapAnnotation(Race = to.hm[,ncol(to.hm)-3],
                           Duplicate = to.hm[,ncol(to.hm)-2],
                           Time = to.hm[,ncol(to.hm)-1],
                           Class = to.hm[,ncol(to.hm)],
                           col = list(Race = colorRace,
                                      Duplicate = colorDuplicate,
                                      Time = colorTime,
                                      Class = colorClass))
  pdf(paste("graph/heatmap", group,".pdf",sep=""),width=20, height=15)
    print(Heatmap(t(to.hm[,1:(ncol(to.hm)-5)]),
                  col = col_fun, column_title = group,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 5),
                  heatmap_legend_param = list(title = "Proportion"),
                  top_annotation = HA))
    
  } else {
  HA <-  HeatmapAnnotation(Race = to.hm[,ncol(to.hm)-3],
                           Duplicate = to.hm[,ncol(to.hm)-2],
                           Time = to.hm[,ncol(to.hm)-1],
                           col = list(Race = colorRace,
                                      Duplicate = colorDuplicate,
                                      Time = colorTime))
  pdf(paste("graph/heatmap", group,".pdf",sep=""),width=15, height=15)
  print(163)
  print(Heatmap(t(to.hm[,1:(ncol(to.hm)-5)]),
                col = col_fun, column_title = group,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = "Proportion"),
                top_annotation = HA))
  print(170)
  }
  
  
  
  print("HEATMAP VDJ : DONE")
  
  
  
  dev.off()
   
}
