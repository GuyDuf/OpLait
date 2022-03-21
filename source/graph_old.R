
#### Chargement des librairies ####

packages = c("ggplot2","viridis","ShortRead", "stringr", "tidyr", "tibble", "ComplexHeatmap")

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

#Clean the names of the results from the igblast csv file
name_clean_up <- function(string){
  out <- substring(string, first=1 ,last = (str_locate(string = string, pattern = ","))[1] -1)
  if(is.na(out)){
    string
  } else {
    out
  }
}

miseEnOrdre <- function(dataset){
  sort(sapply(id_all[dataset], function(x){
    pos = stringr::str_locate(x,"-Ig")
    stringr::str_sub(x, start = pos[1,"end"] + 1 , end = pos[1,"end"] + 3)
  }))
}

heatmapper <- function(data, VDJ, group){
  IGH.df <- data.frame()
  if(VDJ == "V"){
    for(x in data){
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$v_call, Echantillon = rep(x, nrow(igblast[[x]]))))
    }
  }
  if(VDJ == "D"){
    for(x in echantillons){
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$d_call, Echantillon = rep(x, nrow(igblast[[x]]))))
    }
  }
  if(VDJ == "J"){
    for(x in echantillons){
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$j_call, Echantillon = rep(x, nrow(igblast[[x]]))))
    }
  }
  IGH.df <- data.frame(echantillon = IGH.df$Echantillon, IGH = IGH.df$IGH )
  IGH.tb <- as.data.frame.matrix(table(IGH.df))
  IGH.tb <- IGH.tb[1:ncol(IGH.tb)]/rowSums(IGH.tb[,1:ncol(IGH.tb)])
  #heatmap(as.matrix(IGH.tb),main = paste("IGH",VDJ," : ",group, sep = ""))
  as.matrix(IGH.tb)
}

raceName <- function(ids, rep.by, data){
  to.rep <- str_locate(data[ids],"-[:alnum:]+_S[:digit:]+")
  str_sub(data[ids], to.rep[,1]+1, to.rep[,2]) <- rep.by
  data
}
  

#### Arguments ####
args <- commandArgs(trailingOnly=TRUE)
id_all <-  sort(unique(args))



groupement <- list("G1" = grep(pattern = "IgG1", x = id_all),
  "G2" = grep(pattern = "IgG2", x = id_all),
  "G3" = grep(pattern = "IgG3", x = id_all))
#  "GM1" = grep(pattern = "IgGM1", x = id_all),
#  "GM2" = grep(pattern = "IgGM2", x = id_all))


for(group in names(groupement)){
  
  echantillons <- id_all[groupement[[group]]]
  print(group)
  print(echantillons)
  

  IGHV_possible = c("IGHV1-7","IGHV1-10","IGHV1-14","IGHV1-17","IGHV1-21/33","IGHV1-25","IGHV1-27","IGHV1-30","IGHV1-37","IGHV1-39","IGHV1-20","IGHV1-32")
  IGHD_possible = c("IGHD1-1","IGHD1-2/4","IGHD1-3","IGHD1-4","IGHD2-1/2/3/4","IGHD3-1/3/4","IGHD4-1","IGHD5-2","IGHD5-3/4","IGHD6-2","IGHD6-3/4","IGHD7-3","IGHD7-4","IGHD8-2","IGHD9-1/4")
  IGHJ_possible = c("IGHJ1-4","IGHJ1-6","IGHJ2-4")
  
  #### Files ####
  ## Raw reads
  print(221)
  chemin_raw_reads <- sapply(echantillons, function(x) paste("./data/rawReads/", x,"_L001_R1_001.fastq.gz", sep =""),
                             simplify = FALSE, USE.NAMES = TRUE)
  
  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
                          simplify = FALSE, USE.NAMES = TRUE)
  print(227)
  ## IgBlast
  chemins_igblast <- sapply(echantillons, function(x) paste("./output/VDJ_", x, ".csv2", sep = ""),
                            simplify = FALSE, USE.NAMES = TRUE)
  igblast <- sapply(chemins_igblast, function(x) read.csv(x, sep= "\t", na.strings=c("","NA")),
                    simplify = FALSE, USE.NAMES = TRUE)
  print(233)
  nbr_ligne_igblast_raw <- list(rep(1, times = length(echantillons)))
  
  print("igblast 1/2")
  names(nbr_ligne_igblast_raw) <- names(igblast)
  
  
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
  
  
  pdf(paste("graph/graph", group,".pdf",sep=""),width=15, height=15)
  #### Nombre de sequences total ####
  print(267)
  nbr_seq <- sapply(echantillons, function(x) nrow(igblast[[x]]))
  print(echantillons)
  print(nbr_seq)
  nbr_sequence <- data.frame(Nom = echantillons,
                             Nbr_seq = nbr_seq)
  
 nbr_sequence$nom <-  factor(nbr_sequence$Nom, unique(nbr_sequence$nom))
  print(318)

  print(271)
  plt <- ggplot(data = nbr_sequence, aes(x = Nom, y = nbr_seq)) +
    geom_bar(stat = "identity", color = "black", fill = "darkblue")+
    labs(title = "Valid sequences",
         subtitle = group,
         y = "Number of valid sequences",
         x = "Sample") +
    theme(title =element_text(size=12, face='bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text    = element_text(color = "black", face = "bold", size = 14),
          axis.text.x  = element_text(face = "bold", size = 13)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  print(plt)
  
  x = FALSE
  #nchar(igblast[[i]]$sequence) peut etre remplacé par une colonne longueur de seq
  if(x == TRUE){
    for(j in names(igblast)){
      print(ggplot(igblast[[j]], aes(x = (nchar(igblast[[j]]$sequence)))) +
              geom_histogram(binwidth=3, color = "black", fill = "darkblue") + 
              labs(title= paste("Sequences length (N=", nrow(igblast[[j]]) ,")", sep = ""),
                   subtitle = j,
                   y="Number of sequences",
                   x="Length of senquences (nt)") +
              theme(title = element_text(size=12, face='bold'),
                    axis.title.x = element_text(vjust = 0, size = 15),
                    axis.title.y = element_text(vjust = 2, size = 15),
                    axis.text = element_text(color = "black", face = "bold", size = 14),
                    axis.text.x = element_text(face = "bold", size = 13)))
    }
    
    print("Histogrames des longueurs des séquences : DONE")
  }
  #### Violinplot longueur sequence ####
  print(304)
  longueur <- sapply(names(igblast), function(x) nchar(igblast[[x]]$sequence),
                     simplify = FALSE, USE.NAMES = TRUE)
  
  data_longueur <- data.frame()
  
  temp.df <- sapply(names(longueur), function(x) data.frame(longueur = longueur[[x]],
                                                            nom = rep(x, length(longueur[[x]]))),
                    simplify = FALSE, USE.NAMES = TRUE)

  for(i in (1:length(temp.df))){
    data_longueur <- rbind(data_longueur, temp.df[[i]])
  }
  if(str_sub(group, start= 1 , end = 9)=="duplicate"){
 
  data_longueur$nom <-  factor(data_longueur$nom, unique(data_longueur$nom)[c(1,6,11,2,7,12,3,8,13,4,9,14,5,10,15)])
  print(318)
  catego <- str_sub(data_longueur$nom, start = 1, end = 2)
  
  temp_plot <- ggplot(data_longueur, aes(x = nom, y=longueur, fill = catego))+
    geom_violin(trim = FALSE)+
    geom_boxplot(width=0.1, fill = "lightblue")+
    stat_summary(fun=mean, geom="point", shape=23, size=2)+
    labs(title="Longueur des sequences valides",
         y="Longueur des sequences (nt)",
         x="Echantillons")+
    theme(title =element_text(size=12, face='bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text = element_text(color = "black", face = "bold", size = 14),
          axis.text.x = element_text(face = "bold", size = 13))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
    scale_fill_discrete(name = "Sample", labels = c("Original", "Duplicate", "Time"))
  } else {
    catego <- str_sub(str_extract(data_longueur$nom, "-[:digit:]{3,4}"),2,5)
    temp_plot <- ggplot(data_longueur, aes(x = nom, y=longueur, fill = catego))+
      geom_violin(trim = FALSE)+
      geom_boxplot(width=0.1, fill = "lightblue")+
      stat_summary(fun=mean, geom="point", shape=23, size=2)+
      labs(title="Longueur des sequences valides",
           y="Longueur des sequences (nt)",
           x="Echantillons")+
      theme(title =element_text(size=12, face='bold'),
            axis.title.x = element_text(vjust = 0, size = 15),
            axis.title.y = element_text(vjust = 2, size = 15),
            axis.text = element_text(color = "black", face = "bold", size = 14),
            axis.text.x = element_text(face = "bold", size = 13))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
      scale_fill_discrete(name = "Cow #")
    
}
  print(temp_plot)
  print("Longueur Séquences ViolinPlot : DONE")
  print(338)
  #### Barplot IGHJ NO####
  x = FALSE
  if(x == TRUE){
    for(x in names(igblast)){
      temp_plot <- ggplot(igblast[[x]], aes(x = igblast[[x]]$j_call))+ 
        geom_bar(color = "black", fill = "darkblue") +
        labs(title="Distribution des IGHJ",
             subtitle = x, 
             y="Nombre de sequences",
             x="IGHJ") +
        theme(title =element_text(size=12, face='bold'),
              axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text = element_text(color = "black", face = "bold", size = 14),
              axis.text.x = element_text(face = "bold", size = 13))+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
      print(temp_plot)
    }
    
    print("Barplot IGHJ : DONE")
  }
  
  #### Barplot IGHV NO####
  x = FALSE
  if(x == TRUE){
    for(x in names(igblast)){
      temp_plot <- ggplot(igblast[[x]], aes(x = igblast[[x]]$v_call))+ 
        geom_bar(color = "black", fill = "darkblue") +
        labs(title="Distribution des IGHV",
             subtitle = x,
             y="Nombre de sequences",
             x="IGHV")+
        theme(title =element_text(size=12, face='bold'),
              axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text = element_text(color = "black", face = "bold", size = 14),
              axis.text.x = element_text(face = "bold", size = 13))+
        scale_x_discrete(limits = IGHV_possible)+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
      print(temp_plot)}
    
    print("Barplot IGHV : DONE")
  }
  #### Barplot IGHV CDR3 > 40 NO####
  x = FALSE
  if(x==TRUE){
    for(x in names(igblast)){
      temp_plot <- ggplot(igblast[[x]][nchar(igblast[[x]]$cdr3_aa) >= 40,],
                          aes(x = igblast[[x]][nchar(igblast[[x]]$cdr3_aa) >= 40,]$v_call))+ 
        geom_bar(color = "black", fill = "darkblue") +
        labs(title="Distribution des IGHV avec CDR3 >= 40aa",
             subtitle = x,
             y="Nombre de sequences",
             x="IGHV")+
        theme(title =element_text(size=12, face='bold'),
              axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text = element_text(color = "black", face = "bold", size = 14),
              axis.text.x = element_text(face = "bold", size = 13))+
        scale_x_discrete(limits = IGHV_possible)+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
      print(temp_plot)}
    print("Barplot IGHV CDR3 > 40 : DONE")
    
    
     #### Barplot IGHV CDR3 < 40 NO####
     for(x in names(igblast)){
       temp_plot <- ggplot(igblast[[x]][nchar(igblast[[x]]$cdr3_aa) <= 40,],
                           aes(x = igblast[[x]][nchar(igblast[[x]]$cdr3_aa) <= 40,]$v_call))+ 
         geom_bar(color = "black", fill = "darkblue") +
         labs(title="Distribution des IGHV avec CDR3 <= 40aa",
              subtitle = x,
              y="Nombre de sequences",
              x="IGHV")+
         theme(title =element_text(size=12, face='bold'),
               axis.title.x = element_text(vjust = 0, size = 15),
               axis.title.y = element_text(vjust = 2, size = 15),
               axis.text = element_text(color = "black", face = "bold", size = 14),
               axis.text.x = element_text(face = "bold", size = 13))+
         scale_x_discrete(limits = IGHV_possible)+
         theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
       print(temp_plot)}
     print("Barplot IGHV CDR3 < 40 : DONE")
   }
   #### Heatmap IGHV NO####

  #### Barplot IGHD NO####
  x = FALSE
  if(x==TRUE){
    for(x in names(igblast)){
      temp_plot <- ggplot(igblast[[x]], aes(x = igblast[[x]]$d_call))+ 
        geom_bar(color = "black", fill = "darkblue") +
        labs(title="Distribution des IGHD",
             subtitle = x, 
             y="Nombre de sequences",
             x="IGHD") +
        theme(axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text = element_text(color = "black", face = "bold", size = 12),
              axis.text.x = element_text(face = "bold", size = 11))+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
      print(temp_plot)}
    print("Barplot IGHD : DONE")
  }
  #### HEATMAP IGHD NO ####
  
  
  #### Heatchart IGHV vs IGHD par échantillon NO####
  x = FALSE
  if(x){
  for(x in names(igblast)){
    df <- data.frame(IGHV = igblast[[x]]$v_call, 
                     IGHD = igblast[[x]]$d_call)
    temp <- as.data.frame(table(df))
    temp <- complete(temp, IGHV = IGHV_possible, IGHD = IGHD_possible)
    temp[is.na(temp$Freq),]$Freq <- 0
    
    temp_plot <- ggplot(temp, aes(x = IGHV, y = IGHD, fill = Freq)) +
      geom_tile(aes(fill = Freq), colour = 'black')+
      scale_fill_viridis()+
      labs(title="Distribution de IGHV et IGHD",
           subtitle = x)+
      theme(panel.background = element_rect(fill = 'white'))+
      theme(title =element_text(size=12, face='bold'),
            axis.title.x = element_text(vjust = 0, size = 15),
            axis.title.y = element_text(vjust = 2, size = 15),
            axis.text = element_text(color = "black", face = "bold", size = 13),
            axis.text.x = element_text(face = "bold", size = 13))+
      scale_x_discrete(limits = IGHV_possible)+
      scale_y_discrete(limits = IGHD_possible)+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
    print(temp_plot)}
  print("Heatchart IGHV vs IGHD, par séquence : DONE")
  }
  #### Heatchart IGHV vs IGHD tous NO ####
  x=FALSE
  if(x) {
  allseq <- data.frame()
  for(i in (1:length(igblast))){
    allseq <- rbind(allseq, igblast[[i]])
  }
  #heatmap(table(allseq))
  VS.df <- data.frame(IGHV = allseq$v_call, 
                      IGHD = allseq$d_call)
  
  heatmap(table(VS.df), main = paste("IGHV/IGHD : ", group, sep=""))
  temp <- as.data.frame.matrix(table(VS.df))
  
  temp <- tidyr::complete(temp, IGHV = IGHV_possible, IGHD = IGHD_possible)
  temp[is.na(temp$Freq),]$Freq <- 0
  heatmap(as.matrix(IGH.tb),main = paste("Combinaisons IGHV/IGHD : ",group, sep = ""))

  temp
  
  
  
  IGH.df <- data.frame()
      IGH.df <- rbind(IGH.df,data.frame(IGH = igblast[[x]]$v_call, Echantillon = rep(x, nrow(igblast[[x]]))))

  IGH.df <- data.frame(echantillon = IGH.df$Echantillon, IGH = IGH.df$IGH )
  IGH.tb <- as.data.frame.matrix(table(IGH.df))
  IGH.tb <- IGH.tb[1:ncol(IGH.tb)]/rowSums(IGH.tb[,1:ncol(IGH.tb)])
  
  
  
    temp_plot <- ggplot(temp, aes(x = IGHV, y = IGHD, fill = Freq)) +
    geom_tile(aes(fill = Freq), colour = 'black')+
    scale_fill_viridis()+
    labs(title="Distribution de combinaison de IGHV et IGHD",
         subtitle = "tous")+
    theme(panel.background = element_rect(fill = 'white'))+
    theme(title =element_text(size=12, face='bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text = element_text(color = "black", face = "bold", size = 13),
          axis.text.x = element_text(face = "bold", size = 13))+
    scale_x_discrete(limits = IGHV_possible)+
    scale_y_discrete(limits = IGHD_possible)  +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  
  print(temp_plot)}
  print("Heatchart IGHV vs IGHD, tous : DONE")
  
  
  #### Barplot Longueur CDR3 NO ####
  x = FALSE
  if(x == TRUE){
    for(x in names(igblast)){
      temp_plot <- ggplot(igblast[[x]], aes(x = (nchar(igblast[[x]]$cdr3_aa))))+
        geom_histogram(binwidth=1, color = "black", fill = "darkblue")+
        scale_x_continuous(expand = c(0,1)) + 
        scale_y_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.1)))+
        labs(title="Distribution de la longueur de CDR3",
             subtitle = x,
             y="Nombre de sequences",
             x="Longueur de CDR3")+
        theme(title =element_text(size=12, face='bold'),
              axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text = element_text(color = "black", face = "bold", size = 14),
              axis.text.x = element_text(face = "bold", size = 13))+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
      print(temp_plot)}
    print("Barplot Longueur CDR3 : DONE")
  }
  ##### Boxplot Longueur CDR3 NO ####
  
  longueur_cdr3 <- sapply(names(igblast),
                          function(x) nchar(igblast[[x]]$cdr3_aa),
                          simplify = FALSE, USE.NAMES = TRUE)
  data_longueur_cdr3 <- data.frame()
  temp.df <- sapply(names(longueur_cdr3), 
                    function(x) data.frame(longueur_cdr3 = longueur_cdr3[[x]],
                                           nom = rep(x, length(longueur_cdr3[[x]]))),
                    simplify = FALSE, USE.NAMES = TRUE)
  
  
  for(i in (1:length(temp.df))){
    data_longueur_cdr3 <- rbind(data_longueur_cdr3, temp.df[[i]])
  }
  
  data_longueur_cdr3$nom <-  factor(data_longueur_cdr3$nom, unique(data_longueur_cdr3$nom))
  print(318)
  catego <- str_sub(str_extract(data_longueur_cdr3$nom, "-[:digit:]{3,4}"),2,5)
  
  temp_plot <- ggplot(data_longueur_cdr3, aes(y = longueur_cdr3 ,x = nom, fill = catego))+
    geom_violin()+
    geom_boxplot(width=0.1, fill = "light blue")+
    labs(title="Longueur de CDR3",
         subtitle = group,
         y="Longueur CDR3 (aa)",
         x="Echantillon")+
    theme(title =element_text(size=12, face='bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text = element_text(color = "black", face = "bold", size = 14),
          axis.text.x = element_text(face = "bold", size = 13))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  print(temp_plot)
  print("Boxplot Longueur CDR3 : DONE")
  
  #### NGMerge#######
  
  
  
  ##### Loading Files needed for the end
  x=FALSE
  if(x){
  ## Trimmomatic
  chemins_trimmomatic <- sapply(echantillons, function(x) paste("./data/trimmedReads/", x, "_1P.fastq", sep =""),
                                simplify = FALSE, USE.NAMES = TRUE)
  trimmomatic <- sapply(chemins_trimmomatic, function(x) readFastq(x),
                        simplify = FALSE, USE.NAMES = TRUE)
  
  print("trimmomatic loaded")
  
  ## NGmerge
  chemins_ngmerge <- sapply(echantillons, function(x) paste("./data/mergedReads/", x, ".log", sep =""),
                            simplify = FALSE, USE.NAMES = TRUE)
  ngmerge <-  sapply(chemins_ngmerge, function(x) read.csv(x, sep= "\t", na.strings=c("NA")),
                     simplify = FALSE, USE.NAMES = TRUE)
  
  print("NGmerge loaded")
  
  
  #### ViolinPlot Overlap NGMerge YES ####
  
  
  
  
  
  overlap <- sapply(names(ngmerge), function(x) ngmerge[[x]]$OverlapLen,
                    simplify = FALSE, USE.NAMES = TRUE)
  data_overlap <- data.frame(longueur_overlap = NA, nom = NA)
  temp.df <- sapply(names(longueur),
                    function(x) data.frame(longueur_overlap = overlap[[x]],
                                           nom = rep(x, length(overlap[[x]]))),
                    simplify = FALSE, USE.NAMES = TRUE)
  for(i in (1:length(temp.df))) data_overlap <- rbind(data_overlap, temp.df[[i]])
  
  catego <- str_sub(str_extract(data_overlap$nom, "-+...?"),2,4)
  temp_plot <- ggplot(data = data_overlap[!is.na(data_overlap$nom),], mapping = aes(x = nom, y = longueur_overlap))+
    geom_violin(color = "black", fill = "darkblue")+
    geom_boxplot(width=0.1, fill = "light blue")+
    stat_summary(fun=mean, color = "black", geom="point", shape=23, size=2)+
    labs(title = "Longueur de l'overlap lors de NGMerge",
         y = "Longueur de l'overlap (nt)",
         x = "Echantillon")+
    theme(title = element_text(size = 12, face = 'bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text    = element_text(color = "black", face = "bold", size = 14),
          axis.text.x  = element_text(face = "bold", size = 13))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
    scale_fill_discrete(name = "Cow #")
  print(temp_plot)
  print("Boxplot Overlap NGMerge : DONE")
  }
  dev.off()
  print("dev.off")
  remove(igblast)
  print("IgBlast Removed")
  remove(ngmerge)
  print("Ngmerge removed")
  gc()
  print("gc()")
}
sessionInfo()

