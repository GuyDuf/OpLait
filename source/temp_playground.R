for(j in names(igblast)){
  ggplot(igblast[[j]], aes(x = (nchar(igblast[[j]]$sequence_alignment)))) +
    geom_histogram(binwidth=3, color = "black", fill = "darkblue") + 
    labs(title= paste("Sequences length (N=", nrow(igblast[[j]]) ,")", sep = ""),
         subtitle = j,
         y="Number of sequences",
         x="Length of senquences (nt)") +
    theme(title = element_text(size=12, face='bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text = element_text(color = "black", face = "bold", size = 14),
          axis.text.x = element_text(face = "bold", size = 13))
}

print("Histogrames des longueurs des séquences : DONE")





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

if(all.graph){
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


if(all.graph){
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

if(all.graph) {
  allseq <- data.frame()
  for(i in (1:length(igblast))){
    allseq <- rbind(allseq, igblast[[i]])
  }
  #heatmap(table(allseq))
  VS.df <- data.frame(IGHV = allseq$v_call, 
                      IGHD = allseq$d_call)
  
  
  
  
  
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
  
  print(temp_plot)
  print("Heatchart IGHV vs IGHD, tous : DONE")
  
  
  #### Barplot Longueur CDR3 NO ####
  
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
  
  ##### Boxplot Longueur CDR3 NO ####
  