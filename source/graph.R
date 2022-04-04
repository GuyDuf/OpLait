
#### Chargement des librairies ####

packages = c("ggplot2","viridis","ShortRead", "stringr", "tidyr", "tibble", "dplyr")

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


#raceName <- function(ids, rep.by, data){
#  to.rep <- str_locate(data[ids],"-[:alnum:]+_S[:digit:]+")
#  str_sub(data[ids], to.rep[,1]+1, to.rep[,2]) <- rep.by
#  data
#}
  

#### Arguments ####
args <- commandArgs(trailingOnly=TRUE)
args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")

id_all <-  sort(unique(args))

all.graph <- FALSE



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
  
total.reads.lst <- list()
for(group in names(groupement)){
  chemin_raw_reads <- sapply(id_all[groupement[[group]]], function(x) paste("./data/rawReads/", x,"_L001_R1_001.fastq.gz", sep =""),
                             simplify = FALSE, USE.NAMES = TRUE)
  
  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
                          simplify = FALSE, USE.NAMES = TRUE)  
  total.reads.lst[group] <- sum(unlist(nbr_raw_reads))
}
  
total.reads.df1 <- do.call(cbind.data.frame, total.reads.lst)
total.reads.df <- data.frame(nbr.reads = t(total.reads.df1)[,1], class = c("G1","G2","G3","GM1","GM2"))
  
  IGHV_possible = c("IGHV1-7","IGHV1-10","IGHV1-14","IGHV1-17","IGHV1-21/33","IGHV1-25","IGHV1-27","IGHV1-30","IGHV1-37","IGHV1-39","IGHV1-20","IGHV1-32")
  IGHD_possible = c("IGHD1-1","IGHD1-2/4","IGHD1-3","IGHD1-4","IGHD2-1/2/3/4","IGHD3-1/3/4","IGHD4-1","IGHD5-2","IGHD5-3/4","IGHD6-2","IGHD6-3/4","IGHD7-3","IGHD7-4","IGHD8-2","IGHD9-1/4")
  IGHJ_possible = c("IGHJ1-4","IGHJ1-6","IGHJ2-4")
  #safe = igblast.lst
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
      
      igblast.lst[[j]][[i]] <- na.omit(igblast.lst[[j]][[i]])
      
    }
  }
  
  ####Reason to Drop####

  #### Violinplot longueur sequence ####
  longueur.lst <- list()
  for(group in names(igblast.lst)){
    temp <- igblast.lst[[group]]
    for(i in 1:length(temp)){
      temp[i]
     longueur.lst[[group]] <- append(longueur.lst[[group]], sapply(temp[i], function(x) nchar(x$sequence_alignment),
                     simplify = TRUE, USE.NAMES = TRUE))
      }
  }
  
  toute <- list()
  for(i in 1:length(longueur.lst)){
  
    toute[[i]] <- data.frame(longueur = longueur.lst[[i]], class = names(longueur.lst[i]))
  }
  
  longueur.df <- do.call(rbind.data.frame, toute)
  
  nbr.row <- list() 
  for(group in names(groupement)){
    
      temp <- na.omit(igblast.lst[[group]][[1]])
      nbr.row[[group]] <- nrow(temp)
    for(i in 2:length(igblast.lst[[group]])){
      temp <- na.omit(igblast.lst[[group]][[i]])
      nbr.row[[group]] <- nbr.row[[group]] + nrow(temp)
    }
  }
  
  
  percent <- formatC(unlist(nbr.row / total.reads.df1) * 100, digits = 2, format ="f")

  options(scipen = 1000000)
  axis.names =c(paste("G1 \n n = " ,nbr.row[["G1"]],"\n", formatC(unlist(percent[1]), digits = 2, format = "f"),"% of ", total.reads.df1$G1 ," reads", sep = ""),
                paste("G2 \n n = " ,nbr.row[["G2"]],"\n", formatC(unlist(percent[2]), digits = 2, format = "f"),"% of ", total.reads.df1$G2," reads", sep = ""),
                paste("G3 \n n = " ,nbr.row[["G3"]],"\n", formatC(unlist(percent[3]), digits = 2, format = "f"),"% of ", total.reads.df1$G3," reads", sep = ""),
                paste("GM1 \n n = ",nbr.row[["GM1"]],"\n", formatC(unlist(percent[4]), digits = 2, format = "f"),"% of ", total.reads.df1$GM1 ," reads", sep = ""),
                paste("GM2 \n n = ",nbr.row[["GM2"]],"\n", formatC(unlist(percent[5]), digits = 2, format = "f"),"% of ", total.reads.df1$GM2 ," reads", sep = ""))
  
  
  pdf("graph/length.pdf",width=15, height=15)
  
  plt <- ggplot(data = longueur.df, aes(x= class, y = longueur))+
         geom_violin(trim = FALSE, color = "black", bw = 1.75)+
         geom_boxplot(width=0.1, fill = "lightblue")+
         stat_summary(fun=mean, geom="point", shape=23, size=2)+
         labs(title="Valid sequences length",
              y="Sequence length (nt)",
              x="Class")+
         theme(title =element_text(size=12, face='bold'),
               axis.title.x = element_text(vjust = 0, size = 15),
               axis.title.y = element_text(vjust = 2, size = 15),
               axis.text = element_text(color = "black", face = "bold", size = 14),
               axis.text.x = element_text(face = "bold", size = 13))+
         scale_x_discrete(labels = axis.names)

  print(plt)
  ### CDR3 GRAPHS ####
  
  longueur.cdr3.lst <- list()
  for(group in names(igblast.lst)){
    temp <- igblast.lst[[group]]
    for(i in 1:length(temp)){
      temp[i]
      longueur.cdr3.lst[[group]] <- append(longueur.cdr3.lst[[group]], 
                                           sapply(temp[i], function(x) nchar(x$cdr3_aa),
                                           simplify = TRUE, USE.NAMES = TRUE))
    }
  }
  for(group in names(igblast.lst)){
  for(i in 1:length(igblast.lst[[group]])){
    print(group)
    print(i)
    print(cor(nchar(igblast.lst[[group]][[i]]$cdr3_aa), nchar(igblast.lst[[group]][[i]]$sequence_alignment)))}
  }
  
  toute.cdr3 <- list()
  for(i in 1:length(longueur.cdr3.lst)){
    toute.cdr3[[i]] <- data.frame(longueur = longueur.cdr3.lst[[i]], class = names(longueur.cdr3.lst[i]))
  }
  names(toute.cdr3) <- c("G1","G2","G3","GM1","GM2")
  longueur.cdr3.df <- do.call(rbind.data.frame, toute.cdr3)
  
  
  plt <- ggplot(data = longueur.cdr3.df, aes(x= class, y = longueur))+
         geom_violin(trim = FALSE, color = "black", bw = 0.75)+
         geom_boxplot(width=0.1, fill = "lightblue")+
         stat_summary(fun=mean, geom="point", shape=23, size=2)+
         labs(title="Length of CDR3 region",
              y="Length of region (AA)",
              x="Class")+
         theme(title =element_text(size=12, face='bold'),
               axis.title.x = element_text(vjust = 0, size = 15),
               axis.title.y = element_text(vjust = 2, size = 15),
               axis.text = element_text(color = "black", face = "bold", size = 14),
               axis.text.x = element_text(face = "bold", size = 13))+
         scale_x_discrete(labels = axis.names)
         
  print(plt)
  
  dev.off()
  print("ViolinPlot CDR3 : DONE")
         
    
    
    less.40 <- list(); more.40 <- list(); less.40.v <- list(); more.40.v <- list(); moss.40.v <- list()
    for(group in names(igblast.lst)){
      pdf(paste("graph/moss_",group,".pdf",sep=""),15,15)
      #nbr.row <- 0
      for(i in 1:length(igblast.lst[[group]])){
        temp <- na.omit(igblast.lst[[group]][[i]])
        less.40[[group]] <- rbind(as.data.frame(less.40[[group]]),na.omit(temp[nchar(temp$cdr3_aa)<40,]))
        more.40[[group]] <- rbind(as.data.frame(more.40[[group]]),na.omit(temp[nchar(temp$cdr3_aa)>=40,]))
      }
        more.40[[group]] %>% dplyr::count(v_call) %>% mutate(perc = (n / nbr.row[[group]]*100)) -> more.40.v[[group]]
        less.40[[group]] %>% dplyr::count(v_call) %>% mutate(perc = (n / nbr.row[[group]]*100)) -> less.40.v[[group]]
        
        for(V in IGHV_possible){
          if(!(V %in% more.40.v[[group]]$v_call)){
            more.40.v[[group]] <- rbind(more.40.v[[group]],c(v_call = V, n=0, perc=0))
          }
          if(!(V %in% less.40.v[[group]]$v_call)){
            less.40.v[[group]] <- rbind(less.40.v[[group]],c(v_call = V, n=0, perc=0))
          }
        }
        
        more.40.v[[group]]$perc <- as.numeric(more.40.v[[group]]$perc)         
        more.40.v[[group]]$n <- as.numeric(more.40.v[[group]]$n)         
        
        less.40.v[[group]] <- cbind(less.40.v[[group]], moss = rep("< 40aa",n = nrow(less.40.v)))
        more.40.v[[group]] <- cbind(more.40.v[[group]], moss = rep("> 40aa",n = nrow(more.40.v)))
        
        less.40.v[[group]] <- less.40.v[[group]][order(less.40.v[[group]]$v_call),]
        more.40.v[[group]] <- more.40.v[[group]][order(more.40.v[[group]]$v_call),]
        
        less.40.v[[group]]$position <- less.40.v[[group]]$perc/2
        more.40.v[[group]]$position <- less.40.v[[group]]$perc + more.40.v[[group]]$perc/2
        
        too.close <- abs(less.40.v[[group]]$position - more.40.v[[group]]$position) < 1

        more.40.v[[group]][too.close,]$position <- 1.25 
        
        moss.40.v[[group]] <- rbind(more.40.v[[group]], less.40.v[[group]])
        moss.40.v[[group]]$moss <- factor(moss.40.v[[group]]$moss, levels = c("> 40aa","< 40aa"))
        moss.40.v[[group]]$perc <- as.numeric(moss.40.v[[group]]$perc)

        plt <- ggplot(data = moss.40.v[[group]], aes(x = v_call, y = perc, fill = moss))+
                geom_bar(stat ="identity")+
                theme_light()+
                labs(title = paste("IGHV distribution of ", group, sep =""),
                     y = "Reads (%)",
                     x = paste("IGHV (N = ",nbr.row[[group]],")",sep=""),
                     fill = "Lenght of cdr3")+
                scale_fill_manual(values=c("lightblue3","lightcoral"))+
                geom_label(aes(y = position, label = n), color = "black", fill = c(rep("lightblue2",12),rep("lightcoral",12)))
    
            print(plt)
            
        plt <- ggplot(data = less.40.v[[group]], aes(x = v_call, y = perc))+
          geom_bar(stat = "identity")+
          theme_light()+
          labs(title = paste("IGHV distribution of ", group," with cdr3 <40 aa",sep =""),
               y = "Reads (%)",
               x = paste("IGHV (total nbr of cdr3 <40 aa = ", sum(less.40.v[[group]]$n),")",sep=""))+
          geom_label(aes(y = perc/2, label = n), color = "black")
        print(plt)
        
        plt <- ggplot(data = more.40.v[[group]], aes(x = v_call, y = perc))+
          geom_bar(stat = "identity")+
          theme_light()+
          labs(title = paste("IGHV distribution of ", group," with cdr3 >40 aa",sep =""),
               y = "Reads (%)",
               x = paste("IGHV (total nbr of cdr3 >40 aa = ", sum(more.40.v[[group]]$n),")",sep=""))+
          geom_label(aes(y = perc/2, label = n), color = "black")
    print(plt)
    dev.off()
    }
    


    
    
    
     #### Barplot IGHV CDR3 < 40 NO####
     #for(x in names(igblast)){
     #  temp_plot <- ggplot(igblast[[x]][nchar(igblast[[x]]$cdr3_aa) <= 40,],
     #                      aes(x = igblast[[x]][nchar(igblast[[x]]$cdr3_aa) <= 40,]$v_call))+ 
     #    geom_bar(color = "black", fill = "darkblue") +
     #    labs(title="Distribution des IGHV avec CDR3 <= 40aa",
     #         subtitle = x,
     #         y="Nombre de sequences",
     #         x="IGHV")+
     #    theme(title =element_text(size=12, face='bold'),
     #          axis.title.x = element_text(vjust = 0, size = 15),
     #          axis.title.y = element_text(vjust = 2, size = 15),
     #          axis.text = element_text(color = "black", face = "bold", size = 14),
     #          axis.text.x = element_text(face = "bold", size = 13))+
     #    scale_x_discrete(limits = IGHV_possible)+
     #    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
     #  print(temp_plot)}
     #print("Barplot IGHV CDR3 < 40 : DONE")
  

  
  
  #longueur_cdr3 <- sapply(names(igblast),
  #                        function(x) nchar(igblast[[x]]$cdr3_aa),
  #                        simplify = FALSE, USE.NAMES = TRUE)
  #data_longueur_cdr3 <- data.frame()
  #temp.df <- sapply(names(longueur_cdr3), 
  #                  function(x) data.frame(longueur_cdr3 = longueur_cdr3[[x]],
  #                                         nom = rep(x, length(longueur_cdr3[[x]]))),
  #                  simplify = FALSE, USE.NAMES = TRUE)
  
  
  #for(i in (1:length(temp.df))){
  #  data_longueur_cdr3 <- rbind(data_longueur_cdr3, temp.df[[i]])
  #}
  #
  #data_longueur_cdr3$nom <-  factor(data_longueur_cdr3$nom, unique(data_longueur_cdr3$nom))
  #print(318)
  #catego <- str_sub(str_extract(data_longueur_cdr3$nom, "-[:digit:]{3,4}"),2,5)
  #
  #temp_plot <- ggplot(data_longueur_cdr3, aes(y = longueur_cdr3 ,x = nom, fill = catego))+
  #  geom_violin()+
  #  geom_boxplot(width=0.1, fill = "light blue")+
  #  labs(title="Longueur de CDR3",
  #       subtitle = group,
  #       y="Longueur CDR3 (aa)",
  #       x="Echantillon")+
  #  theme(title =element_text(size=12, face='bold'),
  #        axis.title.x = element_text(vjust = 0, size = 15),
  #        axis.title.y = element_text(vjust = 2, size = 15),
  #        axis.text = element_text(color = "black", face = "bold", size = 14),
  #        axis.text.x = element_text(face = "bold", size = 13))+
  #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  #print(temp_plot)
  #print("Boxplot Longueur CDR3 : DONE")
  #dev.off()
  ##### NGMerge#######
  pdf("graph2.pdf",15,15)
  
  
  ##### Loading Files needed for the end
   
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
  
  
  
  
  
#  overlap <- sapply(names(ngmerge), function(x) ngmerge[[x]]$OverlapLen,
#                    simplify = FALSE, USE.NAMES = TRUE)
#  data_overlap <- data.frame(longueur_overlap = NA, nom = NA)
#  temp.df <- sapply(names(longueur),
#                    function(x) data.frame(longueur_overlap = overlap[[x]],
#                                           nom = rep(x, length(overlap[[x]]))),
#                    simplify = FALSE, USE.NAMES = TRUE)
#  for(i in (1:length(temp.df))) data_overlap <- rbind(data_overlap, temp.df[[i]])
#  
#  catego <- str_sub(str_extract(data_overlap$nom, "-+...?"),2,4)
#  temp_plot <- ggplot(data = data_overlap[!is.na(data_overlap$nom),], mapping = aes(x = nom, y = longueur_overlap))+
#    geom_violin(color = "black", fill = "darkblue")+
#    geom_boxplot(width=0.1, fill = "light blue")+
#    stat_summary(fun=mean, color = "black", geom="point", shape=23, size=2)+
#    labs(title = "Longueur de l'overlap lors de NGMerge",
#         y = "Longueur de l'overlap (nt)",
#         x = "Echantillon")+
#    theme(title = element_text(size = 12, face = 'bold'),
#          axis.title.x = element_text(vjust = 0, size = 15),
#          axis.title.y = element_text(vjust = 2, size = 15),
#          axis.text    = element_text(color = "black", face = "bold", size = 14),
#          axis.text.x  = element_text(face = "bold", size = 13))+
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
#    scale_fill_discrete(name = "Cow #")
#  print(temp_plot)
#  print("Boxplot Overlap NGMerge : DONE")
#
#  dev.off()
#  print("dev.off")
#  remove(igblast)
#  print("IgBlast Removed")
#  remove(ngmerge)
#  print("Ngmerge removed")
#  gc()
#  print("gc()")
#
sessionInfo()

