###########################
#                         #
# Guyllaume Dufresne      #
# Graph sur les séquences #
#                         #
###########################



packages = c("ggplot2","viridis","ShortRead", "stringr", "tidyr", "tibble")

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

print("Package loaded!")

name_clean_up <- function(string){
  out <- substring(string, first=1 ,last = (str_locate(string = string, pattern = ","))[1] -1)
  if(is.na(out)){
    string
  } else {
    out
  }
}

args <- commandArgs(trailingOnly=TRUE)
#args <-  c("25-271-IgGM2-1_S201","26-258-IgG1-1_S114","26-258-IgG1-1_S114","26-258-IgG2-1_S158","26-258-IgG2-1_S158","26-258-IgG3-1_S26","26-258-IgG3-1_S26","26-258-IgGM1-1_S70","26-258-IgGM1-1_S70","26-258-IgGM2-1_S202","20-261-IgG2-1_S152","20-261-IgG2-1_S152","20-261-IgG3-1_S20","20-261-IgG3-1_S20","20-261-IgGM1-1_S64","20-261-IgGM1-1_S64","20-261-IgGM2-1_S196","20-261-IgGM2-1_S196","21-262-IgG1-1_S109","21-262-IgG1-1_S109","21-262-IgG2-1_S153","21-262-IgG2-1_S153","22-264-IgGM2-1_S198","22-264-IgGM2-1_S198","23-267-IgG1-1_S111","23-267-IgG1-1_S111","23-267-IgG2-1_S155","23-267-IgG2-1_S155","23-267-IgG3-1_S23","23-267-IgG3-1_S23","23-267-IgGM1-1_S67","23-267-IgGM1-1_S67","23-267-IgGM2-1_S199","23-267-IgGM2-1_S199","6-253-IgG3-1_S6","6-253-IgG3-1_S6","6-253-IgGM1-1_S50","6-253-IgGM1-1_S50","6-253-IgGM2-1_S182","6-253-IgGM2-1_S182","7-603-IgG1-1_S95","7-603-IgG1-1_S95","7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")
#args <- c("7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","26-258-IgGM1-1_S70","35-603-IgG1-1_S123","6-253-IgG3-1_S6","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","20-261-IgGM1-1_S64","6-253-IgG2-1_S138","7-603-IgGM1-1_S51","21-262-IgG2-1_S153","26-258-IgG2-1_S158","23-267-IgG3-1_S23")
id_all <-  sort(unique(args))

groupement <- list("G1" = grep(pattern = "IgG1", x = id_all),
                   "G2" = grep(pattern = "IgG2", x = id_all),
                   "G3" = grep(pattern = "IgG3", x = id_all),
                   "GM1" = grep(pattern = "IgGM1", x = id_all),
                   "GM2" = grep(pattern = "IgGM2", x = id_all))

  igblast.lst <- list()
  ## IgBlast
  for(group in names(groupement)){
    
  echantillons <- id_all[groupement[[group]]]
  chemins_igblast <- sapply(echantillons, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv1", sep = ""),
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
  


  num.row.ig.bef <- list()
  num.row.ig.aft <- list()
  num.row.aft <- c()
  prc.drop.reason.tot <- list()
  temp.prc.drop.reason <- c(0,0,0,0)
  nas <- list()
  nas.all <- list()
  nas.all.2 <- c(0,0,0,0)
  for(j in names(igblast.lst)){
    temp.igblast <- igblast.lst[[j]]
    num.row.ig.bef[j] <- sum(sapply(temp.igblast,nrow))
    for(i in (1:length(temp.igblast))){
      nas["Codon"] <- sum(is.na(temp.igblast[[i]]$stop_codon))
      stop.codon <- nrow(temp.igblast[[i]][temp.igblast[[i]]$stop_codon   == TRUE ,])/total.reads.df[j,1] * 100
      temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$stop_codon),]
      temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$stop_codon == FALSE,]
      nas["Frameshift"] <- sum(is.na(temp.igblast[[i]]$v_frameshift))
      frameshift <- nrow(temp.igblast[[i]][temp.igblast[[i]]$v_frameshift == TRUE ,])/total.reads.df[j,1] * 100
      temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$v_frameshift),]
      temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$v_frameshift == FALSE,]
      nas["Productive"] <- sum(is.na(temp.igblast[[i]]$productive))
      productive <- nrow(temp.igblast[[i]][temp.igblast[[i]]$productive   == FALSE,])/total.reads.df[j,1] * 100
      temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$productive),]
      temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$productive == TRUE,]
      nas["Complete"] <- sum(is.na(temp.igblast[[i]]$complete_vdj))
      complete.vdj <- nrow(temp.igblast[[i]][temp.igblast[[i]]$complete_vdj == FALSE,])/total.reads.df[j,1] * 100
      temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$complete),]
      temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$complete_vdj  == TRUE,]
      
      num.row.aft[i] <- nrow(temp.igblast[[i]])
      temp.prc.drop.reason <- temp.prc.drop.reason + c(stop.codon,frameshift,productive,complete.vdj)
    }
    nas.all[[j]] <- unlist(list(nas["Codon"],nas["Frameshift"],nas["Productive"],nas["Complete"]))
    num.row.ig.aft[[j]] <- sum(unlist(num.row.aft))
    prc.drop.reason.tot[[j]] <- temp.prc.drop.reason
    temp.prc.drop.reason <- c(0,0,0,0)
    
  }

  prc.drop.reason <- cbind(do.call(cbind.data.frame, prc.drop.reason.tot),Reason = c("stop.codon","frameshift","productive","complete.vdj")) 
  num.row.ig <- data.frame(ig.bef = unlist(num.row.ig.bef), ig.aft = unlist(num.row.ig.aft))
  remove(igblast.lst, num.row.aft, num.row.ig.bef,num.row.ig.aft,temp.prc.drop.reason,temp.igblast)
  
  
  pdf("graph/graphIgBlastDropped.pdf",width=10, height=10)
  
  axis.names =c(paste("G1 \n n=", total.reads.df[total.reads.df$class == "G1",]$nbr.reads, sep = ""),
                paste("G2 \n n=", total.reads.df[total.reads.df$class == "G2",]$nbr.reads, sep = ""),
                paste("G3 \n n=", total.reads.df[total.reads.df$class == "G3",]$nbr.reads, sep = ""),
                paste("GM1 \n n=",total.reads.df[total.reads.df$class == "GM1",]$nbr.reads, sep = ""),
                paste("GM2 \n n=",total.reads.df[total.reads.df$class == "GM2",]$nbr.reads, sep = ""))
  

num.row.trim.aft <- list()
for(group in names(groupement)){ 
  
  echantillons <- id_all[groupement[[group]]]
  chemins_trimmomatic <- sapply(echantillons, function(x) paste("./data/trimmedReads/", x, "_1P.fastq", sep =""),
                                simplify = FALSE, USE.NAMES = TRUE)
  trimmomatic <- sapply(chemins_trimmomatic, function(x) readFastq(x),
                        simplify = FALSE, USE.NAMES = TRUE)
  nbr_read_trim <- sapply(names(trimmomatic), function(x) length(trimmomatic[[x]]),
                          simplify = FALSE,
                          USE.NAMES = TRUE)

  num.row.trim.aft[[group]] <- sum(unlist(nbr_read_trim))
}

num.row.trim <- data.frame(trim.bef = total.reads.df$nbr.reads, trim.aft = unlist(num.row.trim.aft))

#### histogramme sequences non utilise NO ####
ng.lst <- list()
num.row.bef <- list()
for(group in names(groupement)){
  print(group)
  echantillons <- id_all[groupement[[group]]]
  chemins_ngmerge <- sapply(echantillons, function(x) paste("./data/mergedReads/", x, ".log", sep =""),
                            simplify = FALSE, USE.NAMES = TRUE)
  ngmerge <-  sapply(chemins_ngmerge, function(x) read.csv(x, sep= "\t", na.strings=c("NA")),
                     simplify = FALSE, USE.NAMES = TRUE)
  
  num.row.bef[group] <- sum(sapply(names(ngmerge), function(x) nrow(ngmerge[[x]]),
                               simplify = TRUE,
                               USE.NAMES = TRUE))
  ngmerge_no_na <- sapply(names(ngmerge), function(x) ngmerge[[x]][!is.na(ngmerge[[x]]$OverlapLen),],
                        simplify = FALSE,
                        USE.NAMES = TRUE)
  num_read_ngmerge <- sapply(names(ngmerge_no_na), function(x) nrow(ngmerge_no_na[[x]]),
                             simplify = FALSE,
                             USE.NAMES = TRUE)

  ng.lst[[group]] <- sum(unlist(num_read_ngmerge))
}

num.row.ng <- data.frame(NG.bef = unlist(num.row.bef), NG.aft = unlist(ng.lst))


num.row.all <- cbind(num.row.trim,num.row.ng,num.row.ig)

prc.lost.trim <-  (1 - (num.row.all$trim.aft/num.row.all$trim.bef)) * 100

prc.lost.ng   <-  (num.row.all$trim.aft - num.row.all$NG.aft)/num.row.all$trim.bef * 100

prc.lost.ig.frame         <- prc.drop.reason[prc.drop.reason$Reason == "frameshift",]
prc.lost.ig.incomplete    <- prc.drop.reason[prc.drop.reason$Reason == "complete.vdj",]
prc.lost.ig.nonProductive <- prc.drop.reason[prc.drop.reason$Reason == "productive",]
prc.lost.ig.Stop          <- prc.drop.reason[prc.drop.reason$Reason == "stop.codon",]
prc.na.ig <- sapply(nas.all, function(x) sum(x)) 

frct_lost <- cbind(trimmomatic = prc.lost.trim,
                   ngmerge     = prc.lost.ng)
row.names(frct_lost) <- c("G1","G2","G3","GM1","GM2")
frct_lost <- as.data.frame(t(frct_lost))
frct_lost <- gather(data=frct_lost, key = Class, value = Percent)
frct_lost <- cbind(frct_lost, Reason = rep(c("Trimmomatic","NGMerge"), 5))

temp <- rbind(gather(data=prc.lost.ig.frame, key=Class, value = Percent)[-6,],
                   gather(data=prc.lost.ig.incomplete, key=Class, value = Percent)[-6,],
                   gather(data=prc.lost.ig.nonProductive, key=Class, value = Percent)[-6,],
                   gather(data=prc.lost.ig.Stop, key=Class, value = Percent)[-6,])
frct_lost <- rbind(frct_lost, cbind(temp, Reason = c(rep("Frameshift",5),rep("Incomplete VDJ",5),rep("Non Productive",5),rep("Stop Codon",5))))
frct_lost$Reason <- factor(frct_lost$Reason, levels = c("Frameshift","Incomplete VDJ","Stop Codon","Non Productive","NGMerge","Trimmomatic"))

frct_lost$Percent <- as.numeric(frct_lost$Percent)

ggplot(data=frct_lost, aes(x = Class, y = Percent, fill = Reason))+
  geom_bar(stat="identity", position = "stack")


plt <- ggplot(data=frct_lost, aes(x = Class, y = Percent, fill = Reason))+
        geom_bar(stat= "identity",position = "stack")+
        theme_light()+
        labs(title = "Reason for drop",
             y = "Percent of sequences dropped",
             x = "Class") +
        ylim(0,100) +
        scale_fill_manual(values=c(magma(10)[3],
                                   magma(10)[4],
                                   magma(10)[5],
                                   magma(10)[6],
                                   "darkseagreen",
                                   "cornflowerblue"))+
        theme(title =element_text(size=12, face='bold'),
              axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text    = element_text(color = "black", face = "bold", size = 14),
              axis.text.x  = element_text(face = "bold", size = 13))+
        scale_x_discrete(labels = axis.names)
print(plt)

plt <- ggplot()+
        geom_bar(data = total.reads.df, aes(x = class, y = nbr.reads), stat = "identity")+
        labs(title = "Valid reads",
              y = 'Number of reads',
              x = "Class")+
        theme_light() + 
        theme(title =element_text(size=12, face='bold'),
              axis.title.x = element_text(vjust = 0, size = 15),
              axis.title.y = element_text(vjust = 2, size = 15),
              axis.text    = element_text(color = "black", face = "bold", size = 14),
              axis.text.x  = element_text(face = "bold", size = 13))+
        scale_x_discrete(labels = axis.names)
print(plt)


prc.lost <- (num.row.all$trim.bef - num.row.all$ig.aft)/num.row.all$trim.bef

prc.lost.df <- data.frame(percent = prc.lost*100, class = c("G1","G2","G3","GM1","GM2"))  


plt<- ggplot(data = prc.lost.df, aes(x = class, y = percent)) +
       geom_bar(stat = "identity", position = "dodge", color = "black")+
       labs(title = "Percent of dropped sequences",
            y = "Dropped sequences (%)",
            x = "Class") +
       theme_light()+
       ylim(0, 100) +
       theme(title =element_text(size=12, face='bold'),
          axis.title.x = element_text(vjust = 0, size = 15),
          axis.title.y = element_text(vjust = 2, size = 15),
          axis.text    = element_text(color = "black", face = "bold", size = 14),
          axis.text.x  = element_text(face = "bold", size = 13))+
  theme(title =element_text(size=12, face='bold'),
        axis.title.x = element_text(vjust = 0, size = 15),
        axis.title.y = element_text(vjust = 2, size = 15),
        axis.text    = element_text(color = "black", face = "bold", size = 14),
        axis.text.x  = element_text(face = "bold", size = 13))+
  scale_x_discrete(labels = axis.names)
print(plt)

dev.off()
sessionInfo()
  
