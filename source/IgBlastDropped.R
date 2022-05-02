
# Author : Guyllaume Dufresne
# Description : Generate graph of dropped reads
# Last update : 01/05/22

packages <- c("ggplot2", "viridis", "ShortRead", "stringr", "tidyr", "tibble", "dplyr","magrittr")

package.check <- lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(x)
    library(x, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
  }
})

print("Package loaded!")

# Get the position for the label
get.pos <- function(readClass){ 
  c(tot.per.ig[readClass] / 2 +
      frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "Trimmomatic", ]$Percent +
      frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "NGMerge", ]$Percent,
    frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "NGMerge", ]$Percent / 2 +
      frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "Trimmomatic", ]$Percent,
    frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "Trimmomatic", ]$Percent / 2)
}

# Make the text for the labels
labelr <- function(readClass){
  c(tot.per.ig[readClass],
    frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "NGMerge", ]$Percent,
    frct.lost.ordered[frct.lost.ordered$Class == readClass & frct.lost.ordered$Reason == "Trimmomatic", ]$Percent)
}



args <- commandArgs(trailingOnly = TRUE)

# For testing
  #args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")

id_all <- sort(unique(args))

# Create the groups of reads
groupement <- list(
  "G1" = grep(pattern = "IgG1", x = id_all),
  "G2" = grep(pattern = "IgG2", x = id_all),
  "G3" = grep(pattern = "IgG3", x = id_all),
  "GM1" = grep(pattern = "IgGM1", x = id_all),
  "GM2" = grep(pattern = "IgGM2", x = id_all)
)


# read the result of igBlast
igblast.lst <- list()
for (group in names(groupement)) {
  echantillons <- id_all[groupement[[group]]]
  chemins_igblast <- sapply(echantillons, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv1", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )
  igblast.lst[[group]] <- sapply(chemins_igblast, function(x) read.csv(x, sep = "\t", na.strings = c("", "NA")),
    simplify = FALSE, USE.NAMES = TRUE
  )
}

# get the total number of reads
total.reads.lst <- list()
for (group in names(groupement)) {
  chemin_raw_reads <- sapply(id_all[groupement[[group]]], function(x) paste("./data/rawReads/", x, "_L001_R1_001.fastq.gz", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )

  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
    simplify = FALSE, USE.NAMES = TRUE
  )
  total.reads.lst[group] <- sum(unlist(nbr_raw_reads))
}

# turn the list obtained into a dataframe
t(do.call(cbind.data.frame, total.reads.lst))[,1] %>%
  data.frame(nbr.reads = ., class = c("G1", "G2", "G3", "GM1", "GM2")) -> total.reads.df

# Get the number of reads before and after igblast
# Also get the reason for dropping the read
num.row.ig.bef <- list()
num.row.ig.aft <- list()
num.row.aft <- c()
prc.drop.reason.tot <- list()
temp.prc.drop.reason <- c(0, 0, 0, 0)

for (j in names(igblast.lst)) {
  # Select a group of reads
  temp.igblast <- igblast.lst[[j]]
  # Get the number of reads in the group
  num.row.ig.bef[j] <- sum(sapply(temp.igblast, nrow))
  
  # remove unwanted reads
  for (i in (1:length(temp.igblast))) {
    # get percent of lost reads from the presence of stop codon
    stop.codon <- nrow(temp.igblast[[i]][temp.igblast[[i]]$stop_codon == TRUE,   ]) / total.reads.df[j, 1] * 100
    # remove reads with stop codon
    temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$stop_codon), ]
    temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$stop_codon == FALSE,]

    # get percent of lost reads from the presence of a frameshift (indel)
    frameshift <- nrow(temp.igblast[[i]][temp.igblast[[i]]$v_frameshift == TRUE,   ]) / total.reads.df[j, 1] * 100
    # remove reads with a frameshit
    temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$v_frameshift), ]
    temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$v_frameshift == FALSE,]

    # get percent of lost reads from IGHV and IGHJ not in frame
    vj_in_frame <- nrow(temp.igblast[[i]][temp.igblast[[i]]$vj_in_frame == FALSE, ]) / total.reads.df[j, 1] * 100
    # remove reads not in frame
    temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$vj_in_frame), ]
    temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$vj_in_frame == TRUE, ]

    # get percent of lost reads form incomplete vdj (lacks the first nt of IGHV or last nt of IGHJ)
    complete.vdj <- nrow(temp.igblast[[i]][temp.igblast[[i]]$complete_vdj == FALSE, ]) / total.reads.df[j, 1] * 100
    # remove incomplete reads
    temp.igblast[[i]] <- temp.igblast[[i]][!is.na(temp.igblast[[i]]$complete), ]
    temp.igblast[[i]] <- temp.igblast[[i]][temp.igblast[[i]]$complete_vdj == TRUE, ]

    # Number of reads after the removal
    num.row.aft[i] <- nrow(temp.igblast[[i]])
    temp.prc.drop.reason <- temp.prc.drop.reason + c(stop.codon, frameshift, vj_in_frame, complete.vdj)
  }
  num.row.ig.aft[[j]] <- sum(unlist(num.row.aft))
  prc.drop.reason.tot[[j]] <- temp.prc.drop.reason
  temp.prc.drop.reason <- c(0, 0, 0, 0)
}

#turn the lists into dataframe
prc.drop.reason <- cbind(do.call(cbind.data.frame, prc.drop.reason.tot), Reason = c("stop.codon", "frameshift", "vj_in_frame", "complete.vdj"))
num.row.ig <- data.frame(ig.bef = unlist(num.row.ig.bef), ig.aft = unlist(num.row.ig.aft))

#save memory
remove(igblast.lst, num.row.aft, num.row.ig.bef, num.row.ig.aft, temp.prc.drop.reason, temp.igblast)

print("removed unwanted reads")
#remove(igblast.lst, num.row.aft, num.row.ig.bef, num.row.ig.aft, temp.prc.drop.reason, temp.igblast)




# get the number reads left after trimmomatic
num.row.trim.aft <- list()
for (group in names(groupement)) {
  echantillons <- id_all[groupement[[group]]]
  chemins_trimmomatic <- sapply(echantillons, function(x) paste("./data/trimmedReads/", x, "_1P.fastq", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )
  trimmomatic <- sapply(chemins_trimmomatic, function(x) readFastq(x),
    simplify = FALSE, USE.NAMES = TRUE
  )
  nbr_read_trim <- sapply(names(trimmomatic), function(x) length(trimmomatic[[x]]),
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  num.row.trim.aft[[group]] <- sum(unlist(nbr_read_trim))
}
# Turn the list into a dataframe
num.row.trim <- data.frame(trim.bef = total.reads.df$nbr.reads, trim.aft = unlist(num.row.trim.aft))


# Get the number of rows befor and after NGmerge
ng.lst <- list()
num.row.bef <- list()
for (group in names(groupement)) {
  echantillons <- id_all[groupement[[group]]]

  chemins_ngmerge <- sapply(echantillons, function(x) paste("./data/mergedReads/", x, ".log", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )
  ngmerge <- sapply(chemins_ngmerge, function(x) read.csv(x, sep = "\t", na.strings = c("NA")),
    simplify = FALSE, USE.NAMES = TRUE
  )

  num.row.bef[group] <- sum(sapply(names(ngmerge), function(x) nrow(ngmerge[[x]]),
    simplify = TRUE,
    USE.NAMES = TRUE
  ))
  ngmerge_no_na <- sapply(names(ngmerge), function(x) ngmerge[[x]][!is.na(ngmerge[[x]]$OverlapLen), ],
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  num_read_ngmerge <- sapply(names(ngmerge_no_na), function(x) nrow(ngmerge_no_na[[x]]),
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  ng.lst[[group]] <- sum(unlist(num_read_ngmerge))
}
num.row.ng <- data.frame(NG.bef = unlist(num.row.bef), NG.aft = unlist(ng.lst))

# Put all the data of number of reads into one dataframe
num.row.all <- cbind(num.row.trim, num.row.ng, num.row.ig)

axis.names <- c(
  paste("G1 \n N = ", num.row.all["G1", "trim.bef"], sep = ""),
  paste("G2 \n N = ", num.row.all["G2", "trim.bef"], sep = ""),
  paste("G3 \n N = ", num.row.all["G3", "trim.bef"], sep = ""),
  paste("GM1 \n N = ", num.row.all["GM1", "trim.bef"], sep = ""),
  paste("GM2 \n N = ", num.row.all["GM2", "trim.bef"], sep = "")
)


# get the percent of lost reads
prc.lost.trim <- (1 - (num.row.all$trim.aft / num.row.all$trim.bef)) * 100

prc.lost.ng <- (num.row.all$trim.aft - num.row.all$NG.aft) / num.row.all$trim.bef * 100

prc.lost.ig.frame <- prc.drop.reason[prc.drop.reason$Reason == "frameshift", ]
prc.lost.ig.incomplete <- prc.drop.reason[prc.drop.reason$Reason == "complete.vdj", ]
prc.lost.ig.vjFrame <- prc.drop.reason[prc.drop.reason$Reason == "vj_in_frame", ]
prc.lost.ig.Stop <- prc.drop.reason[prc.drop.reason$Reason == "stop.codon", ]


frct_lost <- cbind(
  trimmomatic = prc.lost.trim,
  ngmerge = prc.lost.ng
)

row.names(frct_lost) <- c("G1", "G2", "G3", "GM1","GM2")

# make dataframe of reason of drop
as.data.frame(t(frct_lost)) %>%
gather(data = ., key = Class, value = Percent) %>%
cbind(., Reason = rep(c("Trimmomatic", "NGMerge"), 5)) -> frct_lost

frct_lost <- rbind(frct_lost, cbind(rbind(gather(data = prc.lost.ig.frame, key = Class, value = Percent)[-6, ],
                                          gather(data = prc.lost.ig.incomplete, key = Class, value = Percent)[-6, ],
                                          gather(data = prc.lost.ig.vjFrame, key = Class, value = Percent)[-6, ],
                                          gather(data = prc.lost.ig.Stop, key = Class, value = Percent)[-6, ]),
                                    Reason = c(rep("Frameshift", 5),
                                               rep("Incomplete VDJ", 5),
                                               rep("VJ not in frame", 5),
                                               rep("Stop Codon", 5))))
frct_lost$Percent <- as.numeric(frct_lost$Percent)
frct_lost$Reason <- factor(frct_lost$Reason, levels = c("Frameshift", "Incomplete VDJ", "Stop Codon", "VJ not in frame", "NGMerge", "Trimmomatic"))

pdf("graph/graphIgBlastDropped.pdf", width = 10, height = 10)



frct.lost.ordered <- frct_lost %>% arrange(Class, Reason)
frct.lost.ordered$Percent <- as.numeric(frct.lost.ordered$Percent)

# Get the total percent for igblast drop
tot.per.ig <- c(
  G1= sum(frct.lost.ordered[frct.lost.ordered$Class == "G1" & 
                             frct.lost.ordered$Reason != "NGMerge" &
                             frct.lost.ordered$Reason != "Trimmomatic", ]$Percent),
  G2 = sum(frct.lost.ordered[frct.lost.ordered$Class == "G2" & 
                             frct.lost.ordered$Reason != "NGMerge" &
                             frct.lost.ordered$Reason != "Trimmomatic", ]$Percent),
  G3 = sum(frct.lost.ordered[frct.lost.ordered$Class == "G3" & 
                             frct.lost.ordered$Reason != "NGMerge" &
                             frct.lost.ordered$Reason != "Trimmomatic", ]$Percent),
  GM1= sum(frct.lost.ordered[frct.lost.ordered$Class == "GM1" & 
                             frct.lost.ordered$Reason != "NGMerge" &
                             frct.lost.ordered$Reason != "Trimmomatic", ]$Percent),
  GM2= sum(frct.lost.ordered[frct.lost.ordered$Class == "GM2" & 
                             frct.lost.ordered$Reason != "NGMerge" &
                             frct.lost.ordered$Reason != "Trimmomatic", ]$Percent))


# Set the labels for the graph
y.pos <- data.frame(
  Class = sort(rep(c("G1", "G2", "G3", "GM1", "GM2"), 3)),
  Reason = rep(c("IgBlast", "NGMerge", "Trimmomatic"), 5),
  Position = c(get.pos("G1"),get.pos("G2"),get.pos("G3"),get.pos("GM1"),get.pos("GM2")),
  Label = c(labelr("G1"),labelr("G2"),labelr("G3"),labelr("GM1"),labelr("GM2")))


#change the number format for labels
y.pos$Label <- formatC(y.pos$Label, digits = 2, format = "f")
y.pos$Label <- paste(y.pos$Label, "%", sep = "")

#generate the graph with all reasons
plt <- ggplot() +
  geom_bar(data = frct_lost, aes(x = Class, y =Percent, fill = Reason), stat = "identity", position = "stack") +
  theme_light() +
  labs(
    title = "Reason for drop",
    y = "Percent of sequences dropped",
    x = "Class") +
  ylim(0, 100)+
  scale_fill_manual(values = c(
    magma(10)[3],
    magma(10)[4],
    magma(10)[5],
    magma(10)[6],
    "darkseagreen",
    "cornflowerblue"))+
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(vjust = 0, size = 15),
    axis.title.y = element_text(vjust = 2, size = 15),
    axis.text = element_text(color = "black", face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 13))+
  scale_x_discrete(labels = axis.names)+
  geom_label(data = y.pos, aes(x = Class, y = Position, label = Label))

print(plt)



# setup data for graph of only dropped reads (no reason)
prc.lost <- (num.row.all$trim.bef - num.row.all$ig.aft) / num.row.all$trim.bef

prc.lost.df <- data.frame(
  percent = prc.lost * 100,
  class = c("G1", "G2", "G3", "GM1", "GM2"),
  nbr = (num.row.all$trim.bef - num.row.all$ig.aft)
)

prc.lost.df <- data.frame(prc.lost.df,
                          label = paste(prc.lost.df$nbr, "\n ", 
                                        formatC(prc.lost.df$percent, digits = 2, format = "f"),
                                        "%",sep = ""))

# generate graph of only dropped reads (no reason)
plt <- ggplot(data = prc.lost.df, aes(x = class, y = percent)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(
    title = "Percent of dropped sequences",
    y = "Dropped sequences (%)",
    x = "Class") +
  theme_light() +
  ylim(0, 100) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(vjust = 0, size = 15),
    axis.title.y = element_text(vjust = 2, size = 15),
    axis.text = element_text(color = "black", face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 13)) +
  scale_x_discrete(labels = axis.names) +
  geom_label(data = prc.lost.df, aes(x = class, y = percent / 2, label = label))

print(plt)


# setup for graph of number of totals reads
options(scipen = 100000)

#generate graph of total reads
plt <- ggplot() +
  geom_bar(data = total.reads.df, aes(x = class, y = nbr.reads), stat = "identity", color = "black") +
  labs(
    title = "Total reads",
    y = "Number of reads",
    x = "Class") +
  theme_light() +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(vjust = 0, size = 15),
    axis.title.y = element_text(vjust = 2, size = 15),
    axis.text = element_text(color = "black", face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 13)) +
  geom_label(data = total.reads.df, aes(x = class, y = nbr.reads / 2, label = nbr.reads), stat = "identity", color = "black")

print(plt)


dev.off()
sessionInfo()
