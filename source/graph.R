
# Author : Guyllaume Dufresne
# Description : Generate graphs
# Last update : 01/05/22

packages <- c("ggplot2", "viridis", "ShortRead", "stringr", "tidyr", "tibble", "dplyr")

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
print("Package loaded succesfully!")

#### Functions ####

# Only keeps the first result of igBlast
name_clean_up <- function(string) {
  out <- substring(string, first = 1, last = (str_locate(string = string, pattern = ","))[1] - 1)
  if (is.na(out)) {
    string
  } else {
    out
  }
}


#### Arguments ####
args <- commandArgs(trailingOnly = TRUE)
# For testing
  #args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")

id_all <- sort(unique(args))



groupement <- list(
  "G1" = grep(pattern = "IgG1", x = id_all),
  "G2" = grep(pattern = "IgG2", x = id_all),
  "G3" = grep(pattern = "IgG3", x = id_all),
  "GM1" = grep(pattern = "IgGM1", x = id_all),
  "GM2" = grep(pattern = "IgGM2", x = id_all),
  "All" = grep(pattern = "", x = id_all)
)


# Make a list of list of reads
igblast.lst <- list()
for (group in names(groupement)) {
  echantillons <- id_all[groupement[[group]]]
  chemins_igblast <- sapply(echantillons, function(x) paste("./csv/", x, "_length.csv", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )
  igblast.lst[[group]] <- sapply(chemins_igblast, function(x) read.csv(x, na.strings = c("", "NA")),
    simplify = FALSE, USE.NAMES = TRUE
  )
}


# Get the total number of reads in each sample
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

total.reads.df.1 <- do.call(cbind.data.frame, total.reads.lst)
#total.reads.df.2 <- data.frame(nbr.reads = t(total.reads.df.1)[, 1], class = c("G1", "G2", "G3", "GM1", "GM2", "All"))

IGHV_possible <- c("IGHV1-7", "IGHV1-10", "IGHV1-14", "IGHV1-17", "IGHV1-21/33", "IGHV1-25", "IGHV1-27", "IGHV1-30", "IGHV1-37", "IGHV1-39", "IGHV1-20", "IGHV1-32")
IGHD_possible <- c("IGHD1-1", "IGHD1-2/4", "IGHD1-3", "IGHD1-4", "IGHD2-1/2/3/4", "IGHD3-1/3/4", "IGHD4-1", "IGHD5-2", "IGHD5-3/4", "IGHD6-2", "IGHD6-3/4", "IGHD7-3", "IGHD7-4", "IGHD8-2", "IGHD9-1/4")
IGHJ_possible <- c("IGHJ1-4", "IGHJ1-6", "IGHJ2-4")

#### Violinplot longueur sequence ####
longueur.lst <- list()
for (group in names(igblast.lst)) {
  temp <- igblast.lst[[group]]
  for (i in 1:length(temp)) {
    temp[i]
    longueur.lst[[group]] <- append(longueur.lst[[group]], sapply(temp[i], function(x) x$sequence_alignment,
      simplify = TRUE, USE.NAMES = TRUE
    ))
  }
}

toute <- list()
for (i in 1:length(longueur.lst)) {
  toute[[i]] <- data.frame(longueur = longueur.lst[[i]], class = names(longueur.lst[i]))
}

longueur.df <- do.call(rbind.data.frame, toute)

nbr.row <- list()
for (group in names(groupement)) {
  temp <- na.omit(igblast.lst[[group]][[1]])
  nbr.row[[group]] <- nrow(temp)
  for (i in 2:length(igblast.lst[[group]])) {
    temp <- na.omit(igblast.lst[[group]][[i]])
    nbr.row[[group]] <- nbr.row[[group]] + nrow(temp)
  }
}


percent <- formatC(unlist(nbr.row / total.reads.df.1) * 100, digits = 2, format = "f")

options(scipen = 1000000)
axis.names <- c(
  paste("G1 \n n = ", nbr.row[["G1"]], "\n", formatC(unlist(percent[1]), digits = 2, format = "f"),   "% of", "\n", total.reads.df.1$G1, " reads", sep = ""),
  paste("G2 \n n = ", nbr.row[["G2"]], "\n", formatC(unlist(percent[2]), digits = 2, format = "f"),   "% of", "\n", total.reads.df.1$G2, " reads", sep = ""),
  paste("G3 \n n = ", nbr.row[["G3"]], "\n", formatC(unlist(percent[3]), digits = 2, format = "f"),   "% of", "\n", total.reads.df.1$G3, " reads", sep = ""),
  paste("GM1 \n n = ", nbr.row[["GM1"]], "\n", formatC(unlist(percent[4]), digits = 2, format = "f"), "% of", "\n", total.reads.df.1$GM1, " reads", sep = ""),
  paste("GM2 \n n = ", nbr.row[["GM2"]], "\n", formatC(unlist(percent[5]), digits = 2, format = "f"), "% of", "\n", total.reads.df.1$GM2, " reads", sep = ""),
  paste("All \n n = ", nbr.row[["All"]], "\n", formatC(unlist(percent[6]), digits = 2, format = "f"), "% of", "\n", total.reads.df.1$All, " reads", sep = "")
)


#Generate graph of the length of the CDR3 region
longueur.cdr3.lst <- list()
for (group in names(igblast.lst)) {
  temp <- igblast.lst[[group]]
  for (i in 1:length(temp)) {
    temp[i]
    longueur.cdr3.lst[[group]] <- append(
      longueur.cdr3.lst[[group]],
      sapply(temp[i], function(x) x$cdr3_aa,
        simplify = TRUE, USE.NAMES = TRUE
      )
    )
  }
}


toute.cdr3 <- list()
for (i in 1:length(longueur.cdr3.lst)) {
  toute.cdr3[[i]] <- data.frame(longueur = longueur.cdr3.lst[[i]], class = names(longueur.cdr3.lst[i]))
}
names(toute.cdr3) <- c("G1", "G2", "G3", "GM1", "GM2", "All")
longueur.cdr3.df <- do.call(rbind.data.frame, toute.cdr3)


#Generate graph of the cdr3 region of more or less than 40aa
less.40 <- list()
more.40 <- list()
less.40.v <- list()
more.40.v <- list()
moss.40.v <- list()
for (group in names(igblast.lst)) {
  pdf(paste("graph/moss_", group, ".pdf", sep = ""), 15, 15)

  for (i in 1:length(igblast.lst[[group]])) {
    temp <- na.omit(igblast.lst[[group]][[i]])
    less.40[[group]] <- rbind(as.data.frame(less.40[[group]]), na.omit(temp[temp$cdr3_aa < 40, ]))
    more.40[[group]] <- rbind(as.data.frame(more.40[[group]]), na.omit(temp[temp$cdr3_aa >= 40, ]))
  }
  more.40[[group]] %>%
    dplyr::count(v_call) %>%
    mutate(perc = (n / nbr.row[[group]] * 100)) -> more.40.v[[group]]
  less.40[[group]] %>%
    dplyr::count(v_call) %>%
    mutate(perc = (n / nbr.row[[group]] * 100)) -> less.40.v[[group]]

  for (V in IGHV_possible) {
    if (!(V %in% more.40.v[[group]]$v_call)) {
      more.40.v[[group]] <- rbind(more.40.v[[group]], c(v_call = V, n = 0, perc = 0))
    }
    if (!(V %in% less.40.v[[group]]$v_call)) {
      less.40.v[[group]] <- rbind(less.40.v[[group]], c(v_call = V, n = 0, perc = 0))
    }
  }

  less.40.v[[group]]$perc <- as.numeric(less.40.v[[group]]$perc)
  less.40.v[[group]]$n <- as.numeric(less.40.v[[group]]$n)
  more.40.v[[group]]$perc <- as.numeric(more.40.v[[group]]$perc)
  more.40.v[[group]]$n <- as.numeric(more.40.v[[group]]$n)

  less.40.v[[group]] <- cbind(less.40.v[[group]], moss = rep("< 40aa", n = nrow(less.40.v)))
  more.40.v[[group]] <- cbind(more.40.v[[group]], moss = rep("> 40aa", n = nrow(more.40.v)))

  less.40.v[[group]] <- less.40.v[[group]][order(less.40.v[[group]]$v_call), ]
  more.40.v[[group]] <- more.40.v[[group]][order(more.40.v[[group]]$v_call), ]

  less.40.v[[group]]$position <- less.40.v[[group]]$perc / 2
  more.40.v[[group]]$position <- less.40.v[[group]]$perc + more.40.v[[group]]$perc / 2

  too.close <- abs(less.40.v[[group]]$position) < 0.5

  more.40.v[[group]][too.close, ]$position <- 1

  moss.40.v[[group]] <- rbind(more.40.v[[group]], less.40.v[[group]])
  moss.40.v[[group]]$moss <- factor(moss.40.v[[group]]$moss, levels = c("> 40aa", "< 40aa"))
  moss.40.v[[group]]$perc <- as.numeric(moss.40.v[[group]]$perc)

  labels <- paste(moss.40.v[[group]]$n, "\n", formatC(moss.40.v[[group]]$perc, digits = 2, format = "f"), "%", sep = "")
  moss.40.v[[group]]$label <- labels
  labels <- paste(less.40.v[[group]]$n, "\n", formatC(less.40.v[[group]]$perc, digits = 2, format = "f"), "%", sep = "")
  less.40.v[[group]]$label <- labels
  labels <- paste(more.40.v[[group]]$n, "\n", formatC(more.40.v[[group]]$perc, digits = 2, format = "f"), "%", sep = "")
  more.40.v[[group]]$label <- labels

  plt <- ggplot(data = moss.40.v[[group]], aes(x = v_call, y = perc, fill = moss)) +
    geom_bar(stat = "identity") +
    theme_light() +
    labs(
      title = paste("IGHV distribution of ", group, sep = ""),
      y = "Reads (%)",
      x = paste("IGHV (N = ", nbr.row[[group]], ")", sep = ""),
      fill = "Length of cdr3"
    ) +
    scale_fill_manual(values = c("lightblue3", "coral2")) +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(vjust = 0, size = 15),
      axis.title.y = element_text(vjust = 2, size = 15),
      axis.text = element_text(color = "black", face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 14)
    ) +
    geom_label(aes(y = position, label = label), color = "black", fill = c(rep("lightblue2", 12), rep("lightcoral", 12)))
  print(plt)

  plt <- ggplot(data = less.40.v[[group]], aes(x = v_call, y = perc)) +
    geom_bar(stat = "identity") +
    theme_light() +
    labs(
      title = paste("IGHV distribution of ", group, " with cdr3 <40 aa", sep = ""),
      y = "Reads (%)",
      x = paste("IGHV ", "\n", "N = ", nbr.row[[group]], "\n", "total nbr of cdr3 < 40 aa = ", sum(less.40.v[[group]]$n), sep = "")
    ) +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(vjust = 0, size = 15),
      axis.title.y = element_text(vjust = 2, size = 15),
      axis.text = element_text(color = "black", face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 10)
    ) +
    geom_label(aes(y = perc / 2, label = label), color = "black")
  print(plt)

  plt <- ggplot(data = more.40.v[[group]], aes(x = v_call, y = perc)) +
    geom_bar(stat = "identity") +
    theme_light() +
    labs(
      title = paste("IGHV distribution of ", group, " with cdr3 >40 aa", sep = ""),
      y = "Reads (%)",
      x = paste("IGHV ", "\n", "N = ", nbr.row[[group]], "\n", "total nbr of cdr3 > 40 aa = ", sum(more.40.v[[group]]$n), sep = "")
    ) +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(vjust = 0, size = 15),
      axis.title.y = element_text(vjust = 2, size = 15),
      axis.text = element_text(color = "black", face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 10)
    ) +
    geom_label(aes(y = perc / 2, label = label), color = "black")
  print(plt)
  dev.off()
 
}
pdf("graph/length.pdf", width = 15, height = 15)

plt <- ggplot(data = longueur.df, aes(x = class, y = longueur)) +
  geom_violin(trim = TRUE, color = "black", bw = 1.75, fill = "grey90") +
  geom_boxplot(width = 0.1, fill = "lightblue") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2) +
  labs(
    title = "Valid sequences length",
    y = "Sequence length (nt)",
    x = "Class"
  ) +
  theme_light() +
  theme(
    title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(vjust = 0, size = 20),
    axis.title.y = element_text(vjust = 2, size = 20),
    axis.text = element_text(color = "black", face = "bold", size = 17),
    axis.text.x = element_text(face = "bold", size = 16)
  ) +
  scale_x_discrete(labels = axis.names)

print(plt)


plt <- ggplot(data = longueur.cdr3.df, aes(x = class, y = longueur)) +
  geom_violin(trim = TRUE, color = "black", bw = 0.75, fill = "grey90") +
  geom_boxplot(width = 0.1, fill = "lightblue") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2) +
  labs(
    title = "Length of CDR3 region",
    y = "Length of region (AA)",
    x = "Class"
  ) +
  theme_light() +
  theme(
    title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(vjust = 0, size = 20),
    axis.title.y = element_text(vjust = 2, size = 20),
    axis.text = element_text(color = "black", face = "bold", size = 17),
    axis.text.x = element_text(face = "bold", size = 16)
  ) +
  scale_x_discrete(labels = axis.names)

print(plt)


longueur.cdr3.df.more <- longueur.cdr3.df[longueur.cdr3.df$longueur >= 40,]
nbr.more <- list()
for(group in c("G1", "G2", "G3", "GM1", "GM2", "All")){
  nbr.more[[group]] = nrow(longueur.cdr3.df.more[longueur.cdr3.df.more$class == group,])
}
percent.more <- list()
for(group in c("G1", "G2", "G3", "GM1", "GM2", "All")){
  percent.more[[group]] = nbr.more[[group]]/total.reads.lst[[group]]
}



axis.names.more <- c(
  paste("G1 \n n = ",  nbr.more[["G1"]],  "\n",  formatC(unlist(percent.more[["G1"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$G1,  " reads", sep = ""),
  paste("G2 \n n = ",  nbr.more[["G2"]],  "\n",  formatC(unlist(percent.more[["G2"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$G2,  " reads", sep = ""),
  paste("G3 \n n = ",  nbr.more[["G3"]],  "\n",  formatC(unlist(percent.more[["G3"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$G3,  " reads", sep = ""),
  paste("GM1 \n n = ", nbr.more[["GM1"]], "\n", formatC(unlist(percent.more[["GM1"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$GM1, " reads", sep = ""),
  paste("GM2 \n n = ", nbr.more[["GM2"]], "\n", formatC(unlist(percent.more[["GM2"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$GM2, " reads", sep = ""),
  paste("All \n n = ", nbr.more[["All"]], "\n", formatC(unlist(percent.more[["All"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$All, " reads", sep = "")
)

plt <- ggplot(data = longueur.cdr3.df.more, aes(x = class, y = longueur))+
  geom_violin(trim = TRUE, color = "black", bw = 1.75, fill = "grey90")+
  geom_boxplot(width = 0.1, fill = "lightblue") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2) +
  labs(
    title = "Valid sequences length with cdr3 >= 40 AA",
    y = "Sequence length (nt)",
    x = "Class")+
  theme_light() +
  theme(title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 2, size = 20),
        axis.text = element_text(color = "black", face = "bold", size = 17),
        axis.text.x = element_text(face = "bold", size = 16)) +
  scale_x_discrete(labels = axis.names.more)

print(plt)


longueur.cdr3.df.less <- longueur.cdr3.df[longueur.cdr3.df$longueur < 40,]
nbr.less <- list()
for(group in c("G1", "G2", "G3", "GM1", "GM2", "All")){
  nbr.less[[group]] = nrow(longueur.cdr3.df.less[longueur.cdr3.df.less$class == group,])
}
percent.less <- list()
for(group in c("G1", "G2", "G3", "GM1", "GM2", "All")){
  percent.less[[group]] = nrow(nbr.less[[group]])/total.reads.lst[[group]]
}

axis.names.less <- c(
  paste("G1 \n n = ",  nbr.less[["G1"]],  "\n",  formatC(unlist(percent.more[["G1"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$G1,  " reads", sep = ""),
  paste("G2 \n n = ",  nbr.less[["G2"]],  "\n",  formatC(unlist(percent.more[["G2"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$G2,  " reads", sep = ""),
  paste("G3 \n n = ",  nbr.less[["G3"]],  "\n",  formatC(unlist(percent.more[["G3"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$G3,  " reads", sep = ""),
  paste("GM1 \n n = ", nbr.less[["GM1"]], "\n", formatC(unlist(percent.more[["GM1"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$GM1, " reads", sep = ""),
  paste("GM2 \n n = ", nbr.less[["GM2"]], "\n", formatC(unlist(percent.more[["GM2"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$GM2, " reads", sep = ""),
  paste("All \n n = ", nbr.less[["All"]], "\n", formatC(unlist(percent.more[["All"]]), digits = 4, format = "f"),  "% of", "\n", total.reads.df.1$All, " reads", sep = "")
)



plt <- ggplot(data = longueur.cdr3.df.less, aes(x = class, y = longueur))+
  geom_violin(trim = TRUE, color = "black", bw = 1.75, fill = "grey90") +
  geom_boxplot(width = 0.1, fill = "lightblue") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2) +
  labs(
    title = "Valid sequences length with cdr3 < 40 AA",
    y = "Sequence length (AA)",
    x = "Class") +
  theme_light() +
  theme(title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 2, size = 20),
        axis.text = element_text(color = "black", face = "bold", size = 17),
        axis.text.x = element_text(face = "bold", size = 16)) +
  scale_x_discrete(labels = axis.names.less)

print(plt)




#### Use a lot of memory
## NGmerge
#chemins_ngmerge <- sapply(id_all, function(x) paste("./data/mergedReads/", x, ".log", sep =""),
#                          simplify = FALSE, USE.NAMES = TRUE)
#ngmerge <-  sapply(chemins_ngmerge, function(x) read.csv(x, sep= "\t", na.strings=c("NA")),
#                   simplify = FALSE, USE.NAMES = TRUE)

#print("NGmerge loaded")

# Generate graph of verlap of reads in NGMerge

#overlap <- sapply(names(ngmerge), function(x) ngmerge[[x]]$OverlapLen,
#                  simplify = FALSE, USE.NAMES = TRUE)
#data_overlap <- data.frame(longueur_overlap = NA, nom = NA)
#temp.df <- sapply(names(ngmerge),
#                  function(x) data.frame(longueur_overlap = overlap[[x]],
#                                         nom = rep(str_extract(x, "G[123M][12]?"), length(overlap[[x]]))),
#                  simplify = FALSE, USE.NAMES = TRUE)
#
#
#for(i in (1:length(temp.df))) data_overlap <- rbind(data_overlap, temp.df[[i]])

#temp_plot <- ggplot(data = data_overlap[!is.na(data_overlap$nom),], mapping = aes(x = nom, y = longueur_overlap))+
#  geom_violin(trim = TRUE, color = "black", fill = "grey90") +
#  geom_boxplot(width = 0.1, fill = "lightblue") +
#  stat_summary(fun=mean, color = "black", geom="point", shape=23, size=2)+
#  labs(title = "Longueur de l'overlap lors de NGMerge",
#       y = "Longueur de l'overlap (nt)",
#       x = "Echantillon")+
#  theme_light()+
#  theme(title = element_text(size = 12, face = 'bold'),
#        axis.title.x = element_text(vjust = 0, size = 15),
#        axis.title.y = element_text(vjust = 2, size = 15),
#        axis.text    = element_text(color = "black", face = "bold", size = 14),
#        axis.text.x  = element_text(face = "bold", size = 13))

#print(temp_plot)
dev.off()


sessionInfo()
