
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
    library(x, character.only = TRUE)
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


### CDR3 GRAPHS ####

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
  print(231)
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
  print(167)
}

sessionInfo()
