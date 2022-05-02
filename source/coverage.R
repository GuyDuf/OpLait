
#### Chargement des librairies ####

packages <- c("ggplot2", "viridis", "ShortRead", "stringr", "tidyr", "tibble")

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

print("Package loaded!")

args <- commandArgs(trailingOnly = TRUE)
# For testing
# args <- c("7-603-IgG2-1_S139","7-603-IgG2-1_S139","7-603-IgG3-1_S7","7-603-IgG3-1_S7","7-603-IgGM1-1_S51","28-272-IgGM2-1_S204","3-241-IgGM2-1_S179","32-234-IgG1-1_S120","33-241-IgGM1-1_S77","35-603-IgG1-1_S123","36-618-IgGM2-1_S212","39-279-IgG1-1_S127","40-281-IgG3-1_S40","42-637-IgGM2-1_S218","5-253-IgG1-1_S93","6-253-IgG2-1_S138","7-603-IgGM1-1_S51")

id_all <- sort(unique(args))

print(args)
print(id_all)

groupement <- list(
  "G1" = grep(pattern = "IgG1", x = id_all),
  "G2" = grep(pattern = "IgG2", x = id_all),
  "G3" = grep(pattern = "IgG3", x = id_all),
  "GM1" = grep(pattern = "IgGM1", x = id_all),
  "GM2" = grep(pattern = "IgGM2", x = id_all),
  "ALL" = grep(pattern = "Ig", x = id_all)
)

#### Début de la loop ####
for (group in names(groupement)) {
  pdf(paste("graph/coverage", group, ".pdf", sep = ""), width = 10, height = 10)

  id <- id_all[unlist(groupement[group])]
  nbr_echantillon <- length(id)
  print(id)
  nom_echantillon <- unlist(id)
  names(id) <- nom_echantillon # permet de garder le nom de lorsque j'utilise sapply

  IGHV_possible <- c("IGHV1-7", "IGHV1-10", "IGHV1-14", "IGHV1-17", "IGHV1-21/33", "IGHV1-25", "IGHV1-27", "IGHV1-30", "IGHV1-37", "IGHV1-39", "IGHV1-20", "IGHV1-32")
  IGHD_possible <- c("IGHD1-1", "IGHD1-2/4", "IGHD1-3", "IGHD1-4", "IGHD2-1/2/3/4", "IGHD3-1/3/4", "IGHD4-1", "IGHD5-2", "IGHD5-3/4", "IGHD6-2", "IGHD6-3/4", "IGHD7-3", "IGHD7-4", "IGHD8-2", "IGHD9-1/4")
  IGHJ_possible <- c("IGHJ1-4", "IGHJ1-6", "IGHJ2-4")
  print(43)

  # Fonctions #

  name_clean_up <- function(string) {
    temp <- substring(string, first = 1, last = (str_locate(string = string, pattern = ","))[1] - 1)
    if (is.na(temp)) {
      string
    } else {
      temp
    }
  }

  # Chargement des fichiers necessaires #
  ## IgBLAST
  chemins_igblast <- sapply(nom_echantillon, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv2", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )
  igblast <- sapply(chemins_igblast, function(x) read.csv(x, sep = "\t", na.strings = c("", "NA")),
    simplify = FALSE, USE.NAMES = TRUE
  )


  nbr_ligne_igblast_raw <- list(rep(1, times = nbr_echantillon))

  keep <- c("sequence_alignment", "stop_codon", "v_frameshift", "vj_in_frame", "complete_vdj")

  print(67)
  for (i in (1:length(igblast))) {
    igblast[[i]] <- igblast[[i]][, (colnames(igblast[[i]]) %in% keep)]
    nbr_ligne_igblast_raw[[i]] <- nrow(igblast[[i]])
  }

  print(72)
  names(nbr_ligne_igblast_raw) <- names(igblast)

  for (j in names(igblast)) {
    for (i in (1:length(igblast[[j]]))) {
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$stop_codon), ]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$stop_codon == FALSE, ]
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$v_frameshift), ]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$v_frameshift == FALSE, ]
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$vj_in_frame), ]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$vj_in_frame == TRUE, ]
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$complete), ]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$complete_vdj == TRUE, ]
    }
  }


  print("igblast loaded")

  ## Raw reads
  chemin_raw_reads <- sapply(id, function(x) paste("./data/rawReads/", x, "_L001_R1_001.fastq.gz", sep = ""),
    simplify = FALSE, USE.NAMES = TRUE
  )

  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
    simplify = FALSE, USE.NAMES = TRUE
  )

  #### Pour les ajouts provenant de coverage.R
  fastqs.filter <- list(rep(1, times = nbr_echantillon))


  for (j in 1:length(igblast)) {
    fastqs.filter[[j]] <- igblast[[j]]$sequence
  }

  # Initialisation du PDF
  #### Coverage ####
  ## Estimate the number of new sequences we get as we increase the number of sequenced reads
  ## the idea is to know at with depth we will reach saturation
  n.seqs <- c(rep(1, times = length(id)))
  for (n in 1:length(id)) {
    n.seqs[n] <- nrow(igblast[[n]])
    names(n.seqs) <- names(igblast)
  }

  print(n.seqs)

  coverageSeqPerc <- function(fastq) {
    perc.seq <- seq(0.05, 1, 0.05)
    N <- length(fastq)
    rand.n <- 10
    to.ret <- sapply(perc.seq, function(pi) {
      cur.n <- floor(pi * N)
      as.numeric(lapply(1:rand.n, function(ri) {
        length(unique(sread(fastq[sample(1:N, cur.n)])))
      }))
    })
    colnames(to.ret) <- sprintf("%.2f", perc.seq)
    to.ret
  }

  coverageSeqNum <- function(fastq, n.seq.rnd = seq(1000, 5000, 200)) {
    N <- length(fastq)
    rand.n <- 10
    to.ret <- sapply(n.seq.rnd, function(ni) {
      cur.n <- ni
      as.numeric(lapply(1:rand.n, function(ri) {
        length(unique(fastq[sample(1:N, cur.n)]))
      }))
    })
    colnames(to.ret) <- sprintf("%d", n.seq.rnd)
    to.ret
  }



  ### added
  fastqs.filter <- list(rep(1, length(igblast)))
  for (j in 1:length(igblast)) {
    fastqs.filter[[j]] <- igblast[[j]]$sequence
  }


  ###


  set.seed(1)
  n.seq.rnd <- seq(100, min(n.seqs), 200)

  print(643)
  cov.fastq <- lapply(fastqs.filter, function(x) {
    coverageSeqNum(x, n.seq.rnd)
  })

  cov.mat <- sapply(cov.fastq, colMeans)
  matplot(n.seq.rnd, cov.mat, lty = 1, lwd = 4, type = "l", col = 1:length(id), main = group, ylab = "Number of unique sequences", xlab = "number of sequences")
  abline(a = 0, b = 1, lty = 2, col = "gray")

  seq.depth <- seq(0, 200000, 1000)

  mult.total.reads <- median(unlist(nbr_raw_reads) / n.seqs)
  la.list <- c("empty")
  par(mfrow = c(2, 2))
  for (i in 1:ncol(cov.mat)) {
    cur.fit <- tryCatch(nls(y ~ RMAX * (1 - exp(-k * x)),
      data = list(x = n.seq.rnd, y = cov.mat[, i]),
      start = list(RMAX = 50000, k = 0.002)
    ),
    error = function(e) {
      print(e)
      return(names(igblast[i]))
    }
    )
    if (cur.fit != names(igblast[i])) {
      plot(n.seq.rnd, cov.mat[, i], main = paste("Fit ", names(igblast[i]), sep = ""), ylab = "Num unique reads")
      abline(a = 0, b = 1, lty = 2, col = "gray")
      lines(n.seq.rnd, predict(cur.fit, list(x = n.seq.rnd)))
      plot(seq.depth * mult.total.reads, predict(cur.fit, list(x = seq.depth)) * 100 / summary(cur.fit)$coefficients[1, 1],
        type = "l", lwd = 3, xlab = "Number sequenced reads [corrected]", ylab = "Percentage of total reads"
      )
    } else {
      la.list <- append(la.list, cur.fit)
    }
  }



  print("coverage fin")
  dev.off()
}


sessionInfo()
