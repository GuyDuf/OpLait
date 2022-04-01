
#### Chargement des librairies ####

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

args <- commandArgs(trailingOnly=TRUE)
args <- c("26-258-IgGM2-1_S202","1-234-IgG1-1_S89","1-234-IgG1-1_S89","1-234-IgG2-1_S133","1-234-IgG2-1_S133","1-234-IgG3-1_S1","1-234-IgG3-1_S1","1-234-IgGM1-1_S45","1-234-IgGM1-1_S45","1-234-IgGM2-1_S177","1-234-IgGM2-1_S177","10-618-IgG1-1_S98","10-618-IgG1-1_S98","10-618-IgG2-1_S142","10-618-IgG2-1_S142","10-618-IgG3-1_S10","10-618-IgG3-1_S10","10-618-IgGM1-1_S54","10-618-IgGM1-1_S54","10-618-IgGM2-1_S186","30-275-IgG1-1_S118","30-275-IgG1-1_S118","30-275-IgG2-1_S162","30-275-IgG2-1_S162","30-275-IgG3-1_S30","30-275-IgG3-1_S30","30-275-IgGM1-1_S74","30-275-IgGM1-1_S74","30-275-IgGM2-1_S206","30-275-IgGM2-1_S206","31-277-IgG1-1_S119","31-277-IgG1-1_S119","31-277-IgG2-1_S163","31-277-IgG2-1_S163","31-277-IgG3-1_S31","31-277-IgG3-1_S31","31-277-IgGM1-1_S75","31-277-IgGM1-1_S75","31-277-IgGM2-1_S207","31-277-IgGM2-1_S207","40-281-IgG3-1_S40","40-281-IgGM1-1_S84","40-281-IgGM1-1_S84","40-281-IgGM2-1_S216","40-281-IgGM2-1_S216","41-291-IgG1-1_S129","41-291-IgG1-1_S129","41-291-IgG2-1_S173","41-291-IgG2-1_S173","41-291-IgG3-1_S41","41-291-IgG3-1_S41","41-291-IgGM1-1_S85","41-291-IgGM1-1_S85","41-291-IgGM2-1_S217","41-291-IgGM2-1_S217","42-637-IgG1-1_S130","42-637-IgG1-1_S130","42-637-IgG2-1_S174","42-637-IgG2-1_S174","42-637-IgG3-1_S42","42-637-IgG3-1_S42","42-637-IgGM1-1_S86","42-637-IgGM1-1_S86","42-637-IgGM2-1_S218","10-618-IgGM2-1_S186","12-9004-IgGM1-1_S56","14-281-IgGM2-1_S190","17-641-IgG1-1_S105","18-720-IgG3-1_S18","19-257-IgGM2-1_S195","20-261-IgG1-1_S108","21-262-IgG3-1_S21","22-264-IgGM1-1_S66","24-268-IgG1-1_S112","25-271-IgGM1-1_S69","21-262-IgG3-1_S21","21-262-IgGM1-1_S65","21-262-IgGM1-1_S65","21-262-IgGM2-1_S197","21-262-IgGM2-1_S197","22-264-IgG1-1_S110","22-264-IgG1-1_S110","22-264-IgG2-1_S154","22-264-IgG2-1_S154","22-264-IgG3-1_S22","22-264-IgG3-1_S22","22-264-IgGM1-1_S66","37-278-IgG1-1_S125","37-278-IgG1-1_S125","37-278-IgG2-1_S169","37-278-IgG2-1_S169","37-278-IgG3-1_S37","37-278-IgG3-1_S37","37-278-IgGM1-1_S81","37-278-IgGM1-1_S81","37-278-IgGM2-1_S213","37-278-IgGM2-1_S213","38-9004-IgG1-1_S126","38-9004-IgG1-1_S126","38-9004-IgG2-1_S170","38-9004-IgG2-1_S170","38-9004-IgG3-1_S38","38-9004-IgG3-1_S38","38-9004-IgGM1-1_S82","38-9004-IgGM1-1_S82","38-9004-IgGM2-1_S214","38-9004-IgGM2-1_S214","5-253-IgG1-1_S93","5-253-IgG2-1_S137","5-253-IgG2-1_S137","5-253-IgG3-1_S5","5-253-IgG3-1_S5","5-253-IgGM1-1_S49","5-253-IgGM1-1_S49","5-253-IgGM2-1_S181","5-253-IgGM2-1_S181","6-253-IgG1-1_S94","6-253-IgG1-1_S94","6-253-IgG2-1_S138","17-641-IgG2-1_S149","17-641-IgG2-1_S149","17-641-IgG3-1_S17","17-641-IgG3-1_S17","17-641-IgGM1-1_S61","17-641-IgGM1-1_S61","17-641-IgGM2-1_S193","17-641-IgGM2-1_S193","18-720-IgG1-1_S106","18-720-IgG1-1_S106","18-720-IgG2-1_S150","18-720-IgG2-1_S150","18-720-IgG3-1_S18","27-266-IgG1-1_S115","27-266-IgG1-1_S115","27-266-IgG2-1_S159","27-266-IgG2-1_S159","27-266-IgG3-1_S27","27-266-IgG3-1_S27","27-266-IgGM1-1_S71","27-266-IgGM1-1_S71","27-266-IgGM2-1_S203","27-266-IgGM2-1_S203","28-272-IgG1-1_S116","28-272-IgG1-1_S116","28-272-IgG2-1_S160","28-272-IgG2-1_S160","28-272-IgG3-1_S28","28-272-IgG3-1_S28","28-272-IgGM1-1_S72","28-272-IgGM1-1_S72","28-272-IgGM2-1_S204")
id_all <-  sort(unique(args))
nbr_echantillon = 10

file_nbr = seq(from = 1, to = length(id_all), by = nbr_echantillon)
print(file_nbr)

file_nbr =  file_nbr[1:length(file_nbr)]

print(args)
print(id_all)

groupement <- list("G1" = grep(pattern = "IgG1", x = id_all),
                   "G2" = grep(pattern = "IgG2", x = id_all),
                   "G3" = grep(pattern = "IgG3", x = id_all),
                   "GM1" = grep(pattern = "IgGM1", x = id_all),
                   "GM2" = grep(pattern = "IgGM2", x = id_all))

#### Début de la loop ####
for(groupe in groupement){
  
  #print(paste("Génération des graphiques pour les séquences ", y,"-", y + nbr_echantillon-1, sep =""))
  
  id = id_all[unlist(groupe)]
  print(id)
  nom_echantillon  <- unlist(id)
  names(id) <- nom_echantillon # permet de garder le nom de lorsque j'utilise sapply
  
  IGHV_possible = c("IGHV1-7","IGHV1-10","IGHV1-14","IGHV1-17","IGHV1-21/33","IGHV1-25","IGHV1-27","IGHV1-30","IGHV1-37","IGHV1-39","IGHV1-20","IGHV1-32")
  IGHD_possible = c("IGHD1-1","IGHD1-2/4","IGHD1-3","IGHD1-4","IGHD2-1/2/3/4","IGHD3-1/3/4","IGHD4-1","IGHD5-2","IGHD5-3/4","IGHD6-2","IGHD6-3/4","IGHD7-3","IGHD7-4","IGHD8-2","IGHD9-1/4")
  IGHJ_possible = c("IGHJ1-4","IGHJ1-6","IGHJ2-4")
  print(43)
  
  # Fonctions #
  
  name_clean_up <- function(string){
    temp <- substring(string, first=1 ,last = (str_locate(string = string, pattern = ","))[1] -1)
    if(is.na(temp)){
      string
    } else {
      temp
    }
  }
  
  # Chargement des fichiers necessaires # 
  ## IgBLAST
  chemins_igblast <- sapply(nom_echantillon, function(x) paste("./output/VDJ_", x, ".csv_dropped.csv2", sep = ""),
                            simplify = FALSE, USE.NAMES = TRUE)
  igblast <- sapply(chemins_igblast, function(x) read.csv(x, sep= "\t", na.strings=c("","NA")),
                    simplify = FALSE, USE.NAMES = TRUE)
  
  
  nbr_ligne_igblast_raw <- list(rep(1, times = nbr_echantillon))
  
  keep = c("sequence_alignment", "stop_codon","v_frameshift","vj_in_frame",'complete_vdj')
  
  print(67)
  for(i in (1:length(igblast))){
    igblast[[i]] <- igblast[[i]][,(colnames(igblast[[i]]) %in% keep)]
    nbr_ligne_igblast_raw[[i]] <- nrow(igblast[[i]])
  }
  
  print(72)
  names(nbr_ligne_igblast_raw) <- names(igblast)
  
  for(j in names(igblast)){
    for(i in (1:length(igblast[[j]]))){
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$stop_codon),]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$stop_codon == FALSE,]
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$v_frameshift),]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$v_frameshift == FALSE,]
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$vj_in_frame),]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$vj_in_frame == TRUE,]
      igblast[[j]] <- igblast[[j]][!is.na(igblast[[j]]$complete),]
      igblast[[j]] <- igblast[[j]][igblast[[j]]$complete_vdj  == TRUE,]
    }
  }
  
  
  print("igblast loaded")
  
  ## Raw reads
  chemin_raw_reads <- sapply(id, function(x) paste("./data/rawReads/", x,"_L001_R1_001.fastq.gz", sep =""),
                             simplify = FALSE, USE.NAMES = TRUE)
  
  nbr_raw_reads <- sapply(chemin_raw_reads, function(x) length(readFastq(x)),
                          simplify = FALSE, USE.NAMES = TRUE)
  
  #### Pour les ajouts provenant de coverage.R
  fastqs.filter <- list(rep(1,times = nbr_echantillon))
  
  
  for(j in 1:length(igblast)){
    fastqs.filter[[j]] <- igblast[[j]]$sequence
  }
  
  # Initialisation du PDF
  #### Coverage ####
  ## Estimate the number of new sequences we get as we increase the number of sequenced reads
  ## the idea is to know at with depth we will reach saturation
  n.seqs <- c(rep(1, times = length(id)))
  for(n in 1:length(id)){
    n.seqs[n] <- nrow(igblast[[n]])
    names(n.seqs) <- names(igblast)
  }
  
  print(n.seqs)
  
  coverageSeqPerc <- function(fastq){
    perc.seq <- seq(0.05,1,0.05)
    N <- length(fastq)
    rand.n <- 10
    to.ret <- sapply(perc.seq,function(pi){
      cur.n <- floor(pi*N)
      as.numeric(lapply(1:rand.n,function(ri){
        length(unique(sread(fastq[sample(1:N,cur.n)])))
      }))
    })
    colnames(to.ret) <- sprintf("%.2f",perc.seq)
    to.ret
  }
  
  coverageSeqNum <- function(fastq,n.seq.rnd = seq(1000,5000,200)){
    N <- length(fastq)
    rand.n <- 10
    to.ret <- sapply(n.seq.rnd,function(ni){
      cur.n <- ni
      as.numeric(lapply(1:rand.n,function(ri){
        length(unique(fastq[sample(1:N,cur.n)])) 
      }))
    })
    colnames(to.ret) <- sprintf("%d",n.seq.rnd)
    to.ret
  }
  
  #fastqs
  
  #fastqs.filter <- lapply(fastqs,function(fsi){
  #  fsi[nchar(sread(fsi)) >= 400]
  #})
  
  pdf(paste("graph/coverage",y,"-", y+9,".pdf",sep=""),width=10, height=10)
  
  ### added
  fastqs.filter <- list(rep(1,length(igblast)))
  for(j in 1:length(igblast)){
    fastqs.filter[[j]] <- igblast[[j]]$sequence
  }
  
  
  ###
  
  #print(fastqs.filter)
  
  set.seed(1)
  n.seq.rnd <- seq(100,min(n.seqs),200)
  
  print(643)
  cov.fastq <- lapply(fastqs.filter,function(x){
    coverageSeqNum(x,n.seq.rnd)
  })#,mc.cores=5)
  
  cov.mat <- sapply(cov.fastq,colMeans)
  matplot(n.seq.rnd,cov.mat,lty=1,lwd=4,type="l",col = 1:length(id),ylab="Number of unique sequences",xlab="number of sequences")
  abline(a=0,b=1,lty=2,col="gray")
  
  
  legend("topleft",id,lwd=3,col=1:length(id))
  
  seq.depth = seq(0,200000,1000)
  
  mult.total.reads = median(unlist(nbr_raw_reads)/n.seqs)
  la.list = c("empty")
  par(mfrow=c(4,2))
  for(i in 1:ncol(cov.mat)){
    cur.fit <- tryCatch(nls(y ~ RMAX * (1-exp(-k*x)),
                            data=list(x=n.seq.rnd,y=cov.mat[,i]),
                            start=list(RMAX=50000,k=0.002)),
                        error = function(e) {
                          print(e)
                          return(names(igblast[i]))})
    if(cur.fit != names(igblast[i])){
      plot(n.seq.rnd,cov.mat[,i],main="Fit",ylab="Num unique reads")
      #print(plt)
      abline(a=0,b=1,lty=2,col="gray")
      lines(n.seq.rnd,predict(cur.fit,list(x=n.seq.rnd)))
      plot(seq.depth*mult.total.reads,predict(cur.fit,list(x=seq.depth))*100/summary(cur.fit)$coefficients[1,1],
           type="l",lwd=3,xlab="Number sequenced reads [corrected]",ylab="Percentage of total reads")
    }else{
      la.list <- append(la.list,cur.fit)
    }
  }
  
  
  
  print('coverage fin')
  dev.off()
}


sessionInfo()

