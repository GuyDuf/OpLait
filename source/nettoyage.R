
### élimine les séquences trop courte après NGmerge



### Loading Packages
packages = c("ShortRead")

package.check <- lapply(packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }})


###### LONGUEUR MINIMAL

longueur_minimal = 425

#### Loading Files

args <- commandArgs(trailingOnly=TRUE)
id = unique(args)

#### Pour tester
#id <- list("IgG1-1_S35","IgG2-1_S36","IgG3-1_S37","IgM1-1_S38","IgM2-1_S39")
#args=id



chemins <- sapply(id, function(x) paste("./data/mergedReads/", x,".fastq", sep =""),
                           simplify = FALSE, USE.NAMES = TRUE)
names(chemins) <- id

ensemble <- sapply(names(chemins), function(x) sread(readFastq(chemins[[x]])),
                   simplify = FALSE, USE.NAMES = TRUE)


### Remove Reads < 425
ensemble <- lapply(ensemble, function(x) x[width(x) > longueur_minimal])

### Write Files
lapply(names(ensemble), function(x) {
  writeFasta(ensemble[[x]], file = paste("./output/", x, ".fasta", sep=""))})
