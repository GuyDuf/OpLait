
### permet d'installer les packages necessaires

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