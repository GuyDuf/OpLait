setwd("/media/guyllaume/BackUp/PaquetLab/Batch2/data/rawReads/")
namesOfFile <- list.files(path = ".")
View(namesOfFile)
namesOfFile.df <- as.data.frame(strsplit(namesOfFile, split = "-"))
View(namesOfFile.df)
colnames(namesOfFile.df) <- c(1:ncol(namesOfFile.df))
namesOfFile.df <- t(namesOfFile.df)
name_combo <- table(namesOfFile.df[,3],namesOfFile.df[,2])
write.csv(name_combo, file = "freqNom.csv")
View(name_combo)
