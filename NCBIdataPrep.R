# Created by Javad Rahimipour Anaraki on 14/04/19
# Ph.D. Candidate
# Department of Computer Science
# Memorial University of Newfoundland
# jra066 [AT] mun [DOT] ca | www.cs.mun.ca/~jra066

#   input: Data from NCBI data browser
#  output: Four CSV files, data with and without features, features names and corresponding chromosomes

rm(list=ls())
#========================Libraries=========================
library(data.table)
library(dplyr)
library(rstudioapi)

#Getting current folder path
path <- dirname(rstudioapi::getSourceEditorContext()$path)
fileName <- dir(path, '*.soft')

#Find the beginning of the data
linesToSkip <- system(paste0('grep -nr !dataset_table_begin ', paste0(path, '/', fileName)), intern = T)
linesToSkip <- as.numeric(strsplit(linesToSkip, ':')[[1]][1])

#Extract classes and update labels
nClasses <- system(paste0('grep -nr !subset_description ', paste0(path, '/', fileName)), intern = T)
samples <- system(paste0('grep -nr !subset_sample_id ', paste0(path, '/', fileName)), intern = T)
class <- NULL
for (i in 1:length(nClasses)){
  class <- rbind(class, cbind(strsplit(strsplit(samples[i], '= ')[[1]][2], ',')[[1]], i-1))
}
class <- as.data.frame(class)

#Reading the data file and class in
data <- fread(paste0(path, '/', fileName), sep = "\t", header = F, quote = "", skip = linesToSkip, data.table = F)

#Extract features names
cols <- paste(data[2:nrow(data), 1], data[2:nrow(data), 2], sep = "_")

#Extract chromosomes data
idx <- grep("Chromosome annotation", data[1, ])
tmp <- strsplit(data[2:nrow(data), idx], ",")
chromosomes <- NULL
for (i in 1:length(tmp)){
  chromosomes[i] <- tmp[[i]][1]
}

#Choose the required portion of the data
idx <- grep("Gene title", data[1, ])
data <- data[2:nrow(data), 3:(idx[1]+2)]

#Transpose, add column names and labels to the data
data <- as.data.frame(t(data))
colnames(data) <- cols
data$class <- class$V2

#Shuffle the data
data <- data[sample(nrow(data)), ]

#Store data with and without features, features names and corresponding chromosomes
paste0(path, '/', fileName)
write.table(data, file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '.csv'), sep = ",", quote = F, row.names = F, col.names = T)
write.table(data, file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_NoFeature.csv'), sep = ",", quote = F, row.names = F, col.names = F)
write.table(cols, file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_Justfeatures.csv'), sep = ",", quote = F, row.names = F, col.names = F)
write.table(chromosomes, file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_JustChromosomes.csv'), sep = ",", quote = F, row.names = F, col.names = F)
