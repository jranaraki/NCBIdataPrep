# Created by Javad Rahimipour Anaraki on 14/04/19
# Ph.D. Candidate
# Department of Computer Science
# Memorial University of Newfoundland
# jra066 [AT] mun [DOT] ca | www.cs.mun.ca/~jra066

#   input: Data from NCBI data browser
#  output: Four CSV files, data with and without features, features names and corresponding chromosomes

rm(list = ls())
#========================Libraries=========================
list.of.packages <-
  c("rstudioapi",
    "data.table",
    "e1071",
    "dplyr")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

library(data.table)
library(dplyr)
library(rstudioapi)
library(e1071)

#Getting current folder path
path <- dirname(rstudioapi::getSourceEditorContext()$path)
fileName <- dir(path, '*.soft')

#Find the beginning of the data
startLine <-
  system(paste0('grep -nr !dataset_table_begin ', paste0(path, '/', fileName)), intern = T)
startLine <- as.numeric(strsplit(startLine, ':')[[1]][1])

#Find the ending of the data
endLine <-
  system(paste0('grep -nr !dataset_table_end ', paste0(path, '/', fileName)), intern = T)
endLine <- as.numeric(strsplit(endLine, ':')[[1]][1]) - 1

#Extract classes and update labels
nClasses <-
  system(paste0('grep -nr !subset_description ', paste0(path, '/', fileName)), intern = T)
samples <-
  system(paste0('grep -nr !subset_sample_id ', paste0(path, '/', fileName)), intern = T)
class <- NULL
for (i in 1:length(nClasses)) {
  class <-
    rbind(class, cbind(strsplit(strsplit(samples[i], '= ')[[1]][2], ',')[[1]], i))
}
class <- as.data.frame(class)
colnames(class) <- c("V1", "V2")
unq <- as.data.frame(unique(class$V1))
colnames(unq) <- "V1"
if (nrow(class) != nrow(unq)) {
  for (i in 1:nrow(unq)) {
    tmp <- filter(class, V1 == unq[i, "V1"])
    unq[i, "V2"] <- paste(tmp$V2, collapse = '')
  }
  class <- unq
  unq <- NULL
}

#Reading the data file and class in
data <-
  fread(
    paste0(path, '/', fileName),
    sep = "\t",
    header = F,
    quote = "",
    skip = startLine,
    nrows = (endLine - startLine),
    data.table = F,
    showProgress = T
  )

#Extract features names
cols <-
  paste(data[2:nrow(data), 1], data[2:nrow(data), 2], sep = "_")

#Extract chromosomes data
idx <- grep("Chromosome annotation", data[1, ])
tmp <- strsplit(data[2:nrow(data), idx], ",")
chromosomes <- NULL
for (i in 1:length(tmp)) {
  chromosomes[i] <- tmp[[i]][1]
}
tmp <- NULL

#Choose the required portion of the data
idx <- grep("GSM", data[1, ])
data <- data[2:nrow(data), idx]

#Transpose, add column names and labels to the data
data <- as.data.frame(t(data))
data <-
  as.data.frame(apply(data, 2, function(x)
    as.numeric(as.character(x))))
colnames(data) <- cols

#Remove fully null columns
nullFeatures <- sapply(data, function(x)
  all(is.na(x)))
nullFeatures <- which(nullFeatures == T)
if (length(nullFeatures) > 0) {
  data <- data[, -nullFeatures]
}

#Impute the data
if (length(which(is.na(data) == T)) > 0) {
  data <- as.data.frame(impute(data))
}

#Shuffle the data
data$class <- class$V2
data <- data[sample(nrow(data)), ]

#Store data with and without features, features names and corresponding chromosomes
write.table(
  data,
  file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '.csv'),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  data,
  file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_NoFeature.csv'),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)

write.table(
  cols,
  file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_Justfeatures.csv'),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)

write.table(
  chromosomes,
  file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_JustChromosomes.csv'),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)