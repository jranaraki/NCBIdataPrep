# Created by Javad Rahimipour Anaraki on 14/04/19
# Ph.D. in Computer Science
# Department of Computer Science
# Memorial University of Newfoundland
# jra066 [AT] mun [DOT] ca | www.cs.mun.ca/~jra066

# Modified by Majid Afshar on 23/08/19
# Ph.D. candidate in Computer Science
# Department of Computer Science
# Memorial University of Newfoundland
# mman23 [AT] mun [DOT] ca | www.cs.mun.ca/~mman23

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
    rbind(class, cbind(strsplit(strsplit(samples[i], '= ')[[1]][2], ',')[[1]], i - 1))
}
class <- as.data.frame(class)

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
genes <- data[2:nrow(data), 2]

#Extract chromosomes data
idx <- grep("Chromosome annotation", data[1, ])
tmp <- strsplit(data[2:nrow(data), idx], ",")
chromosomes <- NULL
for (i in 1:length(tmp)) {
  chromosomes[i] <- tmp[[i]][1]
}
chromGenes <- as.data.frame(cbind(chromosomes, genes))

#Choose the required portion of the data
idx <- grep("GSM", data[1, ])
data <- data[2:nrow(data), idx]

#Transpose, add column names and labels to the data
data <- as.data.frame(t(data))
data <-
  as.data.frame(apply(data, 2, function(x)
    as.numeric(as.character(x))))
colnames(data) <- cols

#Merge (average) repetitive genes
idx <- NULL
unqGenes <- unique(genes)
newData <- NULL
for (i in 1:length(unqGenes)) {
  sameGenes <- which(genes == unqGenes[i])
  len <- length(sameGenes)

  if (len > 1) {
    newCol <- rowMeans(data[, sameGenes])
    idx <- c(idx, sameGenes[2:len])
  }
  else
    newCol <- data[, sameGenes]
  
  newData <- cbind(newData, newCol)
}

data <- as.data.frame(newData)
colnames(data) <- unqGenes
chromGenes <- chromGenes[-idx, ]

#Remove fully null columns
nullFeatures <- sapply(data, function(x)
  all(is.na(x)))
nullFeatures <- which(nullFeatures == T)
if (length(nullFeatures) > 0) {
  data <- data[, -nullFeatures]
  chromGenes <- chromGenes[-nullFeatures, ]
}

#Remove columns with ####_at_ names
unqGenes <- colnames(data)
rmGenes <- regexpr("_at", unqGenes)
idxrmGenes <- which(rmGenes > 0, arr.ind = T)
if (length(idxrmGenes) > 0) {
  data <- data[, -idxrmGenes]
  chromGenes <- chromGenes[-idxrmGenes, ]
}

#Remove control columns
unqGenes <- colnames(data)
rmGenes <- regexpr("--Control", unqGenes)
idxrmGenes <- which(rmGenes > 0, arr.ind = T)
if (length(idxrmGenes) > 0) {
  data <- data[, -idxrmGenes]
  chromGenes <- chromGenes[-idxrmGenes, ]
}

#Impute the data
data <- as.data.frame(e1071::impute(data))

#Add class column to the data
class <- class[1:nrow(data),]
data <- cbind(data, class[, 2])
colnames(data)[ncol(data)] <- 'class'

#Shuffle the data
data <- data[sample(nrow(data)), ]


#Store data with and without features, features names and corresponding chromosomes
paste0(path, '/', fileName)
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
  chromGenes$genes,
  file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_Justfeatures.csv'),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)
write.table(
  chromGenes$chromosomes,
  file = paste0(path, '/', strsplit(fileName, '\\.')[[1]][1], '_JustChromosomes.csv'),
  sep = ",",
  quote = F,
  row.names = F,
  col.names = F
)

#Print information related to the dataset
cat(
  "# of genes:",
  length(genes),
  "\n# of unique genes:",
  length(unqGenes),
  "\n# of null genes:",
  length(nullFeatures)
)
