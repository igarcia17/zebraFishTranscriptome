###Author: Sumeet pal Singh

setwd("/media/sumeet/Galaxy/RNA_Seq_Public")

library(stringr)

readFile.cleanColName <- function(file, pattern = "SRR[0-9]*|ERR[0-9]*"){
  data1 <- read.table(file = file,
                    sep = "\t",
                    header = T, row.names = 1, 
                    stringsAsFactors = F)
  
  # colnames(data1)
  
  # str_extract_all(colnames(data1)[2:ncol(data1)], "SRR[0-9]*", simplify = TRUE)
  srr.names <- str_extract_all(colnames(data1)[2:ncol(data1)], 
                             pattern = pattern, 
                             simplify = TRUE)
  
  colnames(data1)[2:ncol(data1)] <- srr.names
  # head(data1)
  return(data1)
}

readFiles.cleanColName <- Vectorize(FUN = readFile.cleanColName, 
                                    vectorize.args = "file")

dir(path = "counts", pattern = "_z11_2.tsv", full.names = T)

all.counts <- readFiles.cleanColName(file = dir(path = "counts", 
                                                pattern = "_z11_2.tsv", 
                                                full.names = T))

str(all.counts)
lengths(all.counts)
# combineCountData

for(ncount in 2:length(all.counts)){
  # print(ncount)
  all.counts[[1]] <- cbind(all.counts[[1]], 
                           all.counts[[ncount]][,2:ncol(all.counts[[ncount]])])
  # print(dim(all.counts[[1]]))
}

final.counts <- all.counts[[1]]
dim(final.counts)

write.csv(x = final.counts, file = "final_counts.csv", quote = F)

sra.info <- read.csv(file = "sra_samples.csv", 
                     header = T, 
                     stringsAsFactors = F, row.names = 1)
head(sra.info)

info.final.counts <- sra.info[colnames(final.counts)[2:ncol(final.counts)],]

colnames(final.counts)[2:ncol(final.counts)] <- apply(info.final.counts, 
                                                      1, 
                                                      paste, collapse = "_")

write.csv(x = final.counts, file = "final_counts_info.csv", quote = F)
