#!/usr/bin/Rscript

# diptest, foreach, and doParallel must be installed
#library(diptest)
library(foreach)

# 3 input arguments - input file, output file, number of threads to use
args = commandArgs(trailingOnly = TRUE)
file_in=args[1]
file_out=args[2]
threads=as.numeric(args[3])

# Create Parallel Backend to multithread diptest
cl <- parallel::makeCluster(threads)
doParallel::registerDoParallel(cl)

data <- read.delim(file_in, header=FALSE)

matrix <- foreach(rep=1:10, .packages="diptest", .combine='cbind', .inorder=FALSE) %dopar% { # run 10 instances in parallel
  toReturn <- c()
  for(i in 1:nrow(data)){
    x <- as.numeric(unlist(strsplit(data[i,4],","))) # get comma separated values in c4 of input (pileup of methylation)
    y <- jitter(x,amount=0.1)
    DIP <- dip.test(y)
    toReturn[i] <- DIP$p.value
  }
  toReturn
}
data$V4 <- -log(rowMeans(matrix)+.00001)
write.table(data, file_out, row.names=FALSE, col.names=FALSE, sep="\t")