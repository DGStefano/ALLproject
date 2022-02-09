library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

#leggo la tabella con nome del file passato → file con nscore
peak_tab <- read.table(args[1])

file_name=tail(unlist(strsplit(args[1],"/")),1)

peak_tab <- mutate(peak_tab,ranking_index=ntile(desc(peak_tab$V6),100))
write.table(peak_tab, file=paste(args[2],paste("RI_",file_name,sep=""), sep="/"), col.names = F , row.names=F , quote = F , sep = "\t")
