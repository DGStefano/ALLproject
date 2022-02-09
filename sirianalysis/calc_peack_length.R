#leggo i parametri passati per bash → nome del file
args <- commandArgs(trailingOnly = TRUE)
#leggo la tabella con nome del file passato → file ottenuto da multicovs tra sprted_dedu.bam e blacklist_filtered
peak_tab <- read.table(args[1])

file_name=tail(unlist(strsplit(args[1],"/")),1)

#prendo inizio e fine dei picchi e faccio sottrazione per trovare la lunghezza
starts <- peak_tab[,2]
finish <- peak_tab[,3]
length_peak <- finish-starts

peak_tab <-cbind(peak_tab,length_peak)
#salvataggio nella stessa directory iniziale, il nuovo file sarà length_nomefileiniziale
write.table(peak_tab, file=paste(args[2],paste("length_",file_name,sep=""), sep="/"), col.names = F , row.names=F , quote = F , sep = "\t")