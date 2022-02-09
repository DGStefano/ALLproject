#leggo i parametri passati per bash → nome del file
#args[1] → nome del file
#args[2] → directory di salvataggio
#args[3] → file bam per samtools view -c

args <- commandArgs(trailingOnly = TRUE)

Normalized_peack <- function(reads,length,total_reads) {
  norm_peak <- ((reads / length)*10^6)*10^3/total_reads
  #non salva in notazione scientifica
  format(norm_peak, scientific = F)
  return(norm_peak)
}




#leggo la tabella con nome del file passato → file ottenuto da multicovs tra sprted_dedu.bam e blacklist_filtered
peak_tab <- read.table(args[1])

#array con valori normalizzati
norm_peak_array <- c()

#numero di reads per campione
total_reads <- as.integer( system(paste("samtools view -c ", args[3] , sep=""), intern = T))

for (i in 1:nrow(peak_tab)){
 norm_peak_array[i] <- Normalized_peack(peak_tab[i,10], peak_tab[i,11], total_reads)
}

file_name=tail(unlist(strsplit(args[1],"/")),1)

peak_tab <-cbind(peak_tab,norm_peak_array)
#tabella ordinata dal maggiore al minore per il valore normalizzato
peak_tab <- peak_tab[order(peak_tab$norm_peak_array, decreasing = T),]
write.table(peak_tab, file=paste(args[2],paste("nscore_",file_name,sep=""), sep=""), col.names = F , row.names=F , quote = F , sep = "\t")
