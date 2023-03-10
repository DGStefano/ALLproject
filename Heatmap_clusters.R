library(edgeR)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(rstatix)
library(multcomp)
library(FSA)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gplots)
library(tidyverse)
ES=read.delim("/Volumes/SSD1T/mac/folgiero/data/readcount_ESREC_specific/concatfiles/counonmergedESREC11k/ES_header_count.bed" , header = T)

MUD=read.delim("/Volumes/SSD1T/mac/folgiero/data/readcount_ESREC_specific/concatfiles/counonmergedESREC11k/MUD_header_count.bed" , header = T)

REC=read.delim("/Volumes/SSD1T/mac/folgiero/data/readcount_ESREC_specific/concatfiles/counonmergedESREC11k/REC_header_count.bed" , header = T)

REM=read.delim("/Volumes/SSD1T/mac/folgiero/data/readcount_ESREC_specific/concatfiles/counonmergedESREC11k/REM_header_count.bed" , header = T)

merged_allstatus <- inner_join(inner_join(inner_join(MUD,ES , by=c("chro","start","end")),REM,by=c("chro","start","end")),REC,by=c("chro","start","end"))
merged_allstatus$chro <-as.character(merged_allstatus$chro) 
merged_allstatus$start <-as.character(merged_allstatus$start) 
merged_allstatus$end <-as.character(merged_allstatus$end) 
merged_allstatus$position <- paste(merged_allstatus$chro,merged_allstatus$start,merged_allstatus$end,sep="-")
rownames(merged_allstatus) <- merged_allstatus$position
merged_allstatus <- dplyr::select(merged_allstatus,  select=-c(chro,start,end,position))



group <- rep(c('MUD' , 'ES','REM','ES','REM','REC') , times=c(6,7,1,4,6,8))
# group <- rep(c('MUD' , 'ES','REM','REC') , times=c(6,11,7,8))

# Make DGE list object(EdgeR basic object)
regions <- DGEList(merged_allstatus , group = group)
# Data correction
#common dispersion calculates a common dispersion value for all fragments
regions <- estimateCommonDisp(regions)
#tagwise method calculates gene-specific dispersions
regions <- estimateTagwiseDisp(regions)
#TMM normalization: scaling the fragments count on the bases of the library preparation
regions <- calcNormFactors(regions, method="TMM")
#count per million
cpm_tmm <- cpm(regions, normalized.lib.size=T , log=T)
#normalizazion (z score)
cpm_tmm <- t(scale(t(cpm_tmm)))

#If needed create an heatmap of cmp_tmm matrix

metadata <- data.frame("sample"= colnames(cpm_tmm) ,
                       "group"=rep(c('MUD' , 'ES','REM','ES','REM','REC') , times=c(6,7,1,4,6,8)))
cc <- rep(c("#B8CCE1","#d2eac8", "#f8daac","#d2eac8", "#f8daac","#f1b7b0") , times=c(6,7,1,4,6,8))

# cc <- c("#B8CCE1","#B8CCE1","#B8CCE1","#B8CCE1","#B8CCE1","#B8CCE1","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#d2eac8","#f8daac","#f8daac","#f8daac","#f8daac","#f8daac","#f8daac","#f8daac", "#f1b7b0","#f1b7b0","#f1b7b0","#f1b7b0","#f1b7b0","#f1b7b0","#f1b7b0","#f1b7b0")
#png(filename="/home/stefano/Documents/IFO/folgiero/pictures/heat_11k_ESspecific_add_RECdiff10_logscaled_ward_D2_samplesname_dividedclusters.png", width=7, height=7, units="in", res=600)

bk <- seq(-2,2,0.1)

mycols <- colorRampPalette(colors = c("#3f8ef5","#e8cc81","#ee2326"))(length(bk))
hv <- heatmap.2(as.matrix(cpm_tmm), scale="none" ,
                ColSideColors=cc,
                hclustfun = function(x) hclust(x, method="ward.D2"),
                # dendrogram = c("both"),
                trace = "none",
                #rowsep=c(2463,5071,7643), #cambia i separatori
                sepwidth=c(0,30),
                col = mycols,
                # RowSideColors = rep(c("blue", "pink" , 'red'), c(3784,5128,2143)),
                density.info=c("none"),
                breaks = seq(-2,2.1,0.1),
                margins = c(6,17),
                # Rowv = FALSE, 
                # Colv = FALSE
)
# png("/Users/sdigiove/Documents/PhDwork/Picture/Heatmap_11kafterKMeanscluster2.png" , res = 300 , height = 5 , width = 7 , units = 'in')
#hv
dev.off()


cpm_tmm.df = as.data.frame(cpm_tmm) %>% rownames_to_column(var = "sites")

cpm_tmm1 = cpm_tmm[cpm_tmm.df$sites %in% Cluster1$sites , ]
cpm_tmm0 = cpm_tmm[cpm_tmm.df$sites %in% Cluster2$sites , ]
cpm_tmm2 = cpm_tmm[cpm_tmm.df$sites %in% Cluster3$sites , ]
cpm_tmm3 = cpm_tmm[cpm_tmm.df$sites %in% Cluster4$sites , ]
cpm_tmm4 = cpm_tmm[cpm_tmm.df$sites %in% Cluster5$sites , ]
cpm_tmm5 = cpm_tmm[cpm_tmm.df$sites %in% Cluster6$sites , ]
cpm_tmm6 = cpm_tmm[cpm_tmm.df$sites %in% Cluster7$sites , ]
cpm_tmm7 = cpm_tmm[cpm_tmm.df$sites %in% Cluster8$sites , ]

cpm_tmm_final = rbind(cpm_tmm0 , cpm_tmm1)
cpm_tmm_final = rbind(cpm_tmm_final , cpm_tmm2)

#enhancer

#legend(x="left",y="bottom" ,legend=c("MUD","ES","REM", "REC"), col=c("#B8CCE1","#d2eac8","#f8daac","#f1b7b0"),pch=15)
#dev.off()


#cut branch of heatmap's dendogram
row_clust <- hclust(dist(cpm_tmm), method = 'ward.D2')
plot(row_clust)

#first cut tree in 2 branchs
abline(h = 140, col = "red2", lty = 2, lwd = 2)

clusters_two <- cutree(row_clust,k=2)
C1C2 <- names(clusters_two[clusters_two==1])
C3C4 <- names(clusters_two[clusters_two==2])
#cut tree in four branchs
clusters_four <- cutree(row_clust,k=4,h=140)
C1 <- names(clusters_four[clusters_four==1])
C2 <- names(clusters_four[clusters_four==2])
C3 <- names(clusters_four[clusters_four==3])
C4 <- names(clusters_four[clusters_four==4])

#Define C2 and C4 from C3C4 less C1 and C3, cutree(k=4) divide in a wrong way row_clust
C2 <- setdiff(C1C2,C1)
C4 <- setdiff(C3C4,C3)


#data.frame of different list from branch an save tabs
sites_list <- data.frame("sites"=c(C4))
sites_list <- sites_list %>% tidyr::separate(sites , c("chro","start","end"),sep="-")
# write.table(sites_list,file="/Volumes/SSD1T/mac/folgiero/data/eRNA_cancer/eRNAinclusters/11ksites_RECES_sitesfromcluster_C4.bed" , quote = F , sep="\t" , row.names = F , col.names = F )
sites_list <- data.frame("sites"=c(C1))
sites_list <- sites_list %>% tidyr::separate(sites , c("chro","start","end"),sep="-")

sites_list <- data.frame("sites"=c(C2))
sites_list <- sites_list %>% tidyr::separate(sites , c("chro","start","end"),sep="-")

sites_list <- data.frame("sites"=c(C3))
sites_list <- sites_list %>% tidyr::separate(sites , c("chro","start","end"),sep="-")
#Dunntesst
PT = dunnTest(value ~ status,
              data=final,
              method="bh")
PT$res
print(PT,dunn.test.results=TRUE)