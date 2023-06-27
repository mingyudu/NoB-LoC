############ Draw the heatmap for each bicluster of the proteomics data

### matrix: microbiome data: #sample by #protein
### infer.result: inference result including biclustering of each protein
### prot_df: annotation information of the proteins
### sample_df: annotation information of the samples
### scale: character indicating if the values in the heatmap should be scaled('column') or not('none').
### outfile: output file of the heatmaps

plot_bicluster_heatmap <- function(matrix, infer.result, prot_df, sample_df, 
                                   scale='column', outfile){
  
  best.S = infer.result$best.S
  best.w = infer.result$best.w
  best.c = infer.result$best.c
  best.K = infer.result$best.K
  best.m = infer.result$best.m
  best.n = infer.result$best.n
  best.p = infer.result$best.p
  n.subject <- nrow(matrix); n.protein <- ncol(matrix)
  
  library(pheatmap)
  print('plotting heatmaps for biclusters ...')
  pdf(outfile)
  for (s in 0:best.S){ 
    print(paste0('bicluster ',s))
    y.gs <- as.matrix(matrix[,best.w==s])
    if (s>0){
      tmp <- cbind(best.c[s,], y.gs)
      tmp <- tmp[order(tmp[,1]),]
      y.gs <- tmp[,-1]
      
      
      y.tick <- rep(NA, best.K[s])
      y.tick[1] <- n.subject - best.m[s] # number of inactive samples in s-th protein cluster
      for(i in 2:best.K[s]) y.tick[i] <- y.tick[i-1] + best.n[s, (i-1)]
    }
    
    colnames(y.gs) <- (1:n.protein)[best.w==s]
    
    annot_cols<-data.frame("Protein"=prot_df$Type[as.integer(colnames(y.gs))])
    annot_rows<-data.frame("Temperature"=sample_df$Temperature[as.integer(rownames(y.gs))],
                           "Nutrient"=sample_df$Nutrient[as.integer(rownames(y.gs))])
    rownames(annot_cols)<-colnames(y.gs)
    rownames(annot_rows)<-rownames(y.gs)
    color_list2<-list(Protein=c(Ribo='blue',Chap='red',PS='orange',plim='cyan',nlim='yellow',Phycobili='green',Unknown='white'),
                      Temperature=c('20'='slateblue4','24'='lawngreen','28'='mediumvioletred'),
                      Nutrient=c('N:P=1.7'='darkred','N:P=80'='deepskyblue'))
    
    if (s>0){
      pheatmap(y.gs,scale=scale, fontsize = 6,
               cluster_cols = F,cluster_rows = F,
               annotation_colors = color_list2,
               annotation_row = annot_rows,annotation_col = annot_cols,
               fontsize_row = 6,fontsize_col = 6,
               gaps_row = y.tick,main = paste0('bicluster ',s))
    }else{
      pheatmap(y.gs,scale=scale, fontsize = 6,
               cluster_cols = F,cluster_rows = F,
               annotation_colors = color_list2,
               annotation_row = annot_rows,annotation_col = annot_cols,
               fontsize_row = 6,fontsize_col = 6,
               main = paste0('bicluster ',s))
    }
  }
  dev.off()
}

library(readxl)

##### read in data

# info of protein id and name
protein_id_list = readRDS('./data/protein_id_name.rds')
protein_id = protein_id_list$ID

# proteomics data matrix
data <- read_excel('./data/prot2peaks_data.xlsx')
data<-subset(data,select=-c(1:4))
data<-t(as.matrix(data))
colnames(data)<-protein_id

n.subject <- nrow(data); n.protein <- ncol(data)
rownames(data)<-1:n.subject
colnames(data)<-1:n.protein

##### read in protein group names

# get the protein type data

ribo<-read_excel('./data/Protien_gps_names_PHI.xlsx',sheet=2)
ribo<-ribo$...5[3:57]
chap<-read_excel('./data/Protien_gps_names_PHI.xlsx',sheet=3)
chap<-chap$...5[3:18]
PS<-read_excel('./data/Protien_gps_names_PHI.xlsx',sheet=4)
PS<-PS$...5[3:18]
# protein Q7U8D8 are both chap and PS!
phyco<-read_excel('./data/Protien_gps_names_PHI.xlsx',sheet=5)
phyco<-phyco$...5[3:24]
plim<-read_excel('./data/Protien_gps_names_PHI.xlsx',sheet=6)
plim<-plim$...5[3:24]
nlim<-read_excel('./data/Protien_gps_names_PHI.xlsx',sheet=7)
nlim<-nlim$...5[3:24]

##### annotate proteins with its protein type

protein_id_list = readRDS('./data/protein_id_name.rds')
protein_id = protein_id_list$ID
prot_df<-data.frame(ID=protein_id,Type=rep('Unknown',length(protein_id)),stringsAsFactors=FALSE)
prot_df$Type[prot_df$ID %in% ribo]='Ribo'
prot_df$Type[prot_df$ID %in% chap]='Chap'
prot_df$Type[prot_df$ID %in% PS]='PS'
prot_df$Type[prot_df$ID %in% phyco]='Phycobili'
prot_df$Type[prot_df$ID %in% plim]='plim'
prot_df$Type[prot_df$ID %in% nlim]='nlim'
table(prot_df$Type)

##### annotate samples with temperature and nutrient conditions

sample_df<-data.frame(Num=1:30,
                      Temperature=rep(c(rep('20',2),rep('24',2),rep('28',2)), 
                                      5),
                      Nutrient=rep(c('N:P=1.7','N:P=80'),15))

##### plot the heatmap for biclusters

infer.result = readRDS('./results/infer.result.rds')
outfile <- "./results/biclusters_heatmap.pdf"
plot_bicluster_heatmap(data, infer.result, prot_df, sample_df, scale='column', outfile)

