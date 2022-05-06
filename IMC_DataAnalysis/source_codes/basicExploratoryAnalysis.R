#BEGIN COPYRIGHT NOTICE

#basicExploratoryAnalysis -- (c) 2022 Dimitrios Kleftogiannis -- UiB and CCBIO
#Copyright 2022 University of Bergen (UiB), Norway.
#This Program is free software licensed under the MIT License.
#You may only use the source code in this repository in compliance with the license provided in this repository. 
#For more details, please refer to the file named "LICENSE.md".
#This Program is distributed as a service to the research community and is experimental in nature and may have hazardous properties. 
#The Program is distributed WITHOUT ANY WARRANTY, express or implied. In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A PARTICULAR PURPOSE are excluded. 
#Published reports of research using this code (or a modified version) should cite the relevant article of this code
#Comments and bug reports are welcome.
#Email to dimitrios.kleftogiannis@uib.no
#I would also appreciate hearing about how you used this code, improvements that you have made to it.
#You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

#END COPYRIGHT NOTICE

#UTILITY
#This code implements a basic downstream analysis workflow that allows to import segmented IMC data (exported from histoCAT), 
#produce basic statistics, perform basic clustering analysis to identify phenotypically similar cell subsets, and visualise the results using heatmaps and t-SNE

#RUNNING
#We provide a step-by-step execution example, this is useful for illustration purpose and to reproduce parts of the analysis presented in the paper. 
#Remember that this is a minimal example. 

#DEPENDENCIES
#To run the tutorial without problem, is required to install the following packages. 

# load the packages we need, and few more for more advanced analyses
pkgs2load <- c('readxl','plyr', 'dplyr', 'magrittr', 'tidyr', 'flowCore', 'FlowSOM', 
               'data.table', 'Rtsne', 'ggplot2', 'gridExtra', 
               'ggpubr', 'readr', 'RColorBrewer','openxlsx',
               'matrixStats','reshape2','limma','ggrepel',
               'ConsensusClusterPlus','pheatmap','ggridges','grid','corrplot')
sapply(pkgs2load, require, character.only = TRUE)


#To begin with the following datasets must be downloaded from XXXXX
#OSCC.csv, Osteosarcoma.csv, Placenta.csv and Tonsil.csv

mySamples <- c('OSCC','Osteosarcoma','Placenta','Tonsil')

#STEP1: data intergration
#read and concatenate the segmented single cells into one data frame
AllData <- data.frame()
for(myname in mySamples){
  
  myStr<-paste('Processing sample: ',myname,sep='')
  print(myStr)
  #data and code must be in the same folder, otherwise modify paths accordingly
  inputFileName <- paste(myname,'.csv',sep='')
  dataValues <- as.data.frame(read_csv(inputFileName)) 
  #list the antibodies to be used
  markers <- c("Cell_165HoVEGFR2_Ho165","Cell_163DyTenascinC_Dy163","Cell_168ErKi67_Er168",
               "Cell_175LuCD146_Lu175","Cell_144NdYAP1_Nd144",
               'Cell_141PraSMA_Pr141',
               'Cell_142NdEGFR_Nd142',
               'Cell_143NdPodoplanin_Nd143',
               'Cell_145NdCaveolin_Nd145',
               'Cell_146NdCD16_Nd146',
               'Cell_147SmCD163_Sm147',
               'Cell_158GdECad_Gd158',
               'Cell_149SmCD140b_Sm149',
               'Cell_151EuCD31_Eu151',
               'Cell_152SmFAP_Sm152',
               'Cell_155GdFOXP3_Gd155',
               'Cell_156GdCD4_Gd156',
               'Cell_159TbCD68_Tb159',
               'Cell_161DyCD20_Dy161',
               'Cell_162DyCD8_Dy162',
               'Cell_169TmCollagen1_Tm169',
               'Cell_172YbFSP1_Yb172',
               'Cell_176YbCD90_Yb176',
               'Cell_153EuITGA11_Eu153',
               'Area','Eccentricity','Extent','Number_Neighbors','X_position','Y_position')
  data <- dataValues[,markers]
  #provide the actual names,
  colnames(data) <- c('VEGFR2','TenascinC','Ki67','CD146','YAP1',
                      'aSMA','EGFR','Podoplanin','Caveolin_1',
                      'CD16','CD163','ECadherin','CD140b',
                      'CD31','FAP','FoxP3','CD4','CD68',
                      'CD20','CD8','Collagen_1','FSP_1','CD90','Integrin_a11',
                      'Area','Eccentricity','Extent','No_of_Neighbors','X','Y')
  #store all data for further analysis
  a <- data.frame(File=myname,NumCells=nrow(data),data)
  AllData <- rbind(AllData,a)
}


#STEP2: data tranformation and normalisation 
#here we only provide censoring, other transformation is not performed
a<-as.matrix(AllData[,3:(ncol(AllData)-6)])
myQuant <- apply(a, 2, function(x) quantile(x, probs = .99))
#just censored at the 99th percentile, globally
for(i in 1:nrow(AllData)){
  k<-1
  for(j in 3:(ncol(AllData)-6)){
    
    if( AllData[i,j] > myQuant[k]){
      AllData[i,j] <- myQuant[k]
    }
    k<-k+1
  }
}

#generate the number of cells plot
plot_data <- data.frame(table(AllData$File))
colnames(plot_data)[1] <- 'Sample'
colnames(plot_data)[2] <- 'NumCells'

p1 <- ggplot(plot_data) + aes(x=Sample,y=NumCells,fill=Sample)+
  geom_bar(stat='identity',width = 0.28)+
  #facet_wrap(~variable,scales='free',nrow=4)+
  theme_classic()+
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10),
        axis.title.y = element_text( size = 10 ),
        axis.title.x = element_blank(),strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "none")
myfile<-paste('Summary_NumCells.pdf',sep ='')
pdf(myfile)
print(p1)
dev.off()


#STEP3: cluster the data using FlowSOM and ConsensusClusterPlus
#similar clustering results can be obtained using Phenograph or any other algorithm

#initialise a grid of relatively big size compared to the available data dimemsion
n_clusters_x <- 15
n_clusters_y <- 15
n_meta_clusters <- 200
#use as input 24 antibody stains, shape features could be also included 
selected_markers<- c('VEGFR2','TenascinC','Ki67','CD146','YAP1',
                     'aSMA','EGFR','Podoplanin','Caveolin_1',
                     'CD16','CD163','ECadherin','CD140b',
                     'CD31','FAP','FoxP3','CD4','CD68',
                     'CD20','CD8','Collagen_1','FSP_1','CD90','Integrin_a11')

data1 <- as.matrix(AllData[,selected_markers])
#run flowSOM and then consensusClusterPlus
clustering_data <- flowFrame(data1)
start_time = Sys.time()
fsom1 <- FlowSOM(clustering_data,
                 colsToUse = selected_markers,
                 xdim = n_clusters_x, ydim = n_clusters_y,
                 nClus=n_meta_clusters,
                 maxMeta = n_meta_clusters, 
                 rlen = 10,transform = FALSE)
end_time = Sys.time()
end_time-start_time
codes <- fsom1$map$codes
nmc <- n_meta_clusters
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234,plot='pdf')


#res <- calcICL(mc,plot='pdf')
results <- mc
maxK <- nmc
#STEP4: estimate the quality of metaclustering given the PAC metric. Other metrics can be also used.
############## PAC implementation ##############
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

#plot the proportion of ambiguous clustering (PAC) values
df <- data.frame(PAC)
df$K <- c(2:(length(PAC)+1))
o <- ggplot(df, aes(x=K,y=PAC))+
  geom_point(size = 0.68)+
  geom_line()+
  theme_classic() +
  geom_vline(xintercept = 115,linetype="dotted", 
             color = 'red', size=0.88)+
  geom_hline(yintercept = 0.01,linetype="dotted", 
             color = 'red', size=0.88)+
  ylab('Proportion of ambiguous clustering (PAC) ')+
  xlab('Number of metaclusters ')+
  #scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  theme(axis.text.y = element_text( size = 12 ),
        axis.text.x = element_text( size = 12 ),
        axis.title.y = element_text( size = 12 ),
        axis.title.x = element_text( size = 12 ),strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5))+
  guides(colour = guide_legend(override.aes = list(size=4.8)))

myfile<-paste('Clustering_PAC_A.pdf',sep ='')
pdf(myfile,onefile = TRUE)
print(o)
dev.off()

#after visual inspection it is necessary to set PAC cutoffs and decide about merging clusters
#here we use the value of 115, as it is the first point under 1% PAC, other solutions are also acceptable depending on the level of resolution
nmc <- 115
code_clustering1 <- mc[[nmc]]$consensusClass
clusters <- code_clustering1[fsom1$map$mapping[,1]]
cluster_num <- max(clusters)
#produce Z scored values
data <- AllData[,3:(ncol(AllData)-6)] 
data_Z <- scale(as.matrix(data))
data_Z <- data_Z[,selected_markers]

#identify the z-transformed centroids of clusters
centroids<-matrix(0L, nrow = cluster_num, ncol = (ncol(data_Z)))
for (i in 1:cluster_num) {
  
  b<-which(clusters %in% i)
  d<-data_Z[b,1:ncol(data_Z)]
  #df <- data[b,2:(ncol(data)-6)]
  d<-colMeans(d)
  centroids[i,]<-d
}
#estimate the cluster abundance
clustering_table <- as.numeric(table(clusters))
clustering_table<-100*clustering_table/sum(clustering_table)
clustering_table<-round(clustering_table,2)
#fix the centroids phenograph row and column labels
b1 <- as.vector(paste('C',1:cluster_num,' (', clustering_table,'%)',sep=''))
b2 <- as.vector(paste(1:cluster_num,sep=''))
c1 <- colnames(data_Z)[1:(ncol(data_Z))]
rownames(centroids) <- b1
#and here save the results into an excel file to
colnames(centroids) <- c(c1)
centroids <- as.data.frame(centroids)
write.xlsx(centroids, "FlowSom_115_metaclusters.xlsx",asTable = FALSE,overwrite=TRUE)


################################################################
#at this step it is necessary to read and annotate the clusters using cell type defining markers
#an example of the annotated file is provided in file Clusters.xlsx that can be downloaded from our repository
################################################################
#load the annotated clusters -- specify the unassigned
inputFileName <- "Clusters.xlsx"
manualLabels <- read_excel(inputFileName)
manualLabels <- manualLabels[,c('Cluster','Annotations')]
plot_data <- AllData
myClust <- paste('C',clusters,sep='')
plot_data$Clusters <- myClust
dr <- match(plot_data$Clusters, manualLabels$Cluster)
plot_data$CellType <- manualLabels$Annotations[dr]

#estimate summary statistics about the phenotypic abundance of CAFs
samples_count <- plot_data %>% 
  group_by(CellType,File) %>% 
  summarise(samples = n())%>% 
  mutate(perc = samples/sum(samples))

#generate exploratory plots
samples_count$CellType <- factor(samples_count$CellType,levels=c('CAF-1','CAF-2','CAF-3','CAF-4','CAF-5','CAF-6',
                                                                 'Endothelial','Epithelial','Immune','Unassigned'))

o1 <- ggplot(samples_count, aes(x=CellType,y=perc,group=File,fill=File))+
  geom_bar(stat = 'identity',width=0.25)+
  ylab('Relative percentage (%)')+
  theme_bw() +
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_text( size = 10 ),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

o2 <- ggplot(samples_count, aes(x=CellType,y=(samples),group=File,fill=File))+
  geom_bar(stat = 'identity',width=0.25)+
  ylab('Number of cells')+
  theme_bw() +
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_text( size = 10 ),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

myfile<-paste('/Summary_clusters_relative_perc.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o1)
dev.off()

myfile<-paste('Summary_clusters_counts.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o2)
dev.off()

o3 <- ggplot(samples_count, aes(x=CellType,y=samples,group=CellType,fill=CellType))+
  geom_bar(stat = 'identity',width=0.25)+
  scale_fill_manual(limits = levels(samples_count$CellType),values = c('gold','firebrick1','deepskyblue1','darkorchid2','forestgreen',
                                                                       'darkseagreen1','salmon','darkorange','gray88','gray44'))+
  ylab('Number of cells')+
  theme_bw() +
  theme(axis.text.y = element_text( size = 10 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_text( size = 10 ),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))#+coord_flip()

myfile<-paste('Summary_cellType_counts.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o3)
dev.off()

#STEP5: visualise the clusters using heatmap
data <- AllData[,3:(ncol(AllData)-6)] 
data_Z <- scale(as.matrix(data))
data_Z <- data_Z[,selected_markers]
cluster_names <- unique(plot_data$CellType)

#identify the z-transformed centroids of clusters
centroids<-matrix(0L, nrow = length(cluster_names), ncol = (ncol(data_Z)))
myC <- 1
myClusterNames <- c('CAF-1','CAF-2','CAF-3','CAF-4','CAF-5','CAF-6',
                    'Endothelial','Epithelial','Immune','Unassigned')
for (i in myClusterNames) {
  
  b<-which(plot_data$CellType==i)
  d<-data_Z[b,1:ncol(data_Z)]
  #df <- data[b,2:(ncol(data)-6)]
  d<-colMeans(d)
  centroids[myC,]<-d
  myC <- myC+1
}
#estumate the abundance of annotated cell types
clustering_table <- table(plot_data$CellType)
clustering_table<-100*clustering_table/sum(clustering_table)
clustering_table<-round(clustering_table,2)
clustering_table <- as.data.frame(clustering_table)
b1 <- myClusterNames
b2 <- paste('',clustering_table$Freq,'%',sep='')
c1 <- colnames(data_Z)[1:(ncol(data_Z))]
rownames(centroids) <- b2
colnames(centroids) <- c1

#generate info for the row annotation in heatmap
my_row_annot <- data.frame(CellType=b1)
rownames(my_row_annot)<- b2
colnames(my_row_annot)[1] <- 'CellType'
#specify the colors
my_colour = list(
  CellType = c(`CAF-1`='gold',`CAF-2` = 'firebrick1', `CAF-3` = 'deepskyblue1', `CAF-4` = 'darkorchid2',`CAF-5`='forestgreen',`CAF-6`='darkseagreen1',
               Endothelial='salmon',Epithelial='darkorange',Immune='gray88',Unassigned='gray44')
)

breaksList<-seq(min(centroids),max(centroids),by=0.5)
#plot the cluster centroids
myfile<-paste('Heatmap_cellTypes.pdf',sep ='')
pheatmap(centroids[,1:ncol(data_Z)],labels_row=b2,labels_col = c1,
         annotation_row = my_row_annot,annotation_colors = my_colour,
         cluster_rows=TRUE,cluster_cols = TRUE ,fontsize = 6,
         display_numbers = TRUE, number_color = "black",fontsize_number = 4,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         filename=myfile,
         cutree_rows = 4,cutree_cols=5,
         cellheight = 16, cellwidth = 16,angle_col = 45)

#STEP6: perform dimensionality reduction using t-sne

df <- AllData[,selected_markers]
duplicates <- duplicated(df[,1:ncol(df)]) | duplicated(df[,1:ncol(df)], fromLast = TRUE) 
tsne_out_high<-Rtsne(df[!duplicates,1:ncol(df)],max_iter=3000,perplexity=30,seed=1234,num_threads=4)

tsne_data_3 <- as.data.frame(data_Z[!duplicates,])
tsne_data_3$TSNE1 <- tsne_out_high$Y[,1]
tsne_data_3$TSNE2 <- tsne_out_high$Y[,2]
tsne_data_3$Sample <- AllData[!duplicates,'File']
myLabel <- plot_data$CellType
tsne_data_3$CellType <- myLabel[!duplicates]
tsne_data_3$CellType<- factor(tsne_data_3$CellType,levels=c('CAF-1','CAF-2','CAF-3','CAF-4','CAF-5','CAF-6',
                                                            'Endothelial','Epithelial','Immune','Unassigned'))

#visualise the cell types in 2D
o1 <- ggplot(tsne_data_3,  aes(x = TSNE1, y = TSNE2, color = CellType )) +
  geom_point(size = 0.28) +
  theme_bw() +
  scale_color_manual(limits = levels(tsne_data_3$CellType),values = c('gold','firebrick1','deepskyblue1','darkorchid2','forestgreen',
                                                                      'darkseagreen1','salmon','darkorange','gray88','gray44'))+
   guides(fill = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 6))+
  theme(axis.text.y = element_text( size = 8 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

myfile<-paste('Figure_tsne.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o1)
dev.off()

#similar t-sne but split based on the sample information
o2 <- ggplot(tsne_data_3,  aes(x = TSNE1, y = TSNE2, color = CellType )) +
  geom_point(size = 0.28) +
  theme_bw() +
  scale_color_manual(limits = levels(tsne_data_3$CellType),values = c('gold','firebrick1','deepskyblue1','darkorchid2','forestgreen',
                                                                      'darkseagreen1','salmon','darkorange','gray88','gray44'))+
  facet_wrap(~Sample,ncol=2)+
  guides(fill = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 6))+
  theme(axis.text.y = element_text( size = 8 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))

myfile<-paste('Figure_tsne_facets',sep='')
pdf(myfile,onefile = TRUE)
print(o2)
dev.off()

o3 <- ggplot(tsne_data_3,  aes(x = TSNE1, y = TSNE2, color = Sample )) +
  geom_point(size = 0.28) +
  theme_bw() +
  guides(fill = FALSE)+
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))+
  theme(axis.text.y = element_text( size = 8 ),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12,face='bold',lineheight=1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 18, face = "bold",hjust = 0.5))


myfile<-paste('Figure_tsne_samples.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(o3)
dev.off()

#make a combined tsne and color points based on antibody intensity
myC <- 1
n <- list()
for(myMarker in selected_markers){
  
  myD <- tsne_data_3[,c('Sample','TSNE1','TSNE2',myMarker)]
  colnames(myD)[4] <- 'Intensity'
  n[[myC]] <- ggplot(myD,  aes(x = TSNE1, y = TSNE2, color = Intensity)) +
    geom_point(size = 0.108) +
    theme_bw() +
    #facet_wrap(~ Sample, nrow = 2) +
    scale_color_gradientn("Intensity",colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
    #ylim(-98, 98)+
    #xlim(-98, 98)+
    guides(fill = FALSE)+
    ggtitle(myMarker)+
    scale_shape_manual(values = c(1,5),name='Cell type')+
    guides(shape = guide_legend(override.aes = list(size = 2), ncol = 5))+
    theme(axis.text.y = element_text( size = 8 ),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.8),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text = element_text(size = 12,face='bold',lineheight=1),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 18, face = "bold",hjust = 0.5))
  myC <- myC + 1
  
}
#combine individual plots to one
o4 <- ggarrange(n[[1]],n[[2]],n[[3]],n[[4]],n[[5]],n[[6]],n[[7]],n[[8]],n[[9]],n[[10]],
               n[[11]],n[[12]],n[[13]],n[[14]],n[[15]],n[[16]],n[[17]],n[[18]],n[[19]],n[[20]],
               n[[21]],n[[22]],n[[23]],n[[24]],
               ncol = 3,nrow = 3)

myfile<-paste('Figure_tsne_markers',sep='')
pdf(myfile,onefile = TRUE)
print(o4)
dev.off()


