#BEGIN COPYRIGHT NOTICE

#plotTissueStructures -- (c) 2022 Dimitrios Kleftogiannis -- UiB and CCBIO
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
#The code implements a visualisation function that displays IMC data and annotated cells in a spatially-resolved fashion.
#The code has been originally adapted from cytomapper and EBImage packages.

#RUNNING
#We provide a step-by-step execution example, this is useful for illustration purpose and to reproduce parts of the analysis presented in the paper. 
#Remember that this is a minimal example. The example has to be adapted to reproduce the results presented in the paper for four different tissues.

#DEPENDENCIES
#To run the tutorial without problem, is required to install the following packages:

library(cytomapper)

#To reproduce the example please download the tiff images and masks from our data repository at: XXXXX
#To annotate individual cells from the clustering procedure, code basicExploratoryAnalysis.R must be executed first.
#Data frame named plot_data is required which contains the cell ID and the CellType field which is given from clustering and manual inspection.
#If you face problems to obtain that, please contact Dimitrios.


#############################################################################3
#Example that reproduces OSCC tissue visualisation
#############################################################################
#read the mask
myMask <- loadImages('Low_s0_a2_ac_ilastik_s2_Probabilities_mask.tiff')
myMask <- scaleImages(myMask, 2^16-1)
# Read in images
#here we just read one antibody, repeat the same part if more channels are available
myStack1 <- loadImages('OSCC1 high a11_White147Sm_CD163_imc.tiff')
myStack3 <- loadImages('OSCC1 high a11_White162Dy_CD8_imc.tiff')
#five the appropriate tissue names in mask and channels
mcols(myMask)$ImageNb <- c("OSCC")
mcols(myStack1)$ImageNb <- c("OSCC")
mcols(myStack1)$ImageName <- c("OSCC")
mcols(myStack3)$ImageNb <- c("OSCC")
mcols(myStack3)$ImageName <- c("OSCC")
#name the channels
channelNames(myStack1) <- c("CD163",'a','b')
channelNames(myStack3) <- c("CD8",'a','b')
cur_channel1 <- getChannels(myStack1,'CD163')
cur_channel3 <- getChannels(myStack3,'CD8')
#merge the available stacks in pairs of two
all_stacks <- mergeChannels(cur_channel1,cur_channel3)

#make a basic visualisation of the tissue with two channels
plotPixels(image =all_stacks,colour_by = c('CD163','CD8'),img_id = "ImageNb")

#generate a single cell object
sce <- measureObjects(myMask, all_stacks,img_id = "ImageNb")
#take a copy of the sce object and add manually annotations
A <- sce 
A$CellType <- 'OtherCells'

#remember to load the plot_data from basicExploratoryAnalysis script
load('plot_data.Rdata')
a <- which(plot_data$File=='OSCC')
tmp <- plot_data[a,]
A$CellType <- tmp$CellType
#the number of cells must match the mask and the loaded channels
A$CellNb <- 1:nrow(tmp)
#select the cell types you want to visualise, or you can skip this part if you want to visualise everything
#here we show only CAFs
A$CellType[A$CellType!='CAF-1' &
             A$CellType!='CAF-2' &
             A$CellType!='CAF-3' &
             A$CellType!='CAF-4' &
             A$CellType!='CAF-5' &
             A$CellType!='CAF-6'] <- 'OtherCells'

g1 <- plotCells(myMask, object = A,
                img_id = "ImageNb", cell_id = "CellNb",
                colour_by = c('CD163','CD8'),
                outline_by = "CellType",
                colour = list(CD163 = c("black", "red"),
                              CD8=c('black','blue'),
                              CellType = c(OtherCells = "gray88",
                                           `CAF-1` = "gold",
                                           `CAF-2` = "firebrick1",
                                           `CAF-3`="deepskyblue1",
                                           `CAF-4`="darkorchid2",
                                           `CAF-5`="forestgreen",
                                           `CAF-6`="darkseagreen1")),
                scale_bar = list(length = 10,
                                 cex = 1,
                                 lwidth = 10,
                                 colour = "white",
                                 position = "bottomleft",
                                 margin = c(5,5),
                                 frame = 3),
                image_title = list(text = c('OSCC'),
                                   position = "topleft",
                                   colour = "white",
                                   margin = c(2,10),
                                   font = 1,
                                   cex = 1),
                legend = list(colour_by.title.font = 2,
                              colour_by.title.cex = 1,
                              colour_by.labels.cex = 0.5
                ),
                return_plot = TRUE,thick=TRUE)

myfile<-paste('OSCC_tissuePlot_v1.pdf',sep='')
pdf(myfile,onefile = TRUE)
print(g1)
dev.off()

