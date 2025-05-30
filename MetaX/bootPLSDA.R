##' @param x An object where samples are in rows and features are in columns. 
##' This could be a simple matrix, data frame. 
##' @param y A numeric or factor vector containing the outcome for each sample.
##' @param ncomp The maximal number of component for PLS-DA
##' @param sample A vector contains the sample used for the model
bootPLSDA=function(x,y,ncomp=2,sample=NULL,test=NULL,split=0,
                   method = "repeatedcv",repeats=250,number=7,...)


library(metaX)
# data input 
data <- read.csv("./01_Data/01.MetQuant/meta_intensity_combined.csv",row.names = 1)
data <- data[,c(14:40)]
data <- log2(data)
data <- t(data)

# group input
group <- read.xlsx("./01_Data/01.MetQuant/sam_infor_combined.xlsx")
group <- group[c(1:27),]

# bootPLSDA
target_data <- data[grep("MOLM13_WT|MOLM13_6W",rownames(data)),]
target_group <- group[grep("MOLM13_WT|MOLM13_6W",group$id),3]

my_plsda <- metaX::bootPLSDA(target_data,
                             target_group,
                             method = "LOOCV")
