library(flowCore)
library(flowAI)
library(flowViz)
library(ggcyto)
library(flowStats)
source("flowAutoAnalysisFunctions.R")


fcs.dir <- "C:/Users/veron/Desktop/FlowR/Data/"
#first use a list of files to readFCS and put in a list fList
fList <- lapply(dir(fcs.dir, full.names = TRUE), read.FCS)

#using keywords assign a meaningful/readable name
names(fList) <- sapply(fList, keyword, "GUID")

#change fList to flowSet object type
fSet <- as(fList, "flowSet")

#clean with auto quality control through flowAI
fSet_clean <- fsApply(fSet, flow_auto_qc)

#clean the data with lowMarginCell (default cut off at 1000 for FSC-A)
fSet_clean_marg <- fsApply(fSet_clean, lowMarginCell)

#Transformation
fSet_clean_marg_trans <- fsApply(fSet_clean_marg, transLogicle)

fSet_clean_trans <-fsApply(fSet_clean, FUN = transLogicle, folderName1 = "FlowImagesNonTruncated")

#Singlet Statistics, outputs a list with statistics 
#Because bacteria are very close to the baseline of flow cytometer, 
#I am using transformed data for doublet count and exclusion
SingletStatsSummary <-fsApply(fSet_clean_marg_trans, SingletStats) 
print(SingletStatsSummary)

#exclude the doublets
fSet_clean_marg_trans_singl <- fsApply(fSet_clean_marg_trans, RemoveDoublets)

#outputs a gated image and a count of cells
fSet_clean_marg_trans_singl_gate <- fsApply(fSet_clean_marg_trans_singl, gatingFlowSet)
print(fSet_clean_marg_trans_singl_gate)


#gate and count viable population
#output graph images into a folder

