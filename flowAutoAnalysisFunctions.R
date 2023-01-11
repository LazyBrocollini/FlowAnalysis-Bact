lowMarginCell <- function(flowFrameInput, 
                          lowCutOff=1000,
                          folderName1="FlowImagesLowMargin",
                          resolution1=100){
  # lowMarginCell takes input from the user on the acceptable the cut off for FSC-A below which
  # values should be discarded. The function outputs a flowSet with filtered data.
  # 
  # Default values:
  # lowerCutOff = 1000

  stopifnot("The input for lowCutOff must be a positive integer." = is.numeric(lowCutOff))
  
  #specify the LOW cut offs and the remove from data
  lowMargin <- which(flowFrameInput[,"FSC-A"] < lowCutOff)
  
  saveImage(flowFrameInput, channel1="FSC-H", channel2="SSC-H",folderName=folderName1, resolution=resolution1)
  return (flowFrame_Marg <- flowFrameInput[-lowMargin])
  
  #to QC, un-comment the code below:   
  # # #visualize to QC (redo with ggcyto, when possible)
  # lowCells <- flowFrame[lowMargin]
  # # print(lowCells)
  # # print(typeof(lowCells))
  # # 
  # # print("plotting QC graphs")
  # A <- exprs(lowCells)[, c("FSC-A", "FSC-H")]
  # B <-exprs(flowFrame_Marg)[, c("FSC-A", "FSC-H")]
  # plot(A, pch=".", ylim = c(0, 30000), xlim = c(0, 10000))
  # points(B, pch=".", col="blue")
  # m <- 1:30000
  # y<- m +1
  # x<- rep(1000, each = 30000)
  # lines(x, y, type = "l", lty = 1, col="red")
}


transLogicle <- function(flowFrameInput,
                         folderName1="FlowImagesTrans",
                         resolution1 =256){
  
  #the function transforms the height and area events onto a logicle scale.
  #first the function makes a data based decision on transformation condtions
  #then the data is transformed
  
  trans <- estimateLogicle(flowFrameInput, c("FSC-H", "SSC-H", "FSC-A", "SSC-A"))
  flowFrameInput_trans <- transform(flowFrameInput, trans)
  
  #saving images
  saveImage(flowFrameInput_trans, channel1="FSC-H", channel2="SSC-H",folderName=folderName1, resolution=resolution1)
  return(flowFrameInput_trans)
}


SingletStats <- function(flowFrameInput, 
                         folderName1= "FlowImagesSinglets"){
  #The function takes in the flowFrame object, transforms it to flowSet and then GatingSet
  #We need GatingSet to use it for adding adding the population gate
  #We need flowFrame/flowSet objects to plot the results for QC
  #Returns statistic for the singlet gate and saves an image to the folder of specified name
  
  fSet <- as(flowFrameInput, "flowSet")
  gs <- GatingSet(fSet)
  sg <- gate_singlet(flowFrameInput, height = "FSC-H", area = "FSC-A", prediction_level = 0.99, filterId = "Singlets")
  
  gs_pop_add(gs, sg, name="Singlets")
  recompute(gs)
  
  saveImage(flowFrameInput, gateName = sg, channel1="FSC-H", channel2="FSC-A", folderName = folderName1)
  return(gs_pop_get_stats(gs))
 }



#function remove doublets was taken from http://github.com/saeyslab/PeacoQC, where code is distributed by GPL(<=3) license
RemoveDoublets <- function(ff,
                           channel1="FSC-A",
                           channel2="FSC-H",
                           #in the original code from saeyslab: nmad = 4
                           nmad=1,
                           verbose=FALSE,
                           output="frame"){
  
  if (!is(ff, "flowFrame"))
    stop("ff should be a flowframe.")
  
  # Calculate the ratios
  ratio <- flowCore::exprs(ff)[,channel1] /
    (1+ flowCore::exprs(ff)[,channel2])
  
  # Define the region that is accepted
  r <- stats::median(ratio)
  r_m <- stats::mad(ratio)
  if(verbose) message(paste0("Median ratio: ", r,", width: ", nmad*r_m))
  
  # Make selection
  selection <- ratio < r+nmad*r_m
  
  new_ff <- ff[selection, ]
  
  ##removed the section below, unnecessary for bacteria counting purpose
  # if (!("Original_ID" %in% colnames(flowCore::exprs(new_ff)))){
  #   new_ff <- AppendCellID(new_ff, which(selection))}
  
  if (output == "full"){
    return(
      list("flowframe"=new_ff,
           "indices_doublets"=which(selection == FALSE)))
  } else if (output == "frame"){
    return(new_ff)
  }
}


gatingFlowSet <- function(flowFrameInput,
                          folderName1="FlowImagesGate"){
  #uses the autogating function in flowStats which is implemented as lymphGate
  #then the image of the graph is saved to a folder
  #the count for the gate is returned
  
  flowSetInput <- as(flowFrameInput, "flowSet")
  resultsGate <- flowStats::lymphGate(channels=c("FSC-H", "SSC-H"), flowSetInput)
  
  ##unquote to QC by plotting
  #x<- ggcyto(flowSetInput[[1]], aes("FSC-H", "SSC-H"))+geom_hex(bins=500)+geom_gate(resultsGate)
  # print(x)
  
  gs <- GatingSet(flowSetInput)
  gs_pop_add(gs, resultsGate, name ="autoGate")
  recompute(gs)
  saveImage(flowFrameInput, gateName=resultsGate,channel1="FSC-H", channel2="SSC-H",folderName=folderName1)
  return(gs_pop_get_stats(gs, "autoGate"))
}


saveImage<-function(flowFrameInput, 
                    gateName, 
                    channel1="FSC-H", 
                    channel2="SSC-H", 
                    folderName="FlowImages", 
                    resolution=500){
  # saveImage function allows to save images of graphs made with flowFrame and gates
  # flowFrameInput is data in flowFrame format, 
  # gateName is the variable containing gate information, provide variable with Gate info not a sting of characters
  # channel1 and channel2 are the axis of the graph
  # folderName is a folderName where the images will be saved (if the folder exists, the images will be deposited there)
  # if the folder does not exist, one will be created
  # resolution = number of bins for the image, the bigger the number, the finer the image (more events displayed)
  
  mainDir <- getwd()
  
  if(file.exists(folderName)){
    print(paste(folderName, "is a directory for the images. It already existed."))
  }else{
    dir.create(file.path(mainDir, folderName))
    print(paste(folderName, "is a directory for the images. It was created. "))
  }
  
  ImageOutPut <- paste(mainDir, folderName, sep="/")
  
  #Saving graphs in png format
  flowSetInput <- as(flowFrameInput, "flowSet")
  ch1 <- enquo(channel1)
  ch2 <- enquo(channel2)
  
  if(missing(gateName)){
    x <- ggcyto(flowSetInput[[1]], aes(x=!!ch1, y=!!ch2))+geom_hex(bins=resolution) 
  }else{
    x <- ggcyto(flowSetInput[[1]], aes(x=!!ch1, y=!!ch2))+geom_hex(bins=resolution)+ geom_gate(gateName)
  }
  
  flowFrameName <- paste(keyword(flowFrameInput)$GUID, ".png", sep="")
  
  ggsave(filename=flowFrameName, plot=x, path=ImageOutPut)
}