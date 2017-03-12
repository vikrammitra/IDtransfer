###############################################################################
#                                                                            
#Script used to run pre-processing Identification transfer                  
# Author: Vikram Mitra                                       
# v.mitra@rug.nl                                                                             
###############################################################################

###############################################################################
rm(list = ls(all = TRUE)) # Clean working space
###############################################################################

InputParam           <- vector("list")

## Set Input Parameters for Lib functions and source file##

InputParam[["PreprocessingSourceFile"]] <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/FinalScripts/Pre-process_Source.R'
InputParam[["ID_tranferSourceFile"]] <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/FinalScripts/ID-tranfer_Source.R'
InputParam[["Rlibpath"]] <- 'lib/IDtransferFunc.R'
InputParam[["JAR_filePath"]] <- paste(dirname(InputParam$SourceFile),'lib/jar_libs/',sep="/")
InputParam[["SaveImageName"]] <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/IDtransfer.Rdata'

# Start Pre-processing for Dataset 1 ###################################################

InputParam[["RootWD"]]   <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/Dataset1/'
InputParam[["RawFiles"]] <- 'MSFiles/'
InputParam[["PepXMLFiles"]] <- 'pepxmlFiles/'
InputParam[["QuantFile"]] <- 'Quant/'
InputParam[["ConsensusFilename"]] <- 'Dataset1'
InputParam[["PLOT"]]<- TRUE

setwd(InputParam$RootWD)
source(InputParam$PreprocessingSourceFile)


# Start Pre-processing for Dataset 2 ###################################################

InputParam[["RootWD"]]   <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/Dataset2/'
InputParam[["RawFiles"]] <- 'MSFiles/'
InputParam[["PepXMLFiles"]] <- 'pepxmlFiles/'
InputParam[["QuantFile"]] <- 'Quant/'
InputParam[["ConsensusFilename"]] <- 'Dataset2'
InputParam[["PLOT"]]<- TRUE

setwd(InputParam$RootWD)
source(InputParam$PreprocessingSourceFile)

# Start Pre-processing for Dataset 3 ###################################################

InputParam[["RootWD"]]   <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/Dataset3/'
InputParam[["RawFiles"]] <- 'MSFiles/'
InputParam[["PepXMLFiles"]] <- 'pepxmlFiles/'
InputParam[["QuantFile"]] <- 'Quant/'
InputParam[["ConsensusFilename"]] <- 'Dataset3'
InputParam[["PLOT"]]<- TRUE

setwd(InputParam$RootWD)
source(InputParam$PreprocessingSourceFile)

# Start Pre-processing for Dataset 4 ###################################################

InputParam[["RootWD"]]   <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/Dataset4/'
InputParam[["RawFiles"]] <- 'MSFiles/'
InputParam[["PepXMLFiles"]] <- 'pepxmlFiles/'
InputParam[["QuantFile"]] <- 'Quant/'
InputParam[["ConsensusFilename"]] <- 'Dataset4'
InputParam[["PLOT"]]<- TRUE

setwd(InputParam$RootWD)
source(InputParam$PreprocessingSourceFile)

# Perform ID transfer ##################################################################

###############Comparison1
InputParam[["RefSet"]] <- Dataset1
InputParam[["SampSet"]] <- Dataset2
InputParam[["mzTH"]] <- 0.005
InputParam[["rtTH"]] <- 2
InputParam[["IDtransferFileName"]] <- 'idMat_tabl2_tabl1'

print("Performing Id Tranfer between Ref->Dataset2 and Samp->Dataset1")
source(InputParam$ID_tranferSourceFile)

###############Comparison2
InputParam[["RefSet"]] <- Dataset1
InputParam[["SampSet"]] <- Dataset3
InputParam[["mzTH"]] <- 0.005
InputParam[["rtTH"]] <- 2
InputParam[["IDtransferFileName"]] <- 'idMat_tabl3_tabl1'

print("Performing Id Tranfer between Ref->Dataset3 and Samp->Dataset1")
source(InputParam$ID_tranferSourceFile)

###############Comparison3
InputParam[["RefSet"]] <- Dataset2
InputParam[["SampSet"]] <- Dataset3
InputParam[["mzTH"]] <- 0.005
InputParam[["rtTH"]] <- 2
InputParam[["IDtransferFileName"]] <- 'idMat_tabl3_tabl2'

print("Performing Id Tranfer between Ref->Dataset3 and Samp->Dataset2")
source(InputParam$ID_tranferSourceFile)

###############Comparison3
InputParam[["RefSet"]] <- Dataset3
InputParam[["SampSet"]] <- Dataset4
InputParam[["mzTH"]] <- 0.005
InputParam[["rtTH"]] <- 2
InputParam[["IDtransferFileName"]] <- 'idMat_tabl3_tabl4'

print("Performing Id Tranfer between Ref->Dataset3 and Samp->Dataset4")
source(InputParam$ID_tranferSourceFile)