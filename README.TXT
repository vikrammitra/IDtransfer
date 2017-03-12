========================
BUILD OUTPUT DESCRIPTION
========================

IN ORDER TO RUN THE ID TRANSFER PIPELINE SET THE FOLLOWING PARAMETERS IN THE CONFIG FILE FOR 3 MAIN STEPS 
1.Set Input Parameters for Lib functions and source file locations.
2.Pre-processing STEPS 
3.ID TRANSFER PROCESS

# STEP1. Lib functions and source file locations ###################################################
# SET THE LOCATION OF 'PreprocessingSourceFile' THAT EXTRACTS INOFRMATION FROM ID-FILES [.pepxmlFiles] AND RAW MS FILES
InputParam[["PreprocessingSourceFile"]] <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/FinalScripts/Pre-process_Source.R' 
# SET THE LOCATION OF 'ID_tranferSourceFile' THAT PERFORMS ID-TRANSFER BETWEEN TWO DATASETS
InputParam[["ID_tranferSourceFile"]] <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/FinalScripts/ID-tranfer_Source.R'
# SET THE LOCATION OF R-FUNCTIONS SCRIPT
InputParam[["Rlibpath"]] <- 'lib/IDtransferFunc.R'
InputParam[["JAR_filePath"]] <- paste(dirname(InputParam$SourceFile),'lib/jar_libs/',sep="/")
# SET THE LOCATION TO SAVE R SESSION WORKSPACE
InputParam[["SaveImageName"]] <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/IDtransfer.RData'

# STEP2. SET PARAMETERS FOR Pre-processing STEPS ###################################################
# SET THE PATH OF DATASET TO PRE-PROCESS; CREATE A DIR FOR EACH DATASET
InputParam[["RootWD"]]   <- 'C:/Users/DellPS/Documents/VM_thesis/Chapter4/VM/PeaksHumanProteome_analysis/Dataset1/'
# SET THE PATH OF SUB-DIR CONTAINING ALL RAW MS FILES 
InputParam[["RawFiles"]] <- 'MSFiles/'
# SET THE PATH OF SUB-DIR CONTAINING .pepxml FILE 
InputParam[["PepXMLFiles"]] <- 'pepxmlFiles/'
# SET THE PATH OF SUB-DIR CONTAINING FEATURE QUANT FILE 
InputParam[["QuantFile"]] <- 'Quant/'
# SET DATASET NAME 
InputParam[["ConsensusFilename"]] <- 'Dataset1'
# SET TO TRUE IN CASE PLOTTING IS REQUIRED
InputParam[["PLOT"]]<- TRUE

# STEP3. SET PARAMETERS FOR ID TRANSFER PROCESS ###################################################
# SET THE REFERENCE FILE NAME
InputParam[["RefSet"]] <- Dataset3
# SET THE SAMPLE FILE NAME
InputParam[["SampSet"]] <- Dataset4
# SET THE MZ TOLERANCE IN Da.
InputParam[["mzTH"]] <- 0.005
# SET THE rt TOLERANCE IN Mins.
InputParam[["rtTH"]] <- 2
# SET THE OUTPUT DATASET NAME OF ID TRANSFFRED FILE
InputParam[["IDtransferFileName"]] <- 'idMat_tabl3_tabl4'