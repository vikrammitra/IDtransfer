###############################################################################
#                                                                            
#Source-Script used to run pre-processing Identification transfer                  
# Author: Vikram Mitra                                       
# v.mitra@rug.nl                                                                             
###############################################################################

#Load R functions
cat("\n Starting pre-processing steps for ID transfer workflow..")
source(paste(dirname(InputParam$SourceFile),InputParam$Rlibpath,sep = "/"))

#1.Run Xrea Calculation

cat("\n Running Xrea Calculations..")

#RunXreacalc(InputParam$RawFiles,InputParam$JAR_filePath)

#2.Run pepXML parser
cat("\n Parsing pepXML files..")

pepXML_csvfile<-ParsepepXMLfiles(InputParam$PepXMLFiles,InputParam$JAR_filePath,InputParam$RawFiles)


#3.Create peaklist files from Peaks Feature lists and create consensus feature lists
quantList<- CreatePeakList(InputParam$QuantFile)
noOffiles<- length(quantList)

df<- data.frame()


for(i in 1:noOffiles){
  
  ConsensusMatrix_temp <- merged_pepXML_FeatureMatrix(quantList[i],pepXML_csvfile[i])
  df <- rbind(df,ConsensusMatrix_temp)
  
}

assign(InputParam[["ConsensusFilename"]],df)
df<-NULL

cat(paste("\n Consensus files created for -",InputParam$ConsensusFilename))

########### Plotting Results ###########

if( InputParam$PLOT == TRUE ){
  
  #4.Plotting Xrea V/S intensity
  cat(paste("\n Plotting Xrea vs Intensity for -",InputParam$ConsensusFilename))
  FigName = paste(InputParam[["ConsensusFilename"]],"_Xrea_vs_Intensity.png")
  g<-XreaScatterPlot(pepXML_csvfile,InputParam$RawFiles,InputParam[["ConsensusFilename"]])
  ggsave(FigName,g,height = 5,width = 10,dpi = 720)
  
  
  #5.Plotting Intensity distribution
  cat(paste("\n Plotting Intensity distributions for -",InputParam$ConsensusFilename))
  ecdf_plot(get(InputParam$ConsensusFilename),InputParam$ConsensusFilename,InputParam$RawFiles)

  
}

###########END##############