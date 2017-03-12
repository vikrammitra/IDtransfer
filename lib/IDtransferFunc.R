## -> RunXreacalc - Runs Java Jar lib to calculate Xrea Values for each MS2 Spectrums:
## -> ParsepepXMLfiles - Parses 'XML' file using Java jar lib
## -> CreatePeakList - Parses 'feature.csv' file created using PEAKS,
#                      extracts quantitative values for each feature and creates peaklists and corrects RT and MZ
## -> merged_pepXML_FeatureMatrix - merges peaklists with ID information file
## -> annotateMat - tranfers peptide information between two peaklists using user defined mz and RT thresholds
## -> XreaScatterPlot - Plots Scatter plots and distribution of Feature Intensity in a data set 
## -> ecdf_plot - Plots density distribution  plots of Feature Intensity in a data set of MS1, MS2 and MS2 (identified PSMs)


library(dplyr)
library(splitstackshape)

RunXreacalc<- function(mzXMLpath,JAR_filePath){
  
  mzXML_files<-list.files(mzXMLpath, pattern = ".mzXML", full.names = TRUE, ignore.case = TRUE)
  xrea_runnable_jar<-list.files(JAR_filePath, pattern = "xrea_", full.names = TRUE, ignore.case = TRUE)
  
  cat(paste("\n Total Number of raw Files to process -",length(mzXML_files)))
  currentwd<-getwd()
  ##Xrea computation###
  for(i in 1:length(mzXML_files)){
    
    xrea_jar_path <- paste(xrea_runnable_jar,sep="")
    mzxml_path<- paste(currentwd,"/",mzXML_files[i],sep="")
    xrea_runcommand<-paste("java -jar ",xrea_jar_path," ",mzxml_path,sep="")
    system(xrea_runcommand)
    
    cat(paste(paste(paste("\n\tXrea calculated for",i),"out of"),length(mzXML_files)))
  }
  
}

ParsepepXMLfiles <- function(pepXMLpath,JAR_filePath,rawMZxmlpath){
  
  currentwd<-getwd()
  ##Parse pepxml files##
  pepxml_jar<-list.files(JAR_filePath, pattern = "pepXML", full.names = TRUE, ignore.case = TRUE)
  pepxml_inputfile<-paste(currentwd,"/",list.files(pepXMLpath,pattern =".pep.xml",full.names = TRUE, ignore.case = TRUE),sep="")
  pepxmlparsecommand<-paste("java -jar ",pepxml_jar," ",pepxml_inputfile,sep="")
  system(pepxmlparsecommand)
  pepxml_oututfile<-gsub(".pep.xml", "_pepxml.csv",pepxml_inputfile)
  pepxml_oututfile<-read.csv(pepxml_oututfile,header=T,stringsAsFactors = F)
  pepxml_oututfile <- data.frame(sapply(pepxml_oututfile, as.character))
  listofMSruns<- paste(getwd(),"/",list.files(rawMZxmlpath,pattern ="Raw.csv",full.names = TRUE, ignore.case = TRUE),sep="")
  pepxml_oututfile[pepxml_oututfile == ""]<-NA
  
  mzident_files<- vector("list")
  
  for( i in 1:length(listofMSruns)){
    
    MSfileTocheck<- gsub("_mzRaw.csv",".mzXML",listofMSruns[i])
    temp_pepXMLfile<- pepxml_oututfile[grep(basename(MSfileTocheck),pepxml_oututfile[,"MSRun"]),]
    Xrea_data<-read.csv(listofMSruns[i],header=T)
    mzidenML_csvfile_xrea<-merge(temp_pepXMLfile,Xrea_data,by="Scan",all=T)
    mzidenML_csvfile_xrea<-mzidenML_csvfile_xrea[mzidenML_csvfile_xrea$charge!=1,]
    mzidenML_csvfile_xrea$MZ<- (mzidenML_csvfile_xrea$MZ*mzidenML_csvfile_xrea$charge)- mzidenML_csvfile_xrea$charge
    mzident_files[[listofMSruns[i]]] <- (mzidenML_csvfile_xrea)
    
    
    
  }
  
  return(mzident_files)
  
}

CreatePeakList<- function(quantfilePath){
  
  options(scipen = 999)
  quant_inputfile<-read.csv(paste(paste(getwd(),quantfilePath,sep="/"),"/feature.csv",sep=""),header = T,stringsAsFactors = F)
  No.ofFiles<- length(grep(".m.z",colnames(quant_inputfile)))
  quantDataSets<- vector("list")
  quant_inputfile_mz<- quant_inputfile[,grep(".m.z",colnames(quant_inputfile))]
  quant_inputfile_rt<- quant_inputfile[,grep("^RT.mean$",colnames(quant_inputfile))]
  quant_inputfile_intensity<- quant_inputfile[,grep(".Normalized.Area",colnames(quant_inputfile))]
  quant_inputfile_pept<- quant_inputfile[,grep("Peptide",colnames(quant_inputfile))]
  quant_inputfile_z<- quant_inputfile[,grep("^z$",colnames(quant_inputfile))]
  
  
  
  for(i in 1:No.ofFiles){
    
    
    dataset<-data.frame(cbind(quant_inputfile_pept,quant_inputfile_z,
                              quant_inputfile_mz[,i],quant_inputfile_rt,quant_inputfile_intensity[,i]))
    colnames(dataset)<-paste(c("Peptide","Z","MZ","RT","Intensity"))
    
    dataset<- dataset[dataset$MZ != "-",]
    dataset<- dataset[dataset$Intensity != "-",]
    dataset$Intensity<- as.numeric(as.character(dataset$Intensity))
    dataset$RT<- as.numeric(as.character(dataset$RT))
    dataset$MZ<-as.numeric(as.character(dataset$MZ))*as.numeric(dataset$Z)-as.numeric(dataset$Z)
    dataset<- dataset[dataset$Z != 1,]
    
    quantDataSets[[i]]<- dataset
    
  }
  
  
  return(quantDataSets)
  
  
}

merged_pepXML_FeatureMatrix<-function(featureDataset,idDataset){
  
  ff_data<-do.call("rbind", lapply(featureDataset, data.frame))
  idDataset<-do.call("rbind", lapply(idDataset, data.frame))
  
  ff_csvfile_xrea_ms1<-ff_data[ff_data$Peptide == "",]
  ff_csvfile_xrea_ms1$Z<- NULL
  ff_csvfile_xrea_ms1$MZ<- round(ff_csvfile_xrea_ms1$MZ,2)
  ff_csvfile_xrea_ms1$Peptide<-NULL
  ff_csvfile_xrea_ms1<-ff_csvfile_xrea_ms1 %>% group_by(MZ) %>% summarise_each(funs(median))
  
  
  ff_csvfile_xrea_id<-merge(ff_data,idDataset,by="Peptide")
  ff_csvfile_xrea_id$MSRun<-NULL
  ff_csvfile_xrea_id$Precursor.MZ.<-NULL
  ff_csvfile_xrea_id$Charge<-NULL
  ff_csvfile_xrea_id$charge<-NULL
  ff_csvfile_xrea_id$MZ.y<-NULL
  ff_csvfile_xrea_id$RT.y<-NULL
  colnames(ff_csvfile_xrea_id)[3]<-paste(c("MZ"),sep="")
  colnames(ff_csvfile_xrea_id)[4]<-paste(c("RT"),sep="")
  ff_csvfile_xrea_id$ModStatus<-NULL
  ff_csvfile_xrea_id$Z<-NULL
  ff_csvfile_xrea_id$Scan<-NULL
  
  addCols<- length(ff_csvfile_xrea_id)-length(ff_csvfile_xrea_ms1)
  ff_csvfile_xrea_ms1<-cbind(ff_csvfile_xrea_ms1,matrix(nrow=nrow(ff_csvfile_xrea_ms1),ncol = addCols))
  colnames(ff_csvfile_xrea_ms1)<-paste(c("MZ","RT","Intensity","Xrea","Peptide","PID"),sep="")
  
  
  ff_csvfile_xrea_id$Peptide<- paste(ff_csvfile_xrea_id$Peptide,ff_csvfile_xrea_id$PID,sep="#")
  ff_csvfile_xrea_id$PID<-NULL
  
  ff_csvfile_xrea_id<-ff_csvfile_xrea_id %>% group_by(Peptide) %>% summarise_each(funs(median))
  
  ff_csvfile_xrea_id<- cSplit(ff_csvfile_xrea_id, "Peptide", sep = "#", direction = "wide")
  ff_csvfile_xrea_id$Peptide<- ff_csvfile_xrea_id$Peptide_1
  ff_csvfile_xrea_id$PID<- ff_csvfile_xrea_id$Peptide_2
  ff_csvfile_xrea_id<- ff_csvfile_xrea_id[,-c("Peptide_1","Peptide_2")]
  ff_csvfile_xrea_id<-ff_csvfile_xrea_id[, c("MZ","RT","Intensity","Xrea","Peptide","PID")]
  
  ff_csvfile_xrea<-rbind(ff_csvfile_xrea_ms1,ff_csvfile_xrea_id)
  return(ff_csvfile_xrea)
  
  
}

annotateMat<-function(out_tab1,out_tab2,mztol,rt_tol){
  
  dat1_ff = out_tab1[out_tab1$Peptide !='',]
  dat2_ff = out_tab2[out_tab2$Peptide !='',]
  
  
  commonids = merge(dat1_ff,dat2_ff,by='Peptide') 
  indx = which(abs(commonids$RT.y-commonids$RT.x) <=rt_tol)
  x = as.numeric(as.matrix(commonids[indx,"RT.x"]))
  y = as.numeric(as.matrix(commonids[indx,"RT.y"]))
  fit<-lm(y~x)
  
  commonids_comp =commonids[indx,]
  
  main_dat1 = out_tab1[out_tab1$Peptide !='',]
  main_dat2 = out_tab2[out_tab2$Peptide =='',]
  
  y_hat = main_dat1$RT*fit$coefficients[2]+fit$coefficients[1]
  
  
  ref_matrix1 = main_dat1
  sample_matrix1 = main_dat2 
  ref_matrix1$RT = y_hat
  
  x = NULL
  FP_count=0
  
  total <- dim(ref_matrix1)[1] ##Progress bar set to max rows
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for(i in 1:dim(ref_matrix1)[1]){
    
    xx_Dataset_indx_MZ_Dataset12<-sample_matrix1[which(as.numeric(as.character(sample_matrix1$MZ)) >= as.numeric(as.character(ref_matrix1$MZ[i]))-mztol & 
                                                         as.numeric(as.character(sample_matrix1$MZ)) <= as.numeric(as.character(ref_matrix1$MZ[i]))+mztol),]
    
    xx_Dataset_indx_RT_Dataset12<-xx_Dataset_indx_MZ_Dataset12[which(as.numeric(as.character(xx_Dataset_indx_MZ_Dataset12$RT)) >= as.numeric(as.character(ref_matrix1$RT[i]))-rt_tol & 
                                                                       as.numeric(as.character(xx_Dataset_indx_MZ_Dataset12$RT)) <= as.numeric(as.character(ref_matrix1$RT[i]))+rt_tol),]
    
    
    if(nrow(xx_Dataset_indx_RT_Dataset12)==0){
      xx_Dataset_indx_RT_Dataset12_temp<-data.frame(matrix(NA,1,ncol = length(dat2_ff)))
      names(xx_Dataset_indx_RT_Dataset12_temp)<-colnames(xx_Dataset_indx_RT_Dataset12)
      xx_Dataset_indx_RT_Dataset12<-xx_Dataset_indx_RT_Dataset12_temp
    }
    
    temp_dat = cbind(ref_matrix1[i,],xx_Dataset_indx_RT_Dataset12)
    
    x=rbind(x,temp_dat) 
    
    setTxtProgressBar(pb, i)
  }
  
  return(x)
  
}

ScatterXreaInt<- function(preData,DatasetName){
  
  
  xx<-data.frame(preData[,c("PrecursorIntensity","Status","Xrea.y")])
  colnames(xx)[3]<-"Xrea"
  
  options(scipen = 999)
  
  xx$PrecursorIntensity<- log10(as.numeric(as.character(xx$PrecursorIntensity)))
  xx$Xrea<- (as.numeric(as.character(xx$Xrea)))
  
  p1<-ggplot(xx, aes(Xrea,PrecursorIntensity)) + geom_point(aes(colour=Status),alpha=0.2)+
    scale_color_manual(values=c("darkseagreen", "blue", "red"))+
    guides(fill=FALSE)+
    theme(axis.title=element_text(size=60,face="bold"),axis.text=element_text(size=60))+
    scale_alpha(range = c(0.00, 1), guide = FALSE) +
    xlab("Xrea") +ylab("log10(PrecursorIntensity)")+
    ggtitle(DatasetName) + 
    theme_bw(base_size = 12)+
    theme(legend.position = "none")
  
  pTop <- ggplot(xx, aes(x = Xrea)) +
    geom_density(aes(colour=Status))+scale_color_manual(values=c("darkseagreen", "blue", "red"))+aes(y = ..count..)+theme_bw(base_size = 10)+
    theme(legend.position = "none")
  pRight <- ggplot(xx, aes(x = PrecursorIntensity)) +
    geom_density(aes(colour=Status))+scale_color_manual(values=c("darkseagreen", "blue", "red"))+aes(y = ..count..)+ coord_flip()+theme_bw(base_size = 10)
  
  pEmpty <- ggplot(xx, aes(Xrea,PrecursorIntensity)) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())
  pfinal<- grid.arrange(pTop, pEmpty, p1, pRight,
                        ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  return(pfinal)
}

XreaScatterPlot<- function(pepXMLfile,rawMZfilepath,FigTitle){
  
  pepXML_dataset<-do.call("rbind", lapply(pepXMLfile, data.frame))
  pepXML_dataset$ScanRawFile<- paste(pepXML_dataset$Scan,pepXML_dataset$MSRun,sep="_")
  
  xrea_dataset1_files<-list.files(rawMZfilepath, pattern = "_MS2data.csv", full.names = TRUE, ignore.case = TRUE)
  tables_xrea_dataset1 <- lapply(xrea_dataset1_files, read.csv, header = TRUE)
  combined.df_dataset1 <- do.call(rbind , tables_xrea_dataset1)
  combined.df_dataset1$ScanRawFile<- paste(combined.df_dataset1$Scan,combined.df_dataset1$Scan_rawfile,sep="_")
  
  pepXML_dataset1_append<- merge(pepXML_dataset,combined.df_dataset1,by="ScanRawFile",all=T)
  
  pepXML_dataset1_append[is.na(pepXML_dataset1_append$Peptide),"Status"]<- paste("1.Unidentified")
  pepXML_dataset1_append[!is.na(pepXML_dataset1_append$Peptide),"Status"]<- paste("2.Peaks (PSMs)")
  
  final_xrea_dataset1 <- pepXML_dataset1_append
  
  p<-ScatterXreaInt(final_xrea_dataset1,FigTitle)
  return(p)
  
}

ecdf_plot<-function(DF,plotTitle,mzFilePath){
  
  preData<- DF
  xrea_dataset1_files<-list.files(mzFilePath, pattern = "_MS2data.csv", full.names = TRUE, ignore.case = TRUE)
  tables_xrea_dataset1 <- lapply(xrea_dataset1_files, read.csv, header = TRUE)
  combined.df_dataset1 <- do.call(rbind , tables_xrea_dataset1)
  combined.df_dataset1$ScanMSfile<- paste(combined.df_dataset1$Scan,combined.df_dataset1$Scan_rawfile,sep="_")
  combined.df_dataset1<- combined.df_dataset1[!duplicated(combined.df_dataset1$ScanMSfile),]
  combined.df_dataset1<- combined.df_dataset1[combined.df_dataset1$charge != 1,]
  
  all_msdata<-as.numeric(format(preData$Intensity,scientific=F))
  all_msmsdata<-as.numeric(format(combined.df_dataset1[,"PrecursorIntensity"]))
  all_iddata<-as.numeric(format(preData[!is.na(preData$Peptide),"Intensity"]))
  
  MS1Data<- log10(all_msdata[!is.na(all_msdata)])
  MS2Data<- log10(all_msmsdata[!is.na(all_msmsdata)])
  MS2_id_Data<- log10(all_iddata[!is.na(all_iddata)])
  
  
  MS1Data<-append(MS1Data,MS2Data,MS2_id_Data)
  
  
  ecdf1 <- hist(MS1Data,breaks=seq(min(MS1Data),max(MS1Data),l=sqrt(length(MS1Data))+1),plot = F)
  multiplier1 <- ecdf1$counts / ecdf1$density
  MS1dens <- density(MS1Data)
  MS1dens$y <- MS1dens$y * multiplier1[1]
  
  
  ecdf2 <- hist(MS2Data,breaks=seq(min(MS2Data),max(MS2Data),l=sqrt(length(MS2Data))+1),plot = F)
  multiplier2 <- ecdf2$counts / ecdf2$density
  MS2dens <- density(MS2Data)
  MS2dens$y <- MS2dens$y * multiplier2[1]
  
  ecdf3 <- hist(MS2_id_Data,breaks=seq(min(MS2_id_Data),max(MS2_id_Data),l=sqrt(length(MS2_id_Data))+1),plot = F)
  multiplier3 <- ecdf3$counts / ecdf3$density
  MS2_id_dens <- density(MS2_id_Data)
  MS2_id_dens$y <- MS2_id_dens$y * multiplier3[1]
  
  
  # get the range for the x and y axis 
  xrange <- range(MS1dens$x) 
  yrange <- range(MS1dens$y) 
  
  
  # set up the plot 
  png(paste(plotTitle,"_IntensityDist.png",sep=""),width=900, height =900, res=50)
  
  plot(xrange, yrange, type="n", xlab="Log10(intensity)",main = plotTitle,
       ylab="Counts",cex.lab=2.5, cex.axis=2, cex.main=3) 
  
  
  p<-lines(MS1dens,col="red",ylab='',xlab="",main ="",lwd=5)
  p<-lines(MS2dens,col="green",ylab='',xlab="",main ="",lwd=5,add=T)
  p<-lines(MS2_id_dens,col="orange",ylab='',xlab="",main ="",lwd=5,add=T)
  
  
  legend("topleft", 
         legend = c("MS1", "MS2", "PSMs"), 
         col = c("red","green","orange"), ncol = 1,
         cex = 1.75,bty = "n",lwd=7)
  dev.off()
  return(p)
  
  
}