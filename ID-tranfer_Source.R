# Perform ID transfer ##################################################################

df<- data.frame()
df<- annotateMat(InputParam$RefSet,InputParam$SampSet,InputParam$mzTH,InputParam$rtTH)
assign(InputParam$IDtransferFileName,df)
df<- NULL

print(paste("Finished ID transfer"))

