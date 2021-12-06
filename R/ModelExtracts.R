
Extract_LMER<-function(Mod, OutputFile="OutputFileName.csv"){
datOut = as.data.frame(summary(Mod)$coef)
datOut[,1:4]=round(datOut[,1:4], digits=2)
names(datOut)[names(datOut) == "Pr(>|t|)"] <- "p_value"
datOut$p_value=as.numeric(as.character(datOut$p_value))
names(datOut)[names(datOut) == "t value"] <- "t_value"
datOut<-tibble::rownames_to_column(datOut, "Fixed_Factors")
datOut$Significance[is.na(datOut$p_value)==F]="n.s."
datOut$Significance[datOut$p_value<.1]=" ."
datOut$Significance[datOut$p_value<.05]=" *"
datOut$Significance[datOut$p_value<.01]=" **"
datOut$Significance[datOut$p_value<.001]=" ***"
names(datOut)[names(datOut) == "t_value2"] <- "t value"
names(datOut)[names(datOut) == "Fixed_Factors"] <- "Fixed_Effects"
names(datOut)[names(datOut) == "Std. Error"] <- "SE"
datOut$p_value=format(round(datOut$p_value, digits=3), nsmall=3)
colnames(datOut)[colnames(datOut)=="p_value"]="p value"
colnames(datOut)[colnames(datOut)=="t_value"]="t value"
vc=VarCorr(Mod)
RanEfs=as.data.frame(vc,comp=c("Variance","Std.Dev."),digits=2)
nrows=nrow(datOut)
datOut[nrows+1,1]= "Significance levels: *** p < .001, ** p < .01, * p < .05, . p < .1"

RanEfsString0="Random Intercepts: "
RanEfs$dummyVar="ignore"
RanEfs$dummyVar=RanEfs$var2
RanInts=subset(RanEfs, var1 == "(Intercept)" & is.na(var2)==T)
nIntercepts=nrow(RanInts)
for (i in 1:nIntercepts){
RanEfsString1=RanInts$grp[i]
RanEfsString2= " (Var = "
RanEfsString3=format(round(RanEfs$vcov[i], digits=3), nsmall=3)
RanEfsString4=", SD = "
RanEfsString5=format(round(RanEfs$sdcor[i], digits=3), nsmall=3)
RanEfsString6=")"
if (i==1){
datOut[nrows+3, 1]=paste(RanEfsString0, RanEfsString1, RanEfsString2, RanEfsString3, RanEfsString4, RanEfsString5, RanEfsString6, sep ="")
}else{
datOut[nrows+3, 1]=paste(datOut[nrows+3, 1], "; ", RanEfsString1, RanEfsString2, RanEfsString3, RanEfsString4, RanEfsString5, RanEfsString6, sep ="")
}
}

RanSlopes = subset(RanEfs, var1 != "(Intercept)" & is.na(RanEfs$dummyVar)==T)
nRanSlopes=nrow(RanSlopes)
if (is.na(RanEfs$dummyVar[1])==T){
  for (i in 1: nRanSlopes){
RanSlopes0="Random Slopes: "
RanSlopes1= paste(RanSlopes$var1[i])
RanSlopes2 = " by "
RanSlopes3 = paste(RanSlopes$grp[i])
RanSlopes4= " (Var = "
RanSlopes5=format(round(RanSlopes$vcov[i], digits=3), nsmall=3)
RanSlopes6=", SD = "
RanSlopes7=format(round(RanSlopes$sdcor[i], digits=3), nsmall = 3)
RanSlopes8=")"
if (i==1){
  datOut[nrows+4, 1]=paste(RanSlopes0, RanSlopes1, RanSlopes2, RanSlopes3, RanSlopes4, RanSlopes5, RanSlopes6, RanSlopes7, RanSlopes8, sep ="")
}else{
  datOut[nrows+4, 1]=paste(datOut[nrows+4, 1], "; ", RanSlopes1, RanSlopes2, RanSlopes3, RanSlopes4, RanSlopes5, RanSlopes6, RanSlopes7, RanSlopes8, sep ="")
}
}}

datOut[is.na(datOut)==T]=""
print(datOut)
write.csv(datOut, file = OutputFile)
}



Extract_GLMER<-function(Mod, OutputFile="OutputFileName.csv"){
  datOut = as.data.frame(summary(Mod)$coef)
  datOut[,1:3]=round(datOut[,1:3], digits=2)
  names(datOut)[names(datOut) == "Pr(>|z|)"] <- "p_value"
  datOut$p_value=as.numeric(as.character(datOut$p_value))
  names(datOut)[names(datOut) == "z value"] <- "z_value"
  datOut<-tibble::rownames_to_column(datOut, "Fixed_Factors")
  datOut$Significance[is.na(datOut$p_value)==F]="n.s."
  datOut$Significance[datOut$p_value<.1]=" ."
  datOut$Significance[datOut$p_value<.05]=" *"
  datOut$Significance[datOut$p_value<.01]=" **"
  datOut$Significance[datOut$p_value<.001]=" ***"
  names(datOut)[names(datOut) == "t_value2"] <- "t value"
  names(datOut)[names(datOut) == "Fixed_Factors"] <- "Fixed_Effects"
  names(datOut)[names(datOut) == "Std. Error"] <- "SE"
  datOut$p_value=format(round(datOut$p_value, digits=3), nsmall=3)
  colnames(datOut)[colnames(datOut)=="p_value"]="p value"
  colnames(datOut)[colnames(datOut)=="t_value"]="t value"

  vc=VarCorr(Mod)
  RanEfs=as.data.frame(vc,comp=c("Variance","Std.Dev."),digits=2)
  nrows=nrow(datOut)
  datOut[nrows+1,1]= "Significance levels: *** p < .001, ** p < .01, * p < .05, . p < .1"

  RanEfsString0="Random Intercepts: "
  RanEfs$dummyVar="ignore"
  RanEfs$dummyVar=RanEfs$var2
  RanInts=subset(RanEfs, var1 == "(Intercept)" & is.na(var2)==T)
  nIntercepts=nrow(RanInts)
  for (i in 1:nIntercepts){
    RanEfsString1=RanInts$grp[i]
    RanEfsString2= " (Var = "
    RanEfsString3=format(round(RanEfs$vcov[i], digits=3), nsmall=3)
    RanEfsString4=", SD = "
    RanEfsString5=format(round(RanEfs$sdcor[i], digits=3), nsmall=3)
    RanEfsString6=")"
    if (i==1){
      datOut[nrows+3, 1]=paste(RanEfsString0, RanEfsString1, RanEfsString2, RanEfsString3, RanEfsString4, RanEfsString5, RanEfsString6, sep ="")
    }else{
      datOut[nrows+3, 1]=paste(datOut[nrows+3, 1],";", RanEfsString1, RanEfsString2, RanEfsString3, RanEfsString4, RanEfsString5, RanEfsString6, sep ="")
    }
  }

  RanSlopes = subset(RanEfs, var1 != "(Intercept)" & is.na(RanEfs$dummyVar)==T)
  nRanSlopes=nrow(RanSlopes)
  if (is.na(RanEfs$dummyVar[1])==T){
    for (i in 1: nRanSlopes){
      RanSlopes0="Random Slopes: "
      RanSlopes1= paste(RanSlopes$var1[i])
      RanSlopes2 = " by "
      RanSlopes3 = paste(RanSlopes$grp[i])
      RanSlopes4= " (Var = "
      RanSlopes5=format(round(RanSlopes$vcov[i], digits=3), nsmall=3)
      RanSlopes6=", SD = "
      RanSlopes7=format(round(RanSlopes$sdcor[i], digits=3), nsmall=3)
      RanSlopes8=")"
      if (i==1){
        datOut[nrows+4, 1]=paste(RanSlopes0, RanSlopes1, RanSlopes2, RanSlopes3, RanSlopes4, RanSlopes5, RanSlopes6, RanSlopes7, RanSlopes8, sep ="")
      }else{
        datOut[nrows+4, 1]=paste(datOut[nrows+4, 1], "; ", RanSlopes1, RanSlopes2, RanSlopes3, RanSlopes4, RanSlopes5, RanSlopes6, RanSlopes7, RanSlopes8, sep ="")
      }
    }}

  datOut[is.na(datOut)==T]=""
  print(datOut)
  write.csv(datOut, file = OutputFile)
}
