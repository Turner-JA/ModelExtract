Extract_LMER<-function(Mod, OutputFile="OutputFileName.csv",
                       showFormula=T,
                       showObs=T,
                       showRanInts =T,
                       showRanSlopes=T,
                       showPvalues=T){

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
  names(datOut)[names(datOut) == "Fixed_Factors"] <- "Fixed Effects"
  names(datOut)[names(datOut) == "Std. Error"] <- "SE"
  colnames(datOut)[colnames(datOut)=="t_value"]="t value"
  datOut$p_value=format(round(datOut$p_value, digits=3), nsmall=3)
  datOut$p_value[datOut$p_value=="0.000"]="< .0001"
  datOut$p_value[datOut$p_value=="1.000"]="> .999"
  if (showPvalues==F){
    datOut=subset(datOut, select=-c(p_value))
  }
  colnames(datOut)[colnames(datOut)=="p_value"]="p value"
  colnames(datOut)[colnames(datOut)=="(Intercept)"]="Intercept"
  nrows=nrow(datOut)
  datOut[nrows+1,1]= "Significance levels: *** p < .001, ** p < .01, * p < .05, . p < .1"
  datOut[nrows+2,1]= ""
  vc=VarCorr(Mod)
  RanEfs=as.data.frame(vc,comp=c("Variance","Std.Dev."),digits=2)
  obs.no=nobs(Mod)
  if(showFormula==T){
    formula=format(Mod@call[["formula"]])
    test=as.data.frame(formula)
    nrows2=nrow(test)
    if (nrows2 >1){
      for (i in 2:nrows2){
        test[i,]=str_trim(test[i,])
        formula= paste(formula[1], test[i,], sep = "")
      }
    }
    nrows=nrow(datOut)
    datOut[nrows+1, 1]=paste("Formula: ", formula, sep="")}
  if(showObs==T){
    nrows=nrow(datOut)
    datOut[nrows+1, 1]=paste("No. of observations: ", obs.no, sep="")}
  RanInts=subset(RanEfs, var1 == "(Intercept)" & is.na(var2)==T)
  if (nrow(RanInts)>0 & showRanInts == T){
    nIntercepts=nrow(RanInts)
    nrows=nrow(datOut)
    for (i in 1:nIntercepts){
      RanEfsString1=RanInts$grp[i]
      nobsGrp=(Mod@frame[RanEfsString1])
      nobsGrp= as.character(nlevels(nobsGrp[,1]))
      RanEfsString2=format(round(RanInts$vcov[i], digits=2), nsmall=2)
      RanEfsString3=format(round(RanInts$sdcor[i], digits=2), nsmall=2)
      if (i==1){
        datOut[nrows+1, 1]=paste("Random Intercepts: ", RanEfsString1, " (", nobsGrp, ", Var = ", RanEfsString2, ", SD = ", RanEfsString3, ")", sep ="")
      }else{
        datOut[nrows+1, 1]=paste(datOut[nrows+1, 1], "; ", RanEfsString1, " (", nobsGrp, ", Var = ", RanEfsString2, ", SD = ", RanEfsString3, ")", sep ="")
      }
    }
  }
  RanSlopes = subset(RanEfs, var1 != "(Intercept)" & is.na(var2)==T)
  if (nrow(RanSlopes)>0 & showRanSlopes==T){
    nRanSlopes=nrow(RanSlopes)
    nrows=nrow(datOut)
    for (i in 1: nRanSlopes){
      RanSlopes1= paste(RanSlopes$var1[i])
      RanSlopes2 = paste(RanSlopes$grp[i])
      RanSlopes3=format(round(RanSlopes$vcov[i], digits=2), nsmall=2)
      RanSlopes4=format(round(RanSlopes$sdcor[i], digits=2), nsmall = 2)
      if (i==1){
        datOut[nrows+1, 1]=paste("Random Slopes: ", RanSlopes1, " by ", RanSlopes2, " (Var = ", RanSlopes3, ", SD = ", RanSlopes4, ")", sep ="")
      }else{
        datOut[nrows+1, 1]=paste(datOut[nrows+1, 1], "; ", RanSlopes1, " by ", RanSlopes2, " (Var = ", RanSlopes3, ", SD = ", RanSlopes4, ")", sep ="")
      }
    }}
  datOut[is.na(datOut)==T]=""
  write.xlsx(datOut, file = OutputFile, rownames=F)
}


Extract_GLMER<-function(Mod, OutputFile="OutputFileName.csv",
                        showFormula=T,
                        showObs=T,
                        showRanInts =T,
                        showRanSlopes=T,
                        showPvalues=T){

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
  names(datOut)[names(datOut) == "Fixed_Factors"] <- "Fixed Effects"
  names(datOut)[names(datOut) == "Std. Error"] <- "SE"
  colnames(datOut)[colnames(datOut)=="t_value"]="t value"
  datOut$p_value=format(round(datOut$p_value, digits=3), nsmall=3)
  datOut$p_value[datOut$p_value=="0.000"]="< .0001"
  datOut$p_value[datOut$p_value=="1.000"]="> .999"
  if (showPvalues==F){
    datOut=subset(datOut, select=-c(p_value))
  }
  colnames(datOut)[colnames(datOut)=="p_value"]="p value"
  colnames(datOut)[colnames(datOut)=="(Intercept)"]="Intercept"
  nrows=nrow(datOut)
  datOut[nrows+1,1]= "Significance levels: *** p < .001, ** p < .01, * p < .05, . p < .1"
  datOut[nrows+2,1]= ""
  vc=VarCorr(Mod)
  RanEfs=as.data.frame(vc,comp=c("Variance","Std.Dev."),digits=2)
  obs.no=nobs(Mod)
  if(showFormula==T){
    formula=format(Mod@call[["formula"]])
    test=as.data.frame(formula)
    nrows2=nrow(test)
    if (nrows2 >1){
      for (i in 2:nrows2){
        test[i,]=str_trim(test[i,])
        formula= paste(formula[1], test[i,], sep = "")
      }
    }
    nrows=nrow(datOut)
    datOut[nrows+1, 1]=paste("Formula: ", formula, sep="")}
  if(showObs==T){
    nrows=nrow(datOut)
    datOut[nrows+1, 1]=paste("No. of observations: ", obs.no, sep="")}
  RanInts=subset(RanEfs, var1 == "(Intercept)" & is.na(var2)==T)
  if (nrow(RanInts)>0 & showRanInts == T){
    nIntercepts=nrow(RanInts)
    nrows=nrow(datOut)
    for (i in 1:nIntercepts){
      RanEfsString1=RanInts$grp[i]
      nobsGrp=(Mod@frame[RanEfsString1])
      nobsGrp= as.character(nlevels(nobsGrp[,1]))
      RanEfsString2=format(round(RanInts$vcov[i], digits=2), nsmall=2)
      RanEfsString3=format(round(RanInts$sdcor[i], digits=2), nsmall=2)
      if (i==1){
        datOut[nrows+1, 1]=paste("Random Intercepts: ", RanEfsString1, " (", nobsGrp, ", Var = ", RanEfsString2, ", SD = ", RanEfsString3, ")", sep ="")
      }else{
        datOut[nrows+1, 1]=paste(datOut[nrows+1, 1], "; ", RanEfsString1, " (", nobsGrp, ", Var = ", RanEfsString2, ", SD = ", RanEfsString3, ")", sep ="")
      }
    }
  }
  RanSlopes = subset(RanEfs, var1 != "(Intercept)" & is.na(var2)==T)
  if (nrow(RanSlopes)>0 & showRanSlopes==T){
    nRanSlopes=nrow(RanSlopes)
    nrows=nrow(datOut)
    for (i in 1: nRanSlopes){
      RanSlopes1= paste(RanSlopes$var1[i])
      RanSlopes2 = paste(RanSlopes$grp[i])
      RanSlopes3=format(round(RanSlopes$vcov[i], digits=2), nsmall=2)
      RanSlopes4=format(round(RanSlopes$sdcor[i], digits=2), nsmall = 2)
      if (i==1){
        datOut[nrows+1, 1]=paste("Random Slopes: ", RanSlopes1, " by ", RanSlopes2, " (Var = ", RanSlopes3, ", SD = ", RanSlopes4, ")", sep ="")
      }else{
        datOut[nrows+1, 1]=paste(datOut[nrows+1, 1], "; ", RanSlopes1, " by ", RanSlopes2, " (Var = ", RanSlopes3, ", SD = ", RanSlopes4, ")", sep ="")
      }
    }}
  datOut[is.na(datOut)==T]=""
  write.xlsx(datOut, file = OutputFile, rownames=F)
}


effects_code<-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]
    df[,var]=as.factor(df[,var])
    df[,var]=droplevels(df[,var])
    test=contrasts(df[,var])
    ncats=nrow(test)
    my.coding<-matrix(rep(1/ncats, (ncats*(ncats-1))), ncol=ncats-1)
    contrasts(df[,var])<-test-my.coding
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[,var]))
  }
}

dummy_code<-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]   
    df[,var]=as.factor(df[,var])
    df[,var]=droplevels(df[,var])
    test=contrasts(df[,var])
    ncats=nrow(test)
    collabels=rownames(test)
    collabels=collabels[2:ncats]
    my.coding=contrasts(df[,var])
    my.coding=contr.treatment(ncats)
    colnames(my.coding)=collabels
    contrasts(df[,var])=my.coding
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[,var]))
  }
}

deviation_code<-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]
    df[,var]=as.factor(df[,var])
    df[,var]=droplevels(df[,var])
    test=contrasts(df[,var])
    ncats=nrow(test)
    collabels=rownames(test)
    collabels=collabels[1:ncats-1]
    my.coding=contrasts(df[,var])
    my.coding=contr.sum(ncats)
    colnames(my.coding)=collabels
    contrasts(df[,var])=my.coding
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[,var]))
  }
}

deviation_code2 <-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]
    df[,var]=as.factor(df[,var])
    df[,var]=droplevels(df[,var])
    test=contrasts(df[,var])
    ncats=nrow(test)
    collabels=rownames(test)
    finalcollabel=collabels[ncats]
    df[,var]<-relevel(df[,var], ref=finalcollabel)
    test=contrasts(df[,var])
    ncats=nrow(test)
    collabels=rownames(test)
    finalcollabel=collabels[ncats]
    collabels=collabels[1:ncats-1]
    my.coding=contrasts(df[,var])
    my.coding=contr.sum(ncats)
    colnames(my.coding)=collabels
    contrasts(df[,var])=my.coding
    contrasts(df[,var])
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[,var]))
  }
}
