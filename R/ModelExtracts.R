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
  print(datOut)
  write.csv(datOut, file = OutputFile)
}
