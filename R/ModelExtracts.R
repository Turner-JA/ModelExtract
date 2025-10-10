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
  write.xlsx(datOut, file = OutputFile, row.names=F)
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
  write.xlsx(datOut, file = OutputFile, row.names=F)
}


effects_code<-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]
    df[[var]]=as.factor(df[[var]])
    df[[var]]=droplevels(df[[var]])
    test=contrasts(df[[var]])
    ncats=nrow(test)
    my.coding<-matrix(rep(1/ncats, (ncats*(ncats-1))), ncol=ncats-1)
    contrasts(df[[var]])<-test-my.coding
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[[var]]))
  }
}

dummy_code<-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]   
    df[[var]]=as.factor(df[[var]])
    df[[var]]=droplevels(df[[var]])
    test=contrasts(df[[var]])
    ncats=nrow(test)
    collabels=rownames(test)
    collabels=collabels[2:ncats]
    my.coding=contrasts(df[[var]])
    my.coding=contr.treatment(ncats)
    colnames(my.coding)=collabels
    contrasts(df[[var]])=my.coding
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[[var]]))
  }
}

deviation_code<-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]
    df[[var]]=as.factor(df[[var]])
    df[[var]]=droplevels(df[[var]])
    test=contrasts(df[[var]])
    ncats=nrow(test)
    collabels=rownames(test)
    collabels=collabels[1:ncats-1]
    my.coding=contrasts(df[[var]])
    my.coding=contr.sum(ncats)
    colnames(my.coding)=collabels
    contrasts(df[[var]])=my.coding
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[[var]]))
  }
}

deviation_code2 <-function(df, vars){
  dfname=substitute(df)
  nvars=length(vars)
  for (i in 1:nvars){
    var <-vars[[i]]
    df[[var]]=as.factor(df[[var]])
    df[[var]]=droplevels(df[[var]])
    test=contrasts(df[[var]])
    ncats=nrow(test)
    collabels=rownames(test)
    finalcollabel=collabels[ncats]
    df[[var]]<-relevel(df[[var]], ref=finalcollabel)
    test=contrasts(df[[var]])
    ncats=nrow(test)
    collabels=rownames(test)
    finalcollabel=collabels[ncats]
    collabels=collabels[1:ncats-1]
    my.coding=contrasts(df[[var]])
    my.coding=contr.sum(ncats)
    colnames(my.coding)=collabels
    contrasts(df[[var]])=my.coding
    contrasts(df[[var]])
    assign(paste(dfname),df,envir=.GlobalEnv)
    print(contrasts(df[[var]]))
  }
}


RanSlope_Tester <- function(DF, var, RanIntercepts) {
  # Check if `var` contains "*"; if so, split it
  split_vars <- if (grepl("\\*", var)) unlist(strsplit(var, "\\*")) else var
  
  for (RanIntercept in RanIntercepts) {
    # Construct the formula dynamically based on whether `var` was split
    formula_parts <- c(RanIntercept, split_vars)
    formula_str <- paste("~", paste(sprintf("DF[['%s']]", formula_parts), collapse = " + "))
    test <- as.data.frame(xtabs(as.formula(formula_str), data = DF))
    test <- subset(test, Freq == 0)
    
    # Output colored messages based on `test` result
    if (nrow(test) == 0) {
      message(crayon::green("Random slope of ", var, " justified for ", RanIntercept))
    } else {
      message(crayon::red("Random slope of ", var, " NOT justified for ", RanIntercept, ": lack of variance"))
    }
  }
}

RanSlope_Tester2 <- function(DF, var, RanIntercepts) {
  # Check if `var` contains "*"; if so, split it
  split_vars <- if (grepl("\\*", var)) unlist(strsplit(var, "\\*")) else var
  
  for (RanIntercept in RanIntercepts) {
    # Construct the formula dynamically based on whether `var` was split
    formula_parts <- c(RanIntercept, split_vars)
    formula_str <- paste("~", paste(sprintf("DF[['%s']]", formula_parts), collapse = " + "))
    test <- as.data.frame(xtabs(as.formula(formula_str), data = DF))
    test0s <- subset(test, Freq == 0)
    test0sand1s <- subset(test, Freq != 0 & Freq !=1)
    
    # Output colored messages based on `test` result
    if (nrow(test0s) == 0 & nrow(test0sand1s) !=0) {
      message(crayon::green("Random slope of ", var, " justified for ", RanIntercept))
    } else {
      message(crayon::red("Random slope of ", var, " NOT justified for ", RanIntercept, ": lack of variance"))
    }
  }
}

get_min_ranef_slope <- function(model) {
  vc <- VarCorr(model)$cond
  
  # Collect SDs per grouping factor
  ranef_sd_df <- do.call(rbind, lapply(names(vc), function(g) {
    df <- as.data.frame(attr(vc[[g]], "stddev"))
    
    # If only intercept, skip
    if(nrow(df) <= 1) return(NULL)
    
    # Keep slopes only (drop intercept)
    slopes_df <- df[-1,, drop=FALSE]
    
    # Data frame with all slopes per factor
    data.frame(
      term = rownames(slopes_df),
      stddev = slopes_df[,1],
      group = g,
      row.names = NULL
    )
  }))
  
  # If no slopes exist at all
  if(nrow(ranef_sd_df) == 0) {
    return(list(
      ranef_sd_df = NULL,
      min_ranef_factor = NULL,
      suggestion = "No random slopes in the model to remove."
    ))
  }
  
  # Compute average SD per grouping factor
  avg_sd_per_group <- aggregate(stddev ~ group, data = ranef_sd_df, FUN = mean)
  
  # Identify the factor with the smallest average SD
  min_factor <- avg_sd_per_group[which.min(avg_sd_per_group$stddev), ]
  
  suggestion <- paste0("Suggest removing random slopes for group '", 
                       min_factor$group, "' (average stddev = ", 
                       round(min_factor$stddev, 4), ")")
  
  list(
    ranef_sd_df = ranef_sd_df,       # all slopes
    avg_sd_per_group = avg_sd_per_group,  # average SD per factor
    min_ranef_factor = min_factor,
    suggestion = suggestion
  )
}

RanSlope_Tester_Final <- function(DF, dv, var, RanIntercepts, min_prop = 0.3, min_cluster_n = 5, return_table = FALSE) {
    
    # --- Load required packages ---
    if (!requireNamespace("dplyr", quietly = TRUE) ||
        !requireNamespace("glue", quietly = TRUE) ||
        !requireNamespace("crayon", quietly = TRUE) ||
        !requireNamespace("scales", quietly = TRUE) ||
        !requireNamespace("rlang", quietly = TRUE)) {
        stop("Packages dplyr, glue, crayon, scales, and rlang are required.")
    }
    
    # --- Internal function to check random slopes ---
    RanSlope_Tester12 <- function(DF, dv, var, RanIntercepts, min_prop = 0.3, min_cluster_n = 5, return_table = FALSE) {
        
        # Parse and evaluate var expression safely
        var_expr <- rlang::parse_expr(var)
        evaluated_var <- with(DF, eval(var_expr))
        
        # Add evaluated_var as a temp column to DF (use unique safe name)
        temp_var_name <- ".temp_var_for_ranslope_check"
        DF[[temp_var_name]] <- evaluated_var
        
        var_is_continuous <- is.numeric(DF[[temp_var_name]])
        
        check_dv_variance <- function(dv_vector) {
            if (is.numeric(dv_vector)) {
                return(var(dv_vector, na.rm = TRUE) > 0)
            } else {
                present_levels <- unique(dv_vector[!is.na(dv_vector)])
                return(length(present_levels) > 1)
            }
        }
        
        results <- list()
        
        for (RanIntercept in RanIntercepts) {
            message(crayon::blue("\nChecking random slope for predictor: ", var, " within grouping factor: ", RanIntercept))
            
            cluster_sizes <- DF %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
                dplyr::summarise(n = dplyr::n(), .groups = "drop")
            
            small_clusters <- sum(cluster_sizes$n < min_cluster_n, na.rm = TRUE)
            total_clusters <- nrow(cluster_sizes)
            prop_small_clusters <- small_clusters / total_clusters
            
            if (small_clusters > 0) {
                message(crayon::yellow(glue::glue("⚠️ {small_clusters} out of {total_clusters} groups have fewer than {min_cluster_n} observations.")))
            }
            
            variation_result <- NA
            prop_passing <- NA
            prop_unbalanced <- 0
            
            if (var_is_continuous) {
                overall_sd <- stats::sd(DF[[temp_var_name]], na.rm = TRUE)
                
                variation_check <- DF %>%
                    dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
                    dplyr::summarise(
                        unique_vals = dplyr::n_distinct(.data[[temp_var_name]]),
                        sd_val = stats::sd(.data[[temp_var_name]], na.rm = TRUE),
                        .groups = "drop"
                    ) %>%
                    dplyr::mutate(
                        meaningful_variation = ifelse(is.na(sd_val), FALSE, sd_val > (overall_sd * 0.05))
                    )
                
                n_groups <- nrow(variation_check)
                n_with_variation <- sum(variation_check$meaningful_variation, na.rm = TRUE)
                prop_with_variation <- n_with_variation / n_groups
                prop_passing <- prop_with_variation
                
                message_text <- glue::glue("{n_with_variation} out of {n_groups} ({scales::percent(prop_with_variation)}) groups show meaningful within-cluster variation in {var}.")
                
                if (n_with_variation == 0) {
                    variation_result <- "Impossible"
                    message(crayon::red(message_text))
                } else if (prop_with_variation >= min_prop) {
                    variation_result <- "Recommended"
                    message(crayon::green(message_text))
                } else {
                    variation_result <- "Low Variation"
                    message(crayon::yellow(message_text))
                }
                
            } else {
                counts <- DF %>%
                    dplyr::group_by(dplyr::across(dplyr::all_of(c(RanIntercept, temp_var_name)))) %>%
                    dplyr::summarise(n = dplyr::n(), .groups = "drop")
                
                replicated <- counts %>%
                    dplyr::filter(n >= 2)
                
                replicated_df <- DF %>%
                    dplyr::semi_join(replicated, by = c(RanIntercept, temp_var_name))
                
                dv_variance_check <- replicated_df %>%
                    dplyr::group_by(dplyr::across(dplyr::all_of(c(RanIntercept, temp_var_name)))) %>%
                    dplyr::summarise(has_variance = check_dv_variance(.data[[dv]]), .groups = "drop")
                
                valid_combos <- dv_variance_check %>%
                    dplyr::filter(has_variance)
                
                levels_per_group <- valid_combos %>%
                    dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
                    dplyr::summarise(levels_with_variance = dplyr::n_distinct(.data[[temp_var_name]]), .groups = "drop")
                
                n_groups <- nrow(cluster_sizes)
                n_multilevel <- sum(levels_per_group$levels_with_variance >= 2, na.rm = TRUE)
                prop_multilevel <- n_multilevel / n_groups
                prop_passing <- prop_multilevel
                
                level_props <- DF %>%
                    dplyr::group_by(dplyr::across(dplyr::all_of(c(RanIntercept, temp_var_name)))) %>%
                    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
                    dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
                    dplyr::mutate(prop = n / sum(n)) %>%
                    dplyr::ungroup()
                
                unbalanced_groups <- level_props %>%
                    dplyr::filter(prop < 0.1) %>%
                    dplyr::distinct(dplyr::across(dplyr::all_of(RanIntercept)))
                
                prop_unbalanced <- nrow(unbalanced_groups) / n_groups
                
                if (nrow(unbalanced_groups) > 0) {
                    message(crayon::yellow(glue::glue("⚠️ {nrow(unbalanced_groups)} groups have highly unbalanced category representation (<10% in at least one category).")))
                }
                
                message_text <- glue::glue("{n_multilevel} out of {n_groups} ({scales::percent(prop_multilevel)}) groups have ≥2 levels of {var} with replication and DV variance.")
                
                if (n_multilevel == 0) {
                    variation_result <- "Impossible"
                    message(crayon::red(message_text))
                } else if (prop_multilevel >= min_prop) {
                    variation_result <- "Recommended"
                    message(crayon::green(message_text))
                } else {
                    variation_result <- "Low Variation"
                    message(crayon::yellow(message_text))
                }
            }
            
            small_cluster_thresh <- 0.3
            imbalance_thresh <- 0.3
            
            if (variation_result == "Impossible") {
                overall_recommendation <- "Impossible"
            } else if (prop_small_clusters > small_cluster_thresh || prop_unbalanced > imbalance_thresh || (variation_result == "Low Variation" && prop_passing < (min_prop / 2))) {
                overall_recommendation <- "High Risk"
            } else if (variation_result == "Low Variation") {
                overall_recommendation <- "Low Confidence"
            } else {
                overall_recommendation <- "Recommended"
            }
            
            results[[RanIntercept]] <- list(
                Grouping_Factor = RanIntercept,
                Var_Type = ifelse(var_is_continuous, "Continuous", "Categorical"),
                Total_Groups = total_clusters,
                Problematic_Small_Groups = small_clusters,
                Prop_Small_Clusters = prop_small_clusters,
                Prop_Unbalanced_Groups = prop_unbalanced,
                Proportion_Passing = prop_passing,
                Variation_Result = variation_result,
                Overall_Recommendation = overall_recommendation
            )
        }
        
        summary_table <- dplyr::bind_rows(lapply(results, as.data.frame)) %>%
            dplyr::mutate(
                Overall_Recommendation = factor(
                    Overall_Recommendation,
                    levels = c("Impossible", "High Risk", "Low Confidence", "Recommended")
                )
            )
        
        if (return_table) {
            return(summary_table)
        } else {
            invisible(summary_table)
        }
    }
    
    # --- Helper to parse main effects from interaction ---
    parse_main_effects <- function(var) {
        parts <- unlist(strsplit(var, "\\*"))
        parts <- trimws(parts)
        return(parts)
    }
    
    # --- Main function logic ---
    is_interaction <- grepl("\\*", var)
    main_effects <- if (is_interaction) parse_main_effects(var) else var
    
    # Run diagnostics on main effects
    main_effect_results <- lapply(main_effects, function(mv) {
        RanSlope_Tester12(DF=DF, dv=dv, var=mv, RanIntercepts=RanIntercepts, min_prop=min_prop, min_cluster_n=min_cluster_n, return_table=TRUE)
    })
    
    main_effect_summary <- dplyr::bind_rows(main_effect_results, .id = "MainEffect")
    
    # Run interaction check if applicable
    interaction_summary <- NULL
    if (is_interaction) {
        interaction_summary <- RanSlope_Tester12(DF=DF, dv=dv, var=var, RanIntercepts=RanIntercepts, min_prop=min_prop, min_cluster_n=min_cluster_n, return_table=TRUE)
    }
    
    combined_summary <- main_effect_summary %>% 
        dplyr::mutate(Var_Type = "Main Effect")
    
    if (!is.null(interaction_summary)) {
        interaction_summary <- interaction_summary %>% 
            dplyr::mutate(Var_Type = "Interaction")
        combined_summary <- dplyr::bind_rows(combined_summary, interaction_summary)
    }
    
    # Check warnings for interaction if any main effect is bad
    warnings <- combined_summary %>%
        dplyr::filter(Var_Type == "Main Effect" & Overall_Recommendation %in% c("Impossible", "High Risk")) %>%
        dplyr::select(MainEffect, Grouping_Factor, Overall_Recommendation) %>%
        dplyr::distinct()
    
    if (nrow(warnings) > 0 & is_interaction) {
        warning_msg <- paste0(
            "⚠️ Interaction random slope '", var, "' may NOT be viable because main effect(s) ",
            paste(unique(warnings$MainEffect), collapse = ", "),
            " have low viability (", paste(unique(warnings$Overall_Recommendation), collapse = ", "), ") for grouping factor(s) ",
            paste(unique(warnings$Grouping_Factor), collapse = ", "), "."
        )
        warning(warning_msg)
    }
    
    # Print combined summary
    print(combined_summary)
    
    if (return_table) {
        return(combined_summary)
    } else {
        invisible(combined_summary)
    }
}


