merge_lev <- function(
    df1,
    df2,
    by,
    threshold = 0.9,
    all.x = TRUE,
    show_DF2_by_col = TRUE,
    show_merging_similarity = TRUE
) {

  # Ensure by exists
  stopifnot(all(by %in% names(df1)))
  stopifnot(all(by %in% names(df2)))

  # Convert single column to vector
  if (length(by) == 1) {

    x <- as.character(df1[[by]])
    y <- as.character(df2[[by]])

    dist_matrix <- stringdist::stringdistmatrix(x, y, method = "lv")
    max_len <- outer(nchar(x), nchar(y), pmax)
    similarity_matrix <- 1 - (dist_matrix / max_len)

  } else {

    # --- MULTI COLUMN MODE ---

    similarity_list <- lapply(by, function(col) {

      x <- as.character(df1[[col]])
      y <- as.character(df2[[col]])

      dist_matrix <- stringdist::stringdistmatrix(x, y, method = "lv")
      max_len <- outer(nchar(x), nchar(y), pmax)

      1 - (dist_matrix / max_len)
    })

    # Average similarity across columns
    similarity_matrix <- Reduce("+", similarity_list) / length(similarity_list)
  }

  # Best match index per row
  best_match <- apply(similarity_matrix, 1, function(row) {
    if (all(is.na(row))) return(NA)
    idx <- which.max(row)
    if (row[idx] >= threshold) idx else NA
  })

  # Extract similarity score
  best_similarity <- sapply(seq_along(best_match), function(i) {
    if (!is.na(best_match[i])) similarity_matrix[i, best_match[i]] else NA
  })

  # Rename df2 merge columns
  df2_renamed <- df2

  by2_name <- paste0(by, "_2ndDF")
  names(df2_renamed)[match(by, names(df2_renamed))] <- by2_name

  matched_df2 <- df2_renamed[best_match, , drop = FALSE]

  # Optionally remove df2 merge columns
  if (!show_DF2_by_col) {
    matched_df2 <- matched_df2[, !(names(matched_df2) %in% by2_name), drop = FALSE]
  }

  # Combine results
  result <- dplyr::bind_cols(df1, matched_df2)

  # Add similarity column
  if (show_merging_similarity) {
    result$String_merge_similarity <- best_similarity
  }

  # all.x behaviour
  if (!all.x) {
    keep <- !is.na(best_match)
    result <- result[keep, , drop = FALSE]
  }

  return(result)
}


return_hidden_level_LMER <- function(model, factor_name, transform = NULL) {
library(brms)

  get_terms_chr <- function(formula) {
    tf <- terms(formula)
    fixed_terms <- attr(tf, "term.labels")
    rhs <- deparse(formula[[3]])
    rand_terms <- gregexpr("\\([^\\)]+\\|[^\\)]+\\)", rhs, perl = TRUE)
    rand_matches <- regmatches(rhs, rand_terms)[[1]]
    rand_terms_clean <- gsub("^\\(|\\)$", "", rand_matches)
    all_terms <- c(fixed_terms, rand_terms_clean)
    return(all_terms)
  }

  # Helper functions
  get_nonref <- function(levels, is_dev) {
    if (is_dev) {
      levels[-length(levels)]
    } else {
      levels[-1]
    }
  }
  
    get_reference_level <- function(model, factor_name, tol = 1e-8) {
    if (is.brmsfit(model)==T){
      f <- brms_model[["data"]][[factor_name]] 
    } else if (isLMM(model)==T){
      f <- model@frame[[factor_name]]
    }
    if (!is.factor(f)) stop(factor_name, " is not a factor")
    levels_factor <- levels(f)
    contr <- attr(f, "contrasts")
    if (is.null(contr)) return(levels_factor[1])
    zero_rows <- apply(contr, 1, function(x) all(abs(x) < tol))
    if (any(zero_rows)) return(levels_factor[which(zero_rows)[1]])
    row_scores <- rowSums(contr)
    levels_factor[which.min(row_scores)]
  }
  
# Default to fixef and vcov for other models
coefs <- fixef(model)
coefs_df <- as.data.frame(coefs)
V <- vcov(model, dpar = "mu")

  
  # Extract contrasts for the factor to find hidden reference level
  test <- as.data.frame(attr(model.frame(model)[[factor_name]], "contrasts"))
  cat <- tail(rownames(test), 1)  # Reference category label
  
  # Target coefficient names for non-reference levels of the factor
  target_coef_names <- paste0(factor_name, names(test))
  var_idx <- which(rownames(coefs_df) %in% target_coef_names)
  
  # Build contrast vector L for main effect hidden level
  L <- numeric(length(fixef(model)))
  names(L) <- names(fixef(model))
  L[var_idx] <- -1
  ct <- contest1D(model, L)
  
  # Main effect hidden level data frame
  main_effect_df <- data.frame(
    Term = paste0(factor_name, cat),
    Estimate = ct$Estimate,
    `Std. Error` = ct$`Std. Error`,
    `t value` = ct$`t value`,
    `Pr(>|t|)` = ct$`Pr(>|t|)`,
    check.names = FALSE
  )
  
  # Get original formula and term labels
  formula <- formula(model)
  terms_chr <- attr(terms(formula), "term.labels")
  
  # Find interactions involving factor_name
  relevant_interactions <- terms_chr[grepl(factor_name, terms_chr) & grepl(":", terms_chr)]
  
  # Separate two-way and three-way interactions
  two_way <- relevant_interactions[sapply(strsplit(relevant_interactions, ":"), length) == 2]
  three_way <- relevant_interactions[sapply(strsplit(relevant_interactions, ":"), length) == 3]
  
  # Prepare containers for interaction dfs
  interaction_two_way_dfs <- list()
  interaction_three_way_dfs <- list()
  
  # --- Process two-way interactions ---
  for (interaction in two_way) {
    factors <- unlist(strsplit(interaction, ":"))
    factorA <- factors[1]
    factorB <- factors[2]
    
    devA <- (factorA == factor_name)
    devB <- (factorB == factor_name)
    
    levelsA <- levels(model@frame[[factorA]])
    levelsB <- levels(model@frame[[factorB]])
    
    nonrefA <- get_nonref(levelsA, devA)
    nonrefB <- get_nonref(levelsB, devB)
    
    refA <- get_reference_level(model, factorA)
    refB <- get_reference_level(model, factorB)
    
    # Construct hidden label only when factor_name is first factor (devA)
    if (devA) {
      hidden_label <- paste0(
        factorA, refA, ":",
        factorB, nonrefB
      )
    } else if (devB) {
      # In case factor_name is second factor, adjust accordingly
      hidden_label <- paste0(
        factorA, nonrefA, ":",
        factorB, refB
      )
    } else {
      # If factor_name not in first or second factor, skip
      next
    }
    
    rn <- rownames(coefs_df)
    
    for (k in seq_along(hidden_label)) {
      lab <- hidden_label[k]
      
      var_idx <- which(
        sapply(strsplit(rn, ":"), length) == 2 &   # exactly two-way
          grepl(factorA, rn, fixed = TRUE) &
          grepl(factorB, rn, fixed = TRUE) &
          grepl(nonrefB[k], rn, fixed = TRUE)
      )
      if (length(var_idx) == 0) next
      
      L <- numeric(length(fixef(model)))
      names(L) <- names(fixef(model))
      L[var_idx] <- -1
      ct <- contest1D(model, L)
      
      df <- data.frame(
        Term = lab,
        Estimate = ct$Estimate,
        `Std. Error` = ct$`Std. Error`,
        `t value` = ct$`t value`,
        `Pr(>|t|)` = ct$`Pr(>|t|)`,
        check.names = FALSE
      )
      
      interaction_two_way_dfs[[length(interaction_two_way_dfs) + 1]] <- df
    }
  }
  
  # --- Process three-way interactions ---
  for (interaction in three_way) {
    factors <- unlist(strsplit(interaction, ":"))
    factorA <- factors[1]
    factorB <- factors[2]
    factorC <- factors[3]
    
    devA <- (factorA == factor_name)
    devB <- (factorB == factor_name)
    devC <- (factorC == factor_name)
    
    levelsA <- levels(model@frame[[factorA]])
    levelsB <- levels(model@frame[[factorB]])
    levelsC <- levels(model@frame[[factorC]])
    
    nonrefA <- get_nonref(levelsA, devA)
    nonrefB <- get_nonref(levelsB, devB)
    nonrefC <- get_nonref(levelsC, devC)
    
    refA <- get_reference_level(model, factorA)
    refB <- get_reference_level(model, factorB)
    refC <- get_reference_level(model, factorC)
    
    # Construct hidden label depending on which factor is factor_name
    if (devA) {
      # factor_name is first
      hidden_label <- paste0(
        factorA, refA, ":",
        factorB, nonrefB, ":",
        factorC, nonrefC
      )
    } else if (devB) {
      # factor_name is second
      hidden_label <- paste0(
        factorA, nonrefA, ":",
        factorB, refB, ":",
        factorC, nonrefC
      )
    } else if (devC) {
      # factor_name is third
      hidden_label <- paste0(
        factorA, nonrefA, ":",
        factorB, nonrefB, ":",
        factorC, refC
      )
    } else {
      next
    }
    
    rn <- names(coefs)
    
    for (k in seq_along(hidden_label)) {
      lab <- hidden_label[k]
      
      var_idx <- which(
        sapply(strsplit(rn, ":"), length) == 3 &   # exactly three-way
          grepl(factorA, rn, fixed = TRUE) &
          grepl(factorB, rn, fixed = TRUE) &
          grepl(factorC, rn, fixed = TRUE) &
          grepl(nonrefB[k], rn, fixed = TRUE)
      )
      if (length(var_idx) == 0) next
      
      L <- numeric(length(fixef(model)))
      names(L) <- names(fixef(model))
      L[var_idx] <- -1
      ct <- contest1D(model, L)
      
      df <- data.frame(
        Term = lab,
        Estimate = ct$Estimate,
        `Std. Error` = ct$`Std. Error`,
        `t value` = ct$`t value`,
        `Pr(>|t|)` = ct$`Pr(>|t|)`,
        check.names = FALSE
      )
      
      interaction_three_way_dfs[[length(interaction_three_way_dfs) + 1]] <- df
    }
  }
  
  # Combine all results in order: main → two-way → three-way
  final_df <- rbind(
    main_effect_df,
    if (length(interaction_two_way_dfs) > 0) do.call(rbind, interaction_two_way_dfs) else NULL,
    if (length(interaction_three_way_dfs) > 0) do.call(rbind, interaction_three_way_dfs) else NULL
  )
  
  return(final_df)
}


return_hidden_level_BRMS <- function(model, factor_name, transform) {
  library(brms) 
  
  ## Helpers ##
  get_nonref <- function(levels, is_dev) {
    if (is_dev) {
      levels[-length(levels)]
    } else {
      levels[-1]
    }
  }

  get_reference_level <- function(model, factor_name, tol = 1e-8) {
    if (is.brmsfit(model)==T){
      f <- brms_model[["data"]][[factor_name]] 
    } else if (isLMM(model)==T){
      f <- model@frame[[factor_name]]
    }
    if (!is.factor(f)) stop(factor_name, " is not a factor")
    levels_factor <- levels(f)
    contr <- attr(f, "contrasts")
    if (is.null(contr)) return(levels_factor[1])
    zero_rows <- apply(contr, 1, function(x) all(abs(x) < tol))
    if (any(zero_rows)) return(levels_factor[which(zero_rows)[1]])
    row_scores <- rowSums(contr)
    levels_factor[which.min(row_scores)]
  }
  
  get_terms_chr <- function(formula) {
    tf <- terms(formula)
    fixed_terms <- attr(tf, "term.labels")
    rhs <- deparse(formula[[3]])
    rand_terms <- gregexpr("\\([^\\)]+\\|[^\\)]+\\)", rhs, perl = TRUE)
    rand_matches <- regmatches(rhs, rand_terms)[[1]]
    rand_terms_clean <- gsub("^\\(|\\)$", "", rand_matches)
    all_terms <- c(fixed_terms, rand_terms_clean)
    return(all_terms)
  }
  
  #### Start of actual function ####
  
  # Get fixed effects coefficients and variance-covariance matrix
  coefs <- fixef(model, robust=T)
  coefs_df <-as.data.frame(coefs)
  V <- vcov(model, dpar = "mu")
  dfname = "data"
  
  # Get contrasts for the factor to find the missing (hidden) level
  
  test <- as.data.frame(attr(model[[dfname]][[factor_name]], "contrasts"))
  cat <- tail(rownames(test), 1) # Reference category label
  target_coef_names <- paste0(factor_name, names(test))
  
  var_idx <- which(rownames(coefs_df) %in% target_coef_names)
 
  posterior <- as_draws_matrix(model)
  hidden_draws1 <- posterior[, var_idx, drop = FALSE]
  hidden_draws1=as.data.frame(hidden_draws1)
  hidden_draws1$b_hidden = rowSums(hidden_draws1)*-1
    
    if (transform == F){
      est_hidden = median(hidden_draws1$b_hidden)
      se_hidden <- sd(hidden_draws1$b_hidden)
      ci <- quantile(hidden_draws1$b_hidden, c(.025, .975))
      lower_hidden <- unname(ci[1])
      upper_hidden <- unname(ci[2])
      
      main_effect_df <- data.frame(
        Term = paste(factor_name, cat, sep = ""),
        Estimate = round(est_hidden, digits=2),
        `Est.Error` = round(se_hidden, digits=2),
        `l-95% CI` = round(lower_hidden, digits=2),
        `u-95% CI` = round(upper_hidden, digits=2),
        check.names = FALSE
      )
    } else if (transform == T){
      est_hidden = median(exp(hidden_draws1$b_hidden))
      se_hidden <- sd(exp(hidden_draws1$b_hidden))
      ci <- quantile(exp(hidden_draws1$b_hidden), c(.025, .975))
      lower_hidden <- unname(ci[1])
      upper_hidden <- unname(ci[2])
      
      main_effect_df <- data.frame(
        Term = paste(factor_name, cat, sep = ""),
        Estimate = round(est_hidden, digits=2),
        `Est.Error` = round(se_hidden, digits=2),
        `l-95% CI` = round(lower_hidden, digits=2),
        `u-95% CI` = round(upper_hidden, digits=2),
        check.names = FALSE
      )
    }
    
    formula <- model[["formula"]]
    terms_chr <- get_terms_chr(formula[1]$formula)
  
  # Find all interactions involving factor_name (with at least two factors)
  relevant_interactions <- terms_chr[grepl(factor_name, terms_chr) & grepl(":", terms_chr)]
  
  # Separate two-way and three-way interactions
  two_way <- relevant_interactions[sapply(strsplit(relevant_interactions, ":"), length) == 2]
  three_way <- relevant_interactions[sapply(strsplit(relevant_interactions, ":"), length) == 3]
  
  interaction_two_way_dfs <- list()
  interaction_three_way_dfs <- list()
  
  # ---- Handle two-way interactions ----
  for (interaction in two_way) {
    #interaction = two_way[1]
    factors <- unlist(strsplit(interaction, ":"))
    factorA <- factors[1]
    factorB <- factors[2]
    
    devA <- (factorA == factor_name)
    devB <- (factorB == factor_name)
    
    # Get levels for factors
    levelsA <- levels(model[[dfname]][[factorA]])
    levelsB <- levels(model[[dfname]][[factorB]])
    
    nonrefA <- get_nonref(levelsA, devA)
    nonrefB <- get_nonref(levelsB, devB)
    
    refA <- get_reference_level(model, factorA, devA)
    refB <- get_reference_level(model, factorB, devB)
    
    if (devA) {
      hidden_label <- paste0(
        factorA, refA, ":",
        factorB, nonrefB
      )
    } else if (devB) {
      # In case factor_name is second factor, adjust accordingly
      hidden_label <- paste0(
        factorA, nonrefA, ":",
        factorB, refB
      )
    } else {
      # If factor_name not in first or second factor, skip
      next
    }
      rn <- rownames(coefs_df)
      
      for (k in seq_along(hidden_label)) {
      
      var_idx <- which(
        sapply(strsplit(rn, ":"), length) == 2 &   # exactly two-way
          grepl(factorA, rn, fixed = TRUE) &
          grepl(factorB, rn, fixed = TRUE) &
          grepl(nonrefB[k], rn, fixed = TRUE)
      )
      if (length(var_idx) == 0) next
      ##V_var <- V[var_idx, var_idx, drop = FALSE]
    
    
    hidden_label
      posterior <- as_draws_matrix(model)
      hidden_draws1 <- posterior[, var_idx, drop = FALSE]
      hidden_draws1=as.data.frame(hidden_draws1)
      hidden_draws1$b_hidden = rowSums(hidden_draws1)*-1
      
      if (transform == F){
        est_hidden = median(hidden_draws1$b_hidden)
        se_hidden <- sd(hidden_draws1$b_hidden)
        ci <- quantile(hidden_draws1$b_hidden, c(.025, .975))
        lower_hidden <- unname(ci[1])
        upper_hidden <- unname(ci[2])
        
        interaction_two_way_dfs[[length(interaction_two_way_dfs) + 1]] <- data.frame(
          Term = hidden_label,
          Estimate = round(est_hidden, digits=2),
          `Est.Error` = round(se_hidden, digits=2),
          `l-95% CI` = round(lower_hidden, digits=2),
          `u-95% CI` = round(upper_hidden, digits=2),
          check.names = FALSE
        )
      } else if (transform == T){
        est_hidden = median(exp(hidden_draws1$b_hidden))
        se_hidden <- sd(exp(hidden_draws1$b_hidden))
        ci <- quantile(exp(hidden_draws1$b_hidden), c(.025, .975))
        lower_hidden <- unname(ci[1])
        upper_hidden <- unname(ci[2])
        
        interaction_two_way_dfs[[length(interaction_two_way_dfs) + 1]] <- data.frame(
          Term = hidden_label,
          Estimate = round(est_hidden, digits=2),
          `Est.Error` = round(se_hidden, digits=2),
          `l-95% CI` = round(lower_hidden, digits=2),
          `u-95% CI` = round(upper_hidden, digits=2),
          check.names = FALSE
        )
      }
    } 

    }
  
  
  # ---- Handle three-way interactions ----
  for (interaction in three_way) {
    #interaction = three_way[1]
    factors <- unlist(strsplit(interaction, ":"))
    factorA <- factors[1]
    factorB <- factors[2]
    factorC <- factors[3]
    
    devA <- (factorA == factor_name)
    devB <- (factorB == factor_name)
    devC <- (factorC == factor_name)
    
    levelsA <- levels(model[[dfname]][[factorA]])
    levelsB <- levels(model[[dfname]][[factorB]])
    levelsC <- levels(model[[dfname]][[factorC]])
    
    nonrefA <- get_nonref(levelsA, devA)
    nonrefB <- get_nonref(levelsB, devB)
    nonrefC <- get_nonref(levelsC, devC)
    
    refA <- get_reference_level(model, factorA, devA)
    refB <- get_reference_level(model, factorB, devB)
    refC <- get_reference_level(model, factorC, devC)
    
    # Construct hidden label depending on which factor is factor_name
    if (devA) {
      # factor_name is first
      hidden_label <- paste0(
        factorA, refA, ":",
        factorB, nonrefB, ":",
        factorC, nonrefC
      )
    } else if (devB) {
      # factor_name is second
      hidden_label <- paste0(
        factorA, nonrefA, ":",
        factorB, refB, ":",
        factorC, nonrefC
      )
    } else if (devC) {
      # factor_name is third
      hidden_label <- paste0(
        factorA, nonrefA, ":",
        factorB, nonrefB, ":",
        factorC, refC
      )
    } else {
      next
    }
    
      rn <- rownames(coefs_df)
      
      for (k in seq_along(hidden_label)) {
      var_idx <- which(
        sapply(strsplit(rn, ":"), length) == 3 &   # exactly three-way
          grepl(factorA, rn, fixed = TRUE) &
          grepl(factorB, rn, fixed = TRUE) &
          grepl(factorC, rn, fixed = TRUE) &
          grepl(nonrefB[k], rn, fixed = TRUE)
      )
      
      if (length(var_idx) == 0) next
    
      posterior <- as_draws_matrix(model)
      hidden_draws1 <- posterior[, var_idx, drop = FALSE]
      hidden_draws1=as.data.frame(hidden_draws1)
      hidden_draws1$b_hidden = rowSums(hidden_draws1)*-1
      
      if(transform == F){
        est_hidden = median(hidden_draws1$b_hidden)
        se_hidden <- sd(hidden_draws1$b_hidden)
        ci <- quantile(hidden_draws1$b_hidden, c(.025, .975))
        lower_hidden <- unname(ci[1])
        upper_hidden <- unname(ci[2])
        
        interaction_three_way_dfs[[length(interaction_three_way_dfs) + 1]] <- data.frame(
          Term = hidden_label,
          Estimate = round(est_hidden, digits=2),
          `Est.Error` = round(se_hidden, digits=2),
          `l-95% CI` = round(lower_hidden, digits=2),
          `u-95% CI` = round(upper_hidden, digits=2),
          check.names = FALSE
        ) 
      } else if (transform == T){
        est_hidden = median(exp(hidden_draws1$b_hidden))
        se_hidden <- sd(exp(hidden_draws1$b_hidden))
        ci <- quantile(exp(hidden_draws1$b_hidden), c(.025, .975))
        lower_hidden <- unname(ci[1])
        upper_hidden <- unname(ci[2])
        
        interaction_three_way_dfs[[length(interaction_three_way_dfs) + 1]] <- data.frame(
          Term = hidden_label,
          Estimate = round(est_hidden, digits=2),
          `Est.Error` = round(se_hidden, digits=2),
          `l-95% CI` = round(lower_hidden, digits=2),
          `u-95% CI` = round(upper_hidden, digits=2),
          check.names = FALSE
        )
      }
      }
  }
  
  # Combine all results
  final_df <- rbind(
    main_effect_df,
    if (length(interaction_two_way_dfs) > 0) do.call(rbind, interaction_two_way_dfs) else NULL,
    if (length(interaction_three_way_dfs) > 0) do.call(rbind, interaction_three_way_dfs) else NULL
  )
  
  return(final_df)
}

    
Extract_BRMS <- function(brms_model, fontsize, filename, transform=TRUE){
  model_name <- deparse(substitute(brms_model))
  ## Two helper functions:
  # 1)
  get_sigma2_and_ICC_brms <- function(brms_model) {
    if (!inherits(brms_model, "brmsfit")) {
      stop("This helper expects a 'brmsfit' model.")
    }
    vars_brms <- performance::variance_decomposition(brms_model)
    sig2 <- attr(vars_brms, "var_residual")
    icc <- vars_brms$ICC_decomposed
    return(list(sigma2 = round(as.numeric(sig2), 2),
                ICC = round(as.numeric(icc), 2)))
  }
  
  # 2)
  Effects_brms_word <- function(Both, fontsize, file = filename) {
    
    data = Both
    
    make_ft_with_internal_titles <- function(data, fontsize) {
      
      header_rows <- which(data$Parameter == "Predictors")
      title_rows <- which(data$Parameter %in% c("Fixed Effects", "Random Effects"))
      Observations_row <- which(data$Parameter == "Observations")
      CrI_low_col <- which(names(data) == "CI_low")
      CrI_high_col <- which(names(data) == "CI_high")
      n_cols <- ncol(data)
      all_cols <- seq_len(n_cols)
      cri_cols <- CrI_low_col:CrI_high_col
      other_cols <- setdiff(all_cols, cri_cols)
      
      ft <- flextable(data) %>%
        delete_part(part = "header") %>% 
        fontsize(size = fontsize, part = "body") %>%
        line_spacing(space = 0.3, part = "body") %>%
        valign(valign = "center", part = "body") %>%
        padding(padding.top = fontsize*(6/8), padding.bottom = fontsize*(3/8), part = "body") %>%
        align(j = 1, align = "left", part = "body") %>%
        align(j = 2:ncol(data), align = "center", part = "body") %>%
        merge_at(i = 2, j = 1:ncol(data), part = "body") %>%
        merge_at(i= header_rows, j = CrI_low_col:CrI_high_col, part = "body") %>%
        compose(
          i = header_rows,
          j = CrI_low_col,
          part = "body",
          value = as_paragraph("95% CrI")
        ) %>%
        line_spacing(i = 2, space = 1, part = "body") %>%
        padding(i = 2, padding.top = fontsize*(1/8), padding.bottom = fontsize*(6/8), part = "body") %>% 
        autofit() #%>%
      #set_table_properties(layout = "autofit")
      #set_table_properties(layout="fixed", width=1)
      
      ft <- set_table_properties(ft, layout = "fixed", width = 1)
      
      #colwidths = ft[["body"]][["colwidths"]]
      #colwidths
      #padding = ft
      CrI_low_col_width <- ft[["body"]][["colwidths"]]["CI_low"]
      CrI_dash_col_width <- ft[["body"]][["colwidths"]]["dash"]
      CrI_high_col_width <- ft[["body"]][["colwidths"]]["CI_high"]
      
      ft <- ft %>%
        width(j = CrI_low_col, width = CrI_low_col_width*0.8) %>%
        width(j = CrI_low_col+1, width = CrI_dash_col_width*0.4) %>%
        width(j = CrI_high_col, width = CrI_high_col_width*0.8)
      
      #colwidths = ft[["body"]][["colwidths"]]
      #colwidths
      
      ft <- ft %>%
        align(j = CrI_low_col, align = "right", part = "body") %>%
        align(j = CrI_high_col, align = "left", part = "body") %>%
        align(i=header_rows, j = CrI_low_col, align = "center", part = "body")
      
      ft <- ft %>% 
        padding(j = CrI_low_col, padding.left = 5, padding.right = 1, part = "body") %>%
        padding(j = CrI_low_col+1, padding.left = 1, padding.right = 1, part = "body") %>%
        padding(j = CrI_high_col, padding.left = 1, padding.right = 5, part = "body")
      
      ## ===== INLINE BOLD FOR MODEL & FORMULA =====
      
      formula = brms_model[["formula"]][["formula"]]
      formula_chr = paste(formula[2], formula[1], formula[3])
      
      ft <- compose(
        ft,
        i = 1, j = 1,
        value = as_paragraph(
          as_chunk("Model: ", props = fp_text(bold = TRUE, font.size= fontsize)),
          as_chunk(model_name)
        ),
        part = "body"
      )
      
      ft <- compose(
        ft,
        i = 2, j = 1,
        value = as_paragraph(
          as_chunk("Formula: ", props = fp_text(bold = TRUE, font.size= fontsize)),
          as_chunk(formula_chr)
        ),
        part = "body"
      )
      ## ==========================================
      
      
      
      if(length(header_rows) > 0) {
        ft <- flextable::bold(ft, i = header_rows, bold = TRUE, part = "body")
      }
      
      #if(length(header_rows) > 0) {
      #  ft <- bold(ft, i = header_rows, bold = TRUE, part = "body")
      #}
      
      #class(ft)
      
      std_border <- fp_border(color = "black", width = 1)
      
      #class(ft)
      
      for (i in header_rows) {
        if (i > 1) {
          ft <- border(ft, border.bottom = std_border, part = "body", i = i - 1)
        }
        ft <- border(ft, border.top = std_border, part = "body", i = i)
        ft <- border(ft, border.bottom = std_border, part = "body", i = i)
      }
      
      
      ft <- border(ft, border.bottom = std_border, part = "body", i = Observations_row - 1)
      ft <- border(ft, border.top = std_border, part = "body", i = Observations_row)
      
      if(length(title_rows) > 0) {
        ft <- flextable::bold(ft, i = title_rows, part = "body") %>%
          bg(i = title_rows, bg = "#D3D3D3", part = "body")
        
        for (i in title_rows) {
          if (i > 1) {
            ft <- border(ft, border.bottom = std_border, part = "body", i = i - 1)
          }
          ft <- border(ft, border.top = std_border, part = "body", i = i)
          ft <- border(ft, border.bottom = std_border, part = "body", i = i)
        }
      }
      
      ft <- border(ft, border.bottom = fp_border(style = "none"), i = nrow(data))
      ft <- border(ft, border.top = fp_border(style = "double", width = 1), i = 1)
      
      return(ft)
    }
    doc <- read_docx() %>%
      body_add_flextable(make_ft_with_internal_titles(Both, fontsize))
    
    print(doc, target = file)
    file
  }
  
  message(crayon::green("Analysing posteriors...(this may take a few minutes)"))
  ## ===== Start extraction =====
  posteriors <- describe_posterior(
    brms_model,
    effects = "all",
    component = "all",
    test = "p_direction",
    centrality = "all"
  )
  posteriors$OR = exp(posteriors$Median)
  posteriors$CI_low_OR <- exp(posteriors$CI_low)
  posteriors$CI_high_OR <- exp(posteriors$CI_high)
  
  posteriors <<- posteriors
  #CI <- quantile(draws$b_Time_PointTime2, probs = c(0.025, 0.975))
  #OR_CI= exp(CI)
  #draws$b_Time_PointTime2_OR = exp(draws$b_Time_PointTime2)
  #sd(draws$b_Time_PointTime2_OR)
  #OR_CI2 <- quantile(draws$b_Time_PointTime2_OR, probs = c(0.025, 0.975))
  
  
  library(tidyr)
  draws <- as_draws_df(brms_model)
  #mean(draws$b_Time_PointTime2)
  
  #sd(draws$b_Time_PointTime2)
  
  draws2 <- draws %>%
    dplyr::select(starts_with("b_")) %>%
    dplyr::summarise(across(everything(), sd)) %>%
    pivot_longer(everything(),
                 names_to = "Parameter",
                 values_to = "Est_Error")
  draws2$Est_Error=round(draws2$Est_Error, digits = 2)
  
  draws3 <- draws %>%
    dplyr::select(starts_with("b_")) %>%
    dplyr::mutate(across(everything(), exp)) %>%
    dplyr::summarise(across(everything(), sd)) %>%
    pivot_longer(
      everything(),
      names_to = "Parameter",
      values_to = "Est_Error_OR"
    )
  
  draws3$Est_Error_OR <- round(draws3$Est_Error_OR, 2)
  
  draws = merge(draws2, draws3, by = "Parameter")
  
  message(crayon::green("Extracting fixed effects"))
  ## Fixed Effects
  FEs <- subset(posteriors, Effects == "fixed" & substr(Parameter, 1, 1) == "b")
  FEs[, 3:14] <- round(FEs[, 3:14], 2)
  FEs <- subset(FEs, select = -c(Effects, Mean, MAP, CI, Rhat, ESS))
  names(FEs)[names(FEs) == "Median"] <- "Estimate"
  rownames(FEs)<-NULL
  FEs$ID=rownames(FEs)
  FEs = merge(FEs, draws, by = "Parameter")
  FEs$ID=as.numeric(as.character(FEs$ID))
  FEs <- FEs %>%
    arrange(ID)
  FEs=subset(FEs, select = -(ID))
  FEs$Parameter <- substring(FEs$Parameter, 3)
  FEs$Meaningful <- ""
  FEs$Meaningful[(FEs$CI_low < 0 & FEs$CI_high < 0) |
                   (FEs$CI_low > 0 & FEs$CI_high > 0)] <- "✅"
  
  FEs[] <- lapply(FEs, format, nsmall = 2)
  #transform = TRUE
  if (transform == FALSE){
    FEs <- FEs %>%
      relocate(Parameter, Estimate, Est_Error)
    FEs = subset(FEs, select = -c(OR, CI_low_OR, CI_high_OR, Est_Error_OR))
  } else if (transform == TRUE){
    FEs <- FEs %>%
      relocate(Parameter, OR, Est_Error_OR, CI_low_OR, CI_high_OR)
    FEs = subset(FEs, select = -c(Estimate, Est_Error, CI_low, CI_high))
    names(FEs)[names(FEs)=="OR"]="Estimate"
    names(FEs)[names(FEs)=="Est_Error_OR"]="Est_Error"
    names(FEs)[names(FEs)=="CI_low_OR"]="CI_low"
    names(FEs)[names(FEs)=="CI_high_OR"]="CI_high"
  }
  
  #FEs$CrI = paste(FEs$CI_low, "–", FEs$CI_high, sep="")
  #FEs$CrI = gsub(" ", "", FEs$CrI)
  #FEs$CrI = gsub("–", " – ", FEs$CrI)
  #FEs = subset(FEs, select = -c(CI_low, CI_high))
  #FEs <- FEs %>%
  #  relocate(Parameter, Estimate, Est_Error, CrI)
  FEs$dash = "–"
  FEs <- FEs %>%
    relocate(Parameter, Estimate, Est_Error, CI_low, dash, CI_high)
  
  sapply(FEs, class)
  
  pad_left_to_decimal <- function(vec, pad_char = "  ") {
    vec <- trimws(vec, which = "left")
    before_decimal <- sapply(strsplit(vec, "\\."), `[`, 1)
    lengths_before_decimal <- nchar(before_decimal)
    max_length <- max(lengths_before_decimal, na.rm = TRUE)
    
    mapply(function(x, len) {
      spaces_needed <- max_length - len
      paste0(strrep(pad_char, spaces_needed), x)
    }, vec, lengths_before_decimal, USE.NAMES = FALSE)
  }
  
  # Then apply to your columns:
  FEs$Estimate <- pad_left_to_decimal(FEs$Estimate)
  FEs$Est_Error <- pad_left_to_decimal(FEs$Est_Error)
  FEs$CI_low <- pad_left_to_decimal(FEs$CI_low)
  FEs$CI_high <- pad_left_to_decimal(FEs$CI_high)  # double spaces as in original
  FEs$pd <- pad_left_to_decimal(FEs$pd)
  
  ## Headers (Model & Formula rows left empty for compose())
  header0FE <- data.frame(Parameter = "", Estimate = "", Est_Error="", CI_low = "", dash="", CI_high="", pd="", Meaningful="")
  header0.5FE <- data.frame(Parameter = "", Estimate = "", Est_Error="", CI_low = "", dash="", CI_high="", pd="", Meaningful="")
  header1FE <- data.frame(Parameter = "Fixed Effects", Estimate = "", Est_Error="", CI_low = "", dash="", CI_high="", pd="", Meaningful="")
  if (transform == FALSE){
    header2FE <- data.frame(Parameter = "Predictors", Estimate = "Log Odds", Est_Error ="Est. Error",
                            CI_low = "CrI low", dash="–", CI_high="CrI high", pd="pd", Meaningful="Meaningful")
  } else if (transform == TRUE) {
    header2FE <- data.frame(Parameter = "Predictors", Estimate = "Odds Ratio", Est_Error ="Est. Error",
                            CI_low = "CrI low", dash="–", CI_high="CrI high", pd="pd", Meaningful="Meaningful")
  }
  
  FEs <- rbind(header0FE, header0.5FE, header1FE, header2FE, FEs)
  
  
  message(crayon::green("Extracting random effects"))
  ## Random Effects
  REs = subset(posteriors, Effects == "random" & substr(posteriors$Parameter, 1, 2)=="sd")
  REs[, 3:11] = round(REs[, 3:11], digits=2)
  REs = subset(REs, select = c(Parameter, Median))
  names(REs)[names(REs)=="Median"]="SD"
  sigma2_and_ICC = get_sigma2_and_ICC_brms(brms_model)
  REs$Parameter = substring(REs$Parameter, 4)
  REs <- REs %>%
    dplyr::arrange(!grepl("_Intercept", Parameter))
  REs$Parameter[grepl("_Intercept", REs$Parameter)==T]=paste("τ00", REs$Parameter[grepl("_Intercept", REs$Parameter)==T])
  REs$Parameter[grepl("_Intercept", REs$Parameter)==F]=paste("τ11", REs$Parameter[grepl("_Intercept", REs$Parameter)==F])
  REs$Parameter = gsub("__Intercept", "", REs$Parameter)
  REs$Parameter = gsub("__", ".", REs$Parameter)
  sigma2 = sigma2_and_ICC$sigma2
  ICC = sigma2_and_ICC$ICC
  sigma2row = data.frame(Parameter = "σ2", SD = sigma2)
  ICCrow = data.frame(Parameter = "ICC", SD = ICC)
  REs = rbind(sigma2row, REs)
  REs = rbind(REs, ICCrow)
  REs$SD=format(REs$SD, nsmall =2)
  REs$SD <- pad_left_to_decimal(REs$SD)
  
  rand_eff <- insight::get_data(brms_model, verbose = FALSE)[, insight::find_random(brms_model, split_nested = TRUE, flatten = TRUE), drop = FALSE]
  n_re_grps <- sapply(rand_eff, function(.i) length(unique(.i, na.rm = TRUE)))
  names(n_re_grps) <- sprintf("N %s", names(n_re_grps))
  reNsDF = as.data.frame(n_re_grps)
  names(reNsDF)[names(reNsDF)=="n_re_grps"]="SD"
  reNsDF$Parameter= rownames(reNsDF)
  REs = rbind(REs, reNsDF)
  names(REs)[names(REs)=="SD"]="Estimate"
  for (i in 1: ncol(REs)){
    REs[, i]=format(REs[, i], nsmall=2)
  }

  header1 = data.frame(Parameter = "Random Effects", Estimate = "")
  #header2 = data.frame(Parameter = "Parameter", Estimate = "SD")
  #headers = rbind(header1, header2)
  REs= rbind(header1, REs)
  n_obs <- nobs(brms_model)
  row <- data.frame(Parameter = "Observations", Estimate = n_obs)
  REs = rbind(REs, row)
  rsqdummy <- tryCatch(suppressWarnings(performance::r2(brms_model)),
                       error = function(x) NULL)
  if (!is.null(rsqdummy)) {
    row1 <- data.frame(Parameter = "Marginal R2", 
                       Estimate = round(rsqdummy$R2_Bayes_marginal, digits =3))
    row2 <- data.frame(Parameter = "Conditional R2", Estimate = round(rsqdummy$R2_Bayes, digits = 3))
    r2sdf = rbind(row1, row2)
  }
  REs = rbind(REs, r2sdf)
  REs$Est_Error = ""
  REs$CI_low=""
  REs$dash=""
  REs$CI_high=""
  REs$pd=""
  REs$Meaningful=""
  
  Both <- rbind(FEs, REs)
  rownames(Both)<-NULL
  
  message(crayon::green("Saving output to word..."))
  Effects_brms_word(Both, fontsize, filename)
}
                               
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




RanSlope_Tester_Auto <- function(
  DF, dv, var, RanIntercepts,
  include_lower_order = TRUE,
  verbose = TRUE,
  return_table = FALSE,
  w_small = 0.2, w_unbalanced = 0.2, w_variation = 0.6,
  # --- Customizable thresholds ---
  min_cluster_size = NULL,       # minimum obs per cluster (default: median - SD)
  min_continuous_sd = 0.1,       # minimum meaningful variation as fraction of overall SD
  min_cat_prop = 0.05            # minimum proportion to avoid "highly unbalanced"
) {
  # --- Load required packages ---
  required_pkgs <- c("dplyr", "glue", "crayon", "scales", "rlang")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(glue::glue("Missing required packages: {paste(missing_pkgs, collapse = ', ')}"))
  }
  
  msg <- function(text, color = "white") {
    if (verbose) cat(do.call(crayon::style, list(text, color)), "\n")
  }
  
  # --- Check grouping variables ---
  for (g in RanIntercepts) {
    if (anyNA(DF[[g]])) warning(glue::glue("Grouping variable '{g}' has missing values."))
    if (length(unique(DF[[g]])) < 2) stop(glue::glue("Grouping variable '{g}' has < 2 clusters."))
  }
  
  # --- Core diagnostic function ---
  RanSlope_Tester12 <- function(DF, dv, var, RanIntercepts) {
    temp_var_name <- ".temp_var_for_ranslope_check"
    
    if (grepl("\\*", var)) {
      var_terms <- all.vars(rlang::parse_expr(var))
      DF[[temp_var_name]] <- interaction(DF[var_terms], drop = TRUE)
    } else {
      var_expr <- rlang::parse_expr(var)
      DF[[temp_var_name]] <- with(DF, eval(var_expr))
    }
    
    var_is_continuous <- is.numeric(DF[[temp_var_name]])
    
    check_dv_variance <- function(dv_vector) {
      if (is.numeric(dv_vector)) return(var(dv_vector, na.rm = TRUE) > 0)
      present_levels <- unique(dv_vector[!is.na(dv_vector)])
      if (length(present_levels) <= 1) return(FALSE)
      tab <- table(dv_vector)
      return(all(tab / sum(tab) > 0))
    }
    
    results <- list()
    
    for (RanIntercept in RanIntercepts) {
      msg(glue::glue("\nChecking random slope for {var} within {RanIntercept}"), "blue")
      
      cluster_sizes <- DF %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop")
      
      # --- Determine minimum cluster size ---
      if (is.null(min_cluster_size)) {
        min_cluster_n <- max(2, floor(median(cluster_sizes$n) - sd(cluster_sizes$n)))
      } else {
        min_cluster_n <- min_cluster_size
      }
      
      small_clusters <- sum(cluster_sizes$n < min_cluster_n, na.rm = TRUE)
      total_clusters <- nrow(cluster_sizes)
      prop_small_clusters <- small_clusters / total_clusters
      if (small_clusters > 0)
        msg(glue::glue("⚠️ {small_clusters}/{total_clusters} groups < {min_cluster_n} obs."), "yellow")
      
      prop_passing <- NA
      prop_unbalanced <- 0
      
      if (var_is_continuous) {
        overall_sd <- stats::sd(DF[[temp_var_name]], na.rm = TRUE)
        cluster_sds <- DF %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
          dplyr::summarise(sd_val = stats::sd(.data[[temp_var_name]], na.rm = TRUE), .groups = "drop")
        
        threshold_sd <- max(min_continuous_sd * overall_sd, 1e-8)
        cluster_sds <- cluster_sds %>%
          dplyr::mutate(meaningful_variation = !is.na(sd_val) & sd_val >= threshold_sd)
        
        n_with_variation <- sum(cluster_sds$meaningful_variation, na.rm = TRUE)
        prop_passing <- n_with_variation / total_clusters
        prop_failing <- 1-prop_passing
        #prop_failing <- format(prop_failing, digits = 3)
        
        msg(glue::glue("{n_with_variation}/{total_clusters} ({scales::percent(prop_passing)}) groups show meaningful within-cluster variation: Props_Clusters_Failing = {format(round(prop_failing, 3), nsmall = 3)}"), 
            ifelse(prop_passing == 0, "red", ifelse(prop_passing < 0.5, "yellow", "green")))
        
      } else {
        counts <- DF %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(c(RanIntercept, temp_var_name)))) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop")
        
        replicated <- counts %>% dplyr::filter(n >= 2)
        replicated_df <- DF %>% dplyr::semi_join(replicated, by = c(RanIntercept, temp_var_name))
        
        dv_variance_check <- replicated_df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(c(RanIntercept, temp_var_name)))) %>%
          dplyr::summarise(has_variance = check_dv_variance(.data[[dv]]), .groups = "drop")
        
        valid_combos <- dv_variance_check %>% dplyr::filter(has_variance)
        levels_per_group <- valid_combos %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
          dplyr::summarise(levels_with_variance = dplyr::n_distinct(.data[[temp_var_name]]), .groups = "drop")
        
        n_multilevel <- sum(levels_per_group$levels_with_variance >= 2, na.rm = TRUE)
        prop_passing <- n_multilevel / total_clusters
        prop_failing <- 1-prop_passing
        
        # --- Unbalanced groups based on min_cat_prop ---
        level_props <- counts %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(RanIntercept))) %>%
          dplyr::mutate(prop = n / sum(n)) %>%
          dplyr::ungroup()
        
        unbalanced_groups <- level_props %>%
          dplyr::filter(prop < min_cat_prop) %>%
          dplyr::distinct(dplyr::across(dplyr::all_of(RanIntercept)))
        
        prop_unbalanced <- nrow(unbalanced_groups) / total_clusters
        if (nrow(unbalanced_groups) > 0)
          msg(glue::glue("⚠️ {nrow(unbalanced_groups)} groups are highly unbalanced."), "yellow")
        
        msg(glue::glue("{n_multilevel}/{total_clusters} ({scales::percent(prop_passing)}) groups show ≥2 levels with DV variance: Prop_Clusters_Failing = {format(round(prop_failing, 3), nsmall = 3)}"),
            ifelse(prop_passing == 0, "red", ifelse(prop_passing < 0.5, "yellow", "green")))
      }
      
      # --- Risk score calculation (normalized) ---
      if (prop_failing == 1) {
        risk_score <- 1
      } else {
        risk_score <- w_small * prop_small_clusters +
              w_unbalanced * prop_unbalanced +
              w_variation * prop_failing
        #risk_score <- min(risk_score, 1)  # cap at 1
      }
    
      results[[RanIntercept]] <- list(
        Grouping_Factor = RanIntercept,
        Prop_Small_Groups = prop_small_clusters,
        Prop_Unbalanced = prop_unbalanced,
        Prop_Clusters_Failing = prop_failing,
        Risk_Score = risk_score
      )
    }
    
    dplyr::bind_rows(lapply(results, as.data.frame))
  }
  
  # --- Expand interactions if needed ---
  expand_interactions <- function(var) {
    parts <- trimws(unlist(strsplit(var, "\\*")))
    n <- length(parts)
    if (n == 1 || !include_lower_order) return(var)
    combos <- unlist(lapply(1:(n - 1), function(k) apply(combn(parts, k), 2, paste, collapse = "*")))
    return(c(combos, var))
  }
  
  effects <- if (grepl("\\*", var)) expand_interactions(var) else var
  
  # --- Run diagnostics ---
  all_results <- lapply(effects, function(eff) {
    res <- RanSlope_Tester12(DF, dv, eff, RanIntercepts)
    res$Effect <- eff
    res$Effect_Type <- if (grepl("\\*", eff)) "Interaction" else "Main Effect"
    res
  })
  
  combined <- dplyr::bind_rows(all_results) %>%
    dplyr::select(Effect, Effect_Type, dplyr::everything())
  
  # --- Downgrade interactions if lower-order terms fail ---
  if (any(combined$Effect_Type == "Interaction")) {
    for (eff in combined$Effect[combined$Effect_Type == "Interaction"]) {
      parts <- unlist(strsplit(eff, "\\*"))
      for (grp in unique(combined$Grouping_Factor)) {
        lower_scores <- combined %>%
          dplyr::filter(Effect %in% parts, Grouping_Factor == grp) %>%
          dplyr::pull(Risk_Score)
        if (any(lower_scores == 1)) {
          combined$Risk_Score[combined$Effect == eff & combined$Grouping_Factor == grp] <- 1
        }
      }
    }
  }
  
  # --- Add Overall Recommendation ---
  combined <- combined %>%
    dplyr::mutate(
      Overall_Recommendation = dplyr::case_when(
        Risk_Score == 1 ~ "Impossible",
        Risk_Score > 0.5 ~ "High Risk",
        Risk_Score > 0.25 ~ "Medium Risk",
        TRUE ~ "Low Risk"
      )
    )
  
  # --- Print colored table ---
  if (verbose) {
  # Convert numeric columns to strings to measure width
  combined_str <- combined %>%
    dplyr::mutate(
      Prop_Small_Groups_str = sprintf("%.3f", Prop_Small_Groups),
      Prop_Unbalanced_str = sprintf("%.3f", Prop_Unbalanced),
      Prop_Clusters_Failing_str = sprintf("%.3f", Prop_Clusters_Failing),
      Risk_Score_str = sprintf("%.3f", Risk_Score)
    )

  # Determine max widths (header vs content)
  effect_width <- max(nchar(as.character(combined$Effect)), nchar("Effect"))
  type_width <- max(nchar(as.character(combined$Effect_Type)), nchar("Effect_Type"))
  group_width <- max(nchar(as.character(combined$Grouping_Factor)), nchar("Grouping_Factor"))
  small_width <- max(nchar(combined_str$Prop_Small_Groups_str), nchar("Prop_Small_Groups"))
  unbal_width <- max(nchar(combined_str$Prop_Unbalanced_str), nchar("Prop_Unbalanced"))
  passing_width <- max(nchar(combined_str$Prop_Clusters_Failing_str), nchar("Prop_Clusters_Failing"))
  risk_width <- max(nchar(combined_str$Risk_Score_str), nchar("Risk_Score"))
  rec_width <- max(nchar(as.character(combined$Overall_Recommendation)), nchar("Overall_Recommendation"))
  
  # Header
  cat(sprintf(
    paste0(
      "%-", effect_width, "s  %-", type_width, "s  %-", group_width, "s  ",
      "%-", small_width, "s  %-", unbal_width, "s  %-", passing_width, "s  ",
      "%-", risk_width, "s  %-", rec_width, "s\n"
    ),
    "Effect", "Effect_Type", "Grouping_Factor",
    "Prop_Small_Groups", "Prop_Unbalanced",
    "Prop_Clusters_Failing", "Risk_Score", "Overall_Recommendation"
  ))
  cat(strrep("-", effect_width + type_width + group_width + small_width + unbal_width +
                   passing_width + risk_width + rec_width + 14), "\n")
  
  # Rows
  for (i in 1:nrow(combined)) {
    row <- combined[i, ]
    rec_colored <- switch(
      as.character(row$Overall_Recommendation),
      "Impossible" = crayon::red$bold((row$Overall_Recommendation)),
      "High Risk" = crayon::red(row$Overall_Recommendation),
      "Medium Risk" = crayon::yellow(row$Overall_Recommendation),
      "Low Risk" = crayon::green(row$Overall_Recommendation)
    )
    cat(sprintf(
      paste0(
        "%-", effect_width, "s  %-", type_width, "s  %-", group_width, "s  ",
        "%-", small_width, ".3f  %-", unbal_width, ".3f  %-", passing_width, ".3f  ",
        "%-", risk_width, ".3f  %-", rec_width, "s\n"
      ),
      row$Effect, row$Effect_Type, row$Grouping_Factor,
      row$Prop_Small_Groups, row$Prop_Unbalanced,
      row$Prop_Clusters_Failing, row$Risk_Score,
      rec_colored
    ))
  }
}

    
  if (return_table) return(combined) else invisible(combined)
}




































