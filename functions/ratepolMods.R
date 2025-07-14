

run_iterationMod = function (data_source_run, bin_selection = "first", standardise = FALSE, 
          n_individuals = 150, tranform_to_proportions = TRUE, dissimilarity_coefficient = "euc", 
          time_standardisation = 500, verbose = FALSE, legacySamples=5, lofK = 5) 
{
  
  
  data_subset <- subset_samples(data_source_subset = data_source_run$data, 
                                data_source_bins = data_source_run$bins, bin_selection = bin_selection)
  
  data_subset <- reduce_data_simple(data_source_reduce = data_subset)
  
  if (isTRUE(standardise)) {
    com_data_sums <- rowSums(subset_community(data_source = data_subset), 
                             na.rm = TRUE)
    n_individuals <- min(c(com_data_sums, n_individuals))
    data_subset <- data_subset[com_data_sums >= n_individuals, 
    ]
    data_subset <- reduce_data_simple(data_source_reduce = data_subset)
    data_sd <- standardise_community_data(data_source_standard = data_subset, 
                                          n_individuals = n_individuals)
    if (isTRUE(verbose)) {
      assertthat::assert_that(all(n_individuals == rowSums(subset_community(data_sd), 
                                                           na.rm = TRUE)), msg = paste("Data standardisation was unsuccesfull,", 
                                                                                       "try 'standardise' = FALSE"))
    }
  }
  else {
    data_sd <- data_subset
  }
  
  data_sd <- reduce_data_simple(data_source_reduce = data_sd)
  
    if (isTRUE(tranform_to_proportions)) {
    data_sd_prop <- transform_into_proportions(data_source_trans = data_sd, 
                                               sel_method = "proportions")
  }
  else {
    data_sd_prop <- data_sd
  }
  
  ########################
  # STAPLES CODE INJECTION
  ########################
  
  # calculate dissimilarity matrix. Make sure matrix is ordered old -> young
  # so that legacy and RoC/LoF are calculated in the correct time order
  rownames(data_sd_prop) = data_sd_prop$res_age
  data_sd_prop_ord = data_sd_prop[order(data_sd_prop$res_age, decreasing=TRUE),-(1:3)]
  
  dissMat = as.matrix(vegdist(data_sd_prop_ord, method=dissimilarity_coefficient))
  
  # post-sample dissimilarities to time N-1 and N-2)
  dN = t(sapply(1:nrow(dissMat), function(n){

  if(n == nrow(dissMat)){return(matrix(NA, nrow=legacySamples, ncol=2))}
    
  if((n+legacySamples) < nrow(dissMat)){
  n1 = dissMat[(n+1):(n+legacySamples),n]
  n1gap = as.numeric(names(n1)) - as.numeric(colnames(dissMat)[n])
  } else{ 
  n1 = dissMat[(n+1):nrow(dissMat),n]
  names(n1) = rownames(dissMat)[(n+1):nrow(dissMat)]
  n1gap = as.numeric(names(n1)) - as.numeric(colnames(dissMat)[n])
  n1 = c(n1, rep(NA, legacySamples-length(n1)))
  n1gap = c(n1gap, rep(NA, legacySamples-length(n1gap)))
  }
    
    return(cbind(matrix(n1, ncol=legacySamples), matrix(n1gap, ncol=legacySamples)))
  }))
  colnames(dN) = c(paste0("dn1_", 1:legacySamples),
                   paste0("n1gap_", 1:legacySamples))
  rownames(dN) = rownames(dissMat)
  
  # this does rate of change equivalent to RRatepol:::estimate_dissimilarity_coefficient
  dc_res = c(NA, diag(dissMat[-1,-ncol(dissMat)]))
  
  # use LOF as a measurement of compositional anomaly/novelty
  lof = do.call("rbind", lapply((lofK+1):nrow(dissMat), function(n){
    ecoLOF(dissMat[(1:n),], k=lofK, method=dissimilarity_coefficient)[n,]
  }))
  lof = rbind(matrix(NA, nrow=lofK, ncol=3, dimnames=list(rownames(dissMat)[1:lofK], colnames(lof))), 
              lof)
  
  # diversity estimates
  div = data.frame(H0 = hill_div(t(data_sd_prop_ord), qvalue=0),
                   H1 = hill_div(t(data_sd_prop_ord), qvalue=1),
                   H2 = hill_div(t(data_sd_prop_ord), qvalue=2))
  
  # delta diversity from N-1
  divDelta = rbind(NA,
                   apply(div, 2, diff))
  colnames(divDelta) = paste0("d", colnames(div))
  
  returnData = cbind(dN, roc = dc_res, lof, div, divDelta)
  
  # reorder returnData to young -> old to match RRatePol
  returnData = returnData[nrow(returnData):1,]

  returnData$age = data_sd_prop$res_age
  
  roc_res <- data_sd_prop[seq_along(dc_res), ] %>% dplyr::mutate(dc = dc_res, 
                                                                 age_diff_st = .data$age_diff/time_standardisation, 
                                                                 roc = .data$dc/.data$age_diff_st) %>% 
    dplyr::select("label", "res_age")
  
  roc_res = cbind(roc_res, returnData)
  
  if (nrow(roc_res) < 1) {
    stop("Estimation not succesfull")
  }
  
  return(roc_res)
}

subset_samples = function (data_source_subset, data_source_bins, bin_selection = "first") 
{
  if (is.character(data_source_bins$start)) {
    res <- data_source_bins %>% dplyr::select("label", "age_diff", 
                                              "res_age", "start") %>% dplyr::inner_join(data_source_subset %>% 
                                                                                          tibble::rownames_to_column("start"), by = "start") %>% 
      dplyr::mutate(age_diff = c(diff(.data$age), Inf), 
                    age_diff = ifelse(.data$age_diff == 0, 0.1, .data$age_diff), 
                    res_age = .data$age) %>% dplyr::select(-c("start", 
                                                              "age"))
    return(res)
  }
  res_com <- as.data.frame(matrix(nrow = nrow(data_source_bins), 
                                  ncol = ncol(data_source_subset) - 1, dimnames = list(NULL, 
                                                                                       names(data_source_subset)[2:ncol(data_source_subset)])))
  for (i in 1:nrow(data_source_bins)) {
    subset_w <- data_source_subset[data_source_subset$age >= 
                                     data_source_bins$start[i] & data_source_subset$age < 
                                     data_source_bins$end[i], ]
    if (nrow(subset_w) > 0) {
      if (bin_selection == "random") {
        random_row <- sample(1:nrow(subset_w), 1)
        res_com[i, ] <- subset_w[random_row, -1]
      }
      if (bin_selection == "first") {
        res_com[i, ] <- subset_w[1, -1]
      }
    }
  }
  res <- dplyr::bind_cols(data_source_bins %>% dplyr::select("label", 
                                                             "age_diff", "res_age"), res_com)
  return(res)
}

reduce_data_simple = function (data_source_reduce, ommit_vars = c("label", "res_age", 
                                             "age_diff"), check_taxa = TRUE, check_levels = TRUE) 
{
  data_com <- subset_community(data_source_reduce, ommit_vars = ommit_vars)
  if (isTRUE(check_taxa)) {
    valid_taxa <- (colSums(data_com, na.rm = TRUE) > 0)
    data_source_reduce <- data_source_reduce %>% dplyr::select(dplyr::any_of(c(ommit_vars, 
                                                                               names(valid_taxa[valid_taxa]))))
  }
  if (isTRUE(check_levels)) {
    valid_levels <- (rowSums(data_com, na.rm = TRUE) > 0)
    data_source_reduce <- data_source_reduce[valid_levels, 
    ]
  }
  return(data_source_reduce)
}

subset_community = function (data_source, ommit_vars = c("label", "res_age", "age_diff", 
                                                         "age")) 
{
  data_source %>% dplyr::select(!dplyr::any_of(c(ommit_vars))) %>% 
    return()
}

standardise_community_data = function (data_source_standard, n_individuals = 150) 
{
  data_community <- subset_community(data_source_standard) %>% 
    round()
  n_taxa <- ncol(data_community)
  for (i in 1:nrow(data_community)) {
    select_row <- data_community[i, ]
    vec1 <- vector(, length = 0)
    for (j in 1:n_taxa) {
      v1 <- rep(names(select_row)[j], select_row[j])
      vec1 <- c(vec1, v1)
    }
    rsample <- sample(vec1, size = n_individuals, replace = FALSE)
    sel_names <- table(rsample)
    data_community[i, ] <- rep(0, n_taxa)
    data_community[i, names(sel_names)] <- as.numeric(sel_names)
  }
  data_source_standard[, names(data_community)] <- data_community
  return(data_source_standard)
}

transform_into_proportions = function (data_source_trans, sel_method = c("proportions", "percentages"), 
          verbose = FALSE) 
{
  RUtilpol::check_class("data_source_trans", "data.frame")
  RUtilpol::check_class("sel_method", "character")
  RUtilpol::check_vector_values("sel_method", c("percentages", 
                                                "proportions"))
  sel_method <- match.arg(sel_method)
  RUtilpol::check_class("verbose", "logical")
  if (isTRUE(verbose)) {
    RUtilpol::output_comment("Community data values are being converted to proportions")
  }
  data_com <- subset_community(data_source_trans)
  data_rowsums <- rowSums(data_com, na.rm = TRUE)
  data_com <- data_com/data_rowsums * switch(sel_method, percentages = 100, 
                                             proportions = 1)
  data_source_trans[, names(data_com)] <- data_com
  return(data_source_trans)
}