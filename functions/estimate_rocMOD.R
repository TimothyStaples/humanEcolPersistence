# data_source_community = commData
# data_source_age = ageData
# age_uncertainty = NULL
# smooth_method="shep"
# working_units="MW"
# rand = 29
# use_parallel=FALSE
# dissimilarity_coefficient = "bray"
# bin_size = 500
# time_standardisation = NULL
# bin_selection = "first"
# standardise = TRUE
# n_individuals = 150
# tranform_to_proportions = TRUE
# interest_threshold = NULL
# verbose = TRUE
# number_of_shifts = 1
# smooth_n_points = 9

estimate_rocMOD = function (data_source_community, data_source_age, age_uncertainty = NULL, 
          smooth_method = c("none", "m.avg", "grim", "age.w", "shep"), 
          smooth_n_points = 5, smooth_age_range = 500, smooth_n_max = 9, 
          working_units = c("levels", "bins", "MW"), bin_size = 500, 
          number_of_shifts = 5, bin_selection = c("random", "first"), 
          standardise = FALSE, n_individuals = 150, dissimilarity_coefficient = c("euc", 
                                                                                  "euc.sd", "chord", "chisq", "gower", "bray"), 
          tranform_to_proportions = TRUE, 
          rand = NULL, use_parallel = FALSE, interest_threshold = NULL, 
          time_standardisation = NULL, verbose = FALSE, ...) 
{
  
  # CHECK ARGUMENT CLASS
  assertthat::assert_that(!missing(data_source_community),
                         msg = "Object 'data_source_community' must be included as a 'data.frame'")
  assertthat::assert_that(!missing(data_source_age), msg = "Object 'data_source_age' must be included as a 'data.frame'")
  RUtilpol::check_class("data_source_community", "data.frame")
  RUtilpol::check_class("data_source_age", "data.frame")
  RUtilpol::check_class("age_uncertainty", c("NULL", "matrix"))
  RUtilpol::check_class("working_units", "character")
  RUtilpol::check_vector_values("working_units", c("levels",
                                                  "bins", "MW"))
  working_units <- match.arg(working_units)
  
  # SET TIME BINS
  if (is.null(time_standardisation)) {
    time_standardisation <- bin_size
  }
  
  RUtilpol::check_class("time_standardisation", "numeric")
  RUtilpol::check_if_integer("time_standardisation")

  if (working_units != "levels") {
    RUtilpol::check_class("bin_size", "numeric")
    RUtilpol::check_if_integer("bin_size")
    RUtilpol::check_class("bin_selection", "character")
    RUtilpol::check_vector_values("bin_selection", c("first",
                                                     "random"))
    bin_selection <- match.arg(bin_selection)
    if (working_units == "MW") {
      RUtilpol::check_class("number_of_shifts", "numeric")
      RUtilpol::check_if_integer("number_of_shifts")
    }
  }
  
  RUtilpol::check_class("standardise", "logical")
  
  if (isTRUE(standardise)) {
    RUtilpol::check_class("n_individuals", "numeric")
    RUtilpol::check_if_integer("n_individuals")
  }
  RUtilpol::check_class("tranform_to_proportions", "logical")
  RUtilpol::check_class("interest_threshold", c("NULL", "numeric"))
  RUtilpol::check_class("smooth_method", "character")
  RUtilpol::check_vector_values("smooth_method", c("none",
                                                   "m.avg", "grim", "age.w", "shep"))
  
  smooth_method <- match.arg(smooth_method)
  if (!smooth_method %in% c("none", "shep")) {
    assertthat::assert_that(smooth_n_points%%2 != 0, msg = "'smooth_n_points' must be an odd number")
    if (smooth_method != "m.avg") {
      RUtilpol::check_class("smooth_age_range", "numeric")
      if (smooth_method == "grim") {
        assertthat::assert_that(smooth_n_max%%2 != 0, 
                                msg = "'smooth_n_max' must be an odd number")
        assertthat::assert_that(smooth_n_points < smooth_n_max, 
                                msg = "'smooth_n_max' must be bigger than 'smooth_n_points")
      }
    }
  }
  
  RUtilpol::check_class("dissimilarity_coefficient", "character")
  RUtilpol::check_vector_values("dissimilarity_coefficient", 
                                c("euc", "euc.sd", "chord", "chisq", "gower", "bray"))
  dissimilarity_coefficient <- match.arg(dissimilarity_coefficient)
  
  RUtilpol::check_class("rand", c("NULL", "numeric"))
  if (isFALSE(is.null(rand))) {
    RUtilpol::check_if_integer("rand")
  }
  
  RUtilpol::check_class("use_parallel", c("logical", "numeric"))
  if (is.numeric(use_parallel)) {
    RUtilpol::check_if_integer("use_parallel")
  }
  RUtilpol::check_class("verbose", "logical")
  
  
  # START OF PROCESS ####
  
  # COMMENTS
  start_time <- Sys.time()
  RUtilpol::output_heading(paste("RRatepol started", start_time), 
                           size = "h1")
  if (isFALSE(is.null(age_uncertainty))) {
    RUtilpol::output_comment("'age_uncertainty' will be used for in the RoC estimation")
    if (rand < 100) {
      RUtilpol::output_warning(paste("'age_uncertainty' was selected to be used with low number", 
                                     "of replication. Recommend to increase 'rand'"))
    }
  }
  
  
  switch(working_units, levels = {
    RUtilpol::output_comment("RoC will be estimated between individual subsequent levels")
  }, bins = {
    RUtilpol::output_comment(paste("RoC will be estimated using selective binning with", 
                                   bin_size, "yr time bin"))
  }, MW = {
    RUtilpol::output_comment(paste("RoC will be estimated using 'binning with the mowing window' of", 
                                   bin_size, "yr time bin over", number_of_shifts, "number of window shifts"))
  })
  
  
  if (working_units != "levels") {
    if (bin_selection == "random") {
      RUtilpol::output_comment("Sample will randomly selected for each bin")
      if (rand < 100) {
        RUtilpol::output_warning(paste("'bin_selection' was selected as 'random' with low number", 
                                       "of replication. Recommend to increase 'rand'"))
      }
    }
    else {
      RUtilpol::output_comment("First sample of each time bin will selected")
    }
  }
  
  RUtilpol::output_comment(paste("'time_standardisation' =", 
                                 time_standardisation, ":", "RoC values will be reported as disimilarity per", 
                                 time_standardisation, "years."))
  
  if (working_units != "levels" && time_standardisation != 
      bin_size) {
    RUtilpol::output_comment(paste("RoC values will be reported in different units than size of bin.", 
                                   "Recommend to keep 'time_standardisation'", "and 'bin_size' as same values"))
  }
  
  if (isTRUE(standardise)) {
    RUtilpol::output_comment(paste("Data will be standardise in each Working unit to", 
                                   n_individuals, "or the lowest number detected in dataset"))
    if (rand < 100) {
      RUtilpol::output_warning(paste("'standardise' was selected as 'TRUE' with low number of replication.", 
                                     "Recommend to increase 'rand'"))
    }
  }
  
  
  # extract data
  data_extract <- RRatepol:::extract_data(data_community_extract = data_source_community, 
                               data_age_extract = data_source_age, age_uncertainty = age_uncertainty, 
                               verbose = verbose)
  
  
  if (ncol(data_extract$community) == 1 && isTRUE(tranform_to_proportions)) {
    RUtilpol::output_warning(msg = paste("Community data has only 1 variable and `tranform_to_proportions`", 
                                         "is set to `TRUE`.", "This will result in 0 RoC.", 
                                         "Therefore, `tranform_to_proportions` will be set to `FALSE`"))
    tranform_to_proportions <- FALSE
  }
  
  
  # smooth data
  if (smooth_method != "none") {
    data_smooth <- RRatepol:::smooth_community_data(data_source_smooth = data_extract, 
                                         smooth_method = smooth_method, smooth_n_points = smooth_n_points, 
                                         smooth_n_max = smooth_n_max, smooth_age_range = smooth_age_range, 
                                         round_results = standardise, verbose = verbose)
  }
  else {
    data_smooth <- data_extract
  }
  
  data_work <- RRatepol:::reduce_data(data_source_reduce = data_smooth)
  
  if (isTRUE(verbose)) {
    RRatepol:::check_data(data_work)
  }
  if (is.null(age_uncertainty) && isFALSE(standardise) && isFALSE(is.null(rand)) && 
      (working_units == "levels" || bin_selection == "first")) {
    if (isTRUE(verbose)) {
      RUtilpol::output_comment(msg = paste("There is no need for randomisation.", 
                                           "Changing `rand` to NULL"))
    }
    rand <- NULL
  }
  
  data_prepared <- RRatepol:::prepare_data(data_source_prep = data_work, 
                                working_units = working_units, bin_size = bin_size, number_of_shifts = number_of_shifts, 
                                rand = rand)
  data_to_run <- RUtilpol::flatten_list_by_one(data_prepared)
  if (isTRUE(verbose)) {
    RUtilpol::output_heading(msg = "Start of estimation", 
                             size = "h2")
    RUtilpol::output_comment(msg = paste("Number of estimation set to", 
                                         length(data_to_run)))
  }
  
  if (isTRUE(use_parallel)) {
    if (methods::is(use_parallel, "numeric")) {
      n_cores <- as.numeric(use_parallel)
    }
    else {
      n_cores <- parallel::detectCores()
    }
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, {
      library("tidyverse")
      library("RRatepol")
    })
  }
  else {
    cl <- NULL
  }
  
  # modified run_interation injected
  result_list <- pbapply::pblapply(X = data_to_run, FUN = run_iterationMod, 
                                   cl = cl, bin_selection = bin_selection, standardise = standardise, 
                                   n_individuals = n_individuals, tranform_to_proportions = tranform_to_proportions, 
                                   dissimilarity_coefficient = dissimilarity_coefficient, 
                                   time_standardisation = time_standardisation, verbose = verbose, ...)
  if (isFALSE(is.null(cl))) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  gc(verbose = FALSE)
  
  return(result_list)
  
  # RRatepol summarises the data including upper and lower CIs across iterations
  # custom function to do the equivalent across more columns
  template = result_list[[1]]
  
  # group data by variable
  dataGroup = lapply(2:ncol(result_list[[1]]), function(n){
          sapply(result_list, function(x){x[[n]]})
  })
  names(dataGroup) = colnames(template)[-1]
  
  # estimate median, 2.5% and 97.5% quantiles
  # most of these, including ROC, are not normally distributed?
  dataProc = do.call("cbind", lapply(seq_along(dataGroup), function(n){
    x = t(apply(dataGroup[[n]], 1, function(y){c(median(y, na.rm=TRUE),
                                               quantile(y, 0.025, na.rm=TRUE),
                                               quantile(y, 0.975, na.rm=TRUE))}))
    colnames(x) = paste0(names(dataGroup)[n], c("_med", "_lci", "_uci"))
    return(x)
  }))
  
  end_time <- Sys.time()
  time_duration <- end_time - start_time
  RUtilpol::output_heading(paste("RRatepol finished", end_time, 
                                 "taking", round(time_duration, 2), units(time_duration)), 
                           size = "h1")
  return(dataProc)
}
