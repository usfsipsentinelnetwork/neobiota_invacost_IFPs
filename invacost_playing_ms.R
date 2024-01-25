library(invacost)
library(dplyr)
data(invacost)

# all, worldwide
#x2 <- expandYearlyCosts(
#  invacost %>%
#    filter(!(Probable_starting_year_adjusted %in% c(NA,""))) %>%
#    filter(!(Probable_ending_year_adjusted %in% c(NA,""))),
#  startcolumn='Probable_starting_year_adjusted',
#  endcolumn='Probable_ending_year_adjusted')
#z1<-summarizeCosts(x2)
#z1$average.cost.per.period

# North America
#x2.na <- expandYearlyCosts(
#  invacost %>%filter(grepl("North America", Geographic_region)) %>%
#    filter(!(Probable_starting_year_adjusted %in% c(NA,""))) %>%
#    filter(!(Probable_ending_year_adjusted %in% c(NA,""))),
#  startcolumn='Probable_starting_year_adjusted',
#  endcolumn='Probable_ending_year_adjusted')
#z1.na<-summarizeCosts(x2.na)
#z1.na$average.cost.per.period

########################

#invacost$Geographic_region %>% unique

#invacost$Environment %>% unique
#invacost$Environment_IAS %>% unique
#invacost$Impacted_sector %>% unique
#invacost$Overlap %>% summary
#invacost$HabitatVerbatim %>% unique
#invacost$Habitat %>% table
#invacost$Type_of_cost_merged %>% unique
#invacost$Type_of_cost %>% unique


##################
# MODEL FUNCTION - removed GAM

modelCosts2 <- function(costdb,
                       cost.column = "Cost_estimate_per_year_2017_USD_exchange_rate",
                       year.column = "Impact_year",
                       cost.transf = "log10",
                       in.millions = TRUE,
                       confidence.interval = 0.95,
                       minimum.year = 1960, 
                       maximum.year = max(costdb[, year.column]), 
                       final.year = max(costdb[, year.column]), 
                       # models = c("ols.linear", 
                       #            "ols.quadratic",
                       #            "robust.linear",
                       #            "robust.quadratic",
                       #            "gam",
                       #            "mars",
                       #            "quantile"),
                       incomplete.year.threshold = NULL, # Changed default behaviour 2020.11.18
                       incomplete.year.weights = NULL,
                       gam.k = -1,
                       mars.nprune = NULL,
                       ...
)
{
  
  # Argument checking -------------------------------------------------------
  if(nrow(costdb) == 0)
  {
    stop("costdb is an empty table.\n")
  }
  dots <- list(...)
  # Checking if deprecated mars.nk argument was provided
  if("mars.nk" %in% names(dots))
  {
    stop("Argument mars.nk was specified. If you are looking to reduce model size, please use mars.nprune instead of mars.nk.")
  }
  
  if(any(is.na(costdb[, cost.column])))
  {
    warning("There were NA values in the cost column, they have been removed.\n")
    costdb <- costdb[-which(is.na(costdb[, cost.column])), ]
  }
  
  
  
  # if(any(!(models %in% c("ols.linear", 
  #                        "ols.quadratic", 
  #                        "gam",
  #                        "mars",
  #                        "quantile",
  #                        "robust.linear",
  #                        "robust.quadratic"))))
  # {
  #   stop(paste0("Inadequate model(s) specified:'",
  #               paste(models[which(!(models %in% c("ols.linear", 
  #                                                  "ols.quadratic", 
  #                                                  "gam",
  #                                                  "mars",
  #                                                  "quantile",
  #                                                  "robust.linear",
  #                                                  "robust.quadratic")))],
  #                     collapse = "', '"),
  #               "', please choose among 'ols.linear', 'ols.quadratic', 'robust.linear', 'robust.quadratic', 'gam', 'mars' and 'quantile'"))
  # }
  
  if(maximum.year < minimum.year)
  {
    stop("maximum.year is lower than minimum.year.\n")
  }
  
  if(is.null(incomplete.year.threshold))
  {
    incomplete.year.threshold <- maximum.year + 1
  }
  
  if(is.null(cost.transf))
  {
    cost.transf <- "none"
  }
  
  if(any(costdb[, year.column] < minimum.year))
  {
    warning(paste0("There are ",  length(unique(costdb$Cost_ID[which(costdb[, year.column] < minimum.year)])),
                   " cost values for periods earlier than ",
                   minimum.year, ", which will be removed.\n"))
    costdb <- costdb[-which(costdb[, year.column] < minimum.year), ]
  }
  
  if(any(costdb[, year.column] > maximum.year))
  {
    warning(paste0("There are cost values for periods later than ",
                   maximum.year, ": ",
                   length(unique(costdb$Cost_ID[which(costdb[, year.column] > maximum.year)])),
                   " different cost estimate(s).\nTheir values later than ",
                   maximum.year,
                   " will be removed.\n"))
    costdb <- costdb[-which(costdb[, year.column] > maximum.year), ]
  }
  
  if(nrow(costdb) == 0)
  {
    stop("There are no costs in the database that are between minimum.year and maximum.year")
  }
  
  
  
  if(final.year > maximum.year |
     final.year > max(costdb[, year.column], na.rm = TRUE))
  {
    warning(paste0("The final year is beyond the range of data,
    which creates an extrapolation situation. Be careful, the models included 
    here may not be realistic for extrapolations, because they do not include
    underlying drivers which my result in non-linearities over time."))
  }
  
  parameters <- list(cost.transformation = cost.transf,
                     incomplete.year.threshold = incomplete.year.threshold,
                     in.millions = in.millions,
                     confidence.interval = confidence.interval,
                     minimum.year = min(costdb[, year.column], na.rm = TRUE), 
                     maximum.year = max(costdb[, year.column], na.rm = TRUE), 
                     final.year = final.year,
                     gam.k = gam.k)
  
  
  yeargroups <- dplyr::group_by(costdb,
                                get(year.column)) 
  
  yearly.cost <- dplyr::summarise(yeargroups, 
                                  Annual.cost = sum(get(cost.column)))
  names(yearly.cost)[1] <- "Year"
  
  if(!is.null(incomplete.year.weights))
  {
    if(!all(yearly.cost$Year %in% names(incomplete.year.weights)))
    {
      stop("The vector provided in incomplete.year.weights does not have all the
           years in the range of data.")
    } else
    {
      incomplete.year.weights <- incomplete.year.weights[names(incomplete.year.weights) %in%
                                                           yearly.cost$Year]
      incomplete.year.weights <- incomplete.year.weights[match(yearly.cost$Year,
                                                               names(incomplete.year.weights))]
    }
  }
  
  if(in.millions)
  {
    yearly.cost$Annual.cost <- yearly.cost$Annual.cost / 1e6
  }
  
  
  
  if(cost.transf != "none")
  {
    yearly.cost$transf.cost <- do.call(cost.transf, list(yearly.cost$Annual.cost))
  } else
  {
    yearly.cost$transf.cost <- yearly.cost$Annual.cost
  }
  
  if(any(yearly.cost[, "Year"] >= incomplete.year.threshold))
  {
    message(paste0(length(which(yearly.cost[, "Year"] >= incomplete.year.threshold)),
                   " years will not be included in model calibrations because\n",
                   "they occurred later than incomplete.year.threshold (", incomplete.year.threshold,
                   ")\n"))
    yearly.cost$Calibration <- ifelse(yearly.cost$Year < incomplete.year.threshold,
                                      "Included", "Excluded")
    yearly.cost.calibration <- yearly.cost[-which(yearly.cost[, "Year"] >= incomplete.year.threshold), ]
    if(!is.null(incomplete.year.weights))
    {
      incomplete.year.weights <- incomplete.year.weights[-which(names(incomplete.year.weights) >= incomplete.year.threshold)]
    }
  } else
  {
    yearly.cost$Calibration <- "Included"
    yearly.cost.calibration <- yearly.cost
  }
  # For nicer graphs
  yearly.cost$Calibration <- factor(yearly.cost$Calibration, levels = c("Excluded",
                                                                        "Included"))
  if(!is.null(incomplete.year.weights))
  {
    yearly.cost.calibration <- data.frame(yearly.cost.calibration,
                                          incomplete.year.weights = incomplete.year.weights)
  }
  
  model.RMSE <- array(NA, dim = c(9, 2),
                      dimnames = list(c("ols.linear",
                                        "ols.quadratic",
                                        "robust.linear",
                                        "robust.quadratic",
                                        "mars",
                                        "gam",
                                        "qt0.1",
                                        "qt0.5",
                                        "qt0.9"),
                                      c("RMSE.calibration", 
                                        "RMSE.alldata")))
  
  # Prediction years correspond to the entire range provided by the user
  prediction.years <- data.frame(Year = minimum.year:maximum.year)
  
  
  # Ordinary Least Square (OLS) regression ----------------------------------
  
  message("\n --- Computing OLS regressions\n")
  # Ordinary least square - linear effect
  ols.linear <- lm(transf.cost ~ Year, data = yearly.cost.calibration,
                   weights = incomplete.year.weights)
  
  pred.ols.linear <- predict(ols.linear, 
                             prediction.years,
                             interval = "confidence",
                             level = confidence.interval)
  rownames(pred.ols.linear) <- prediction.years[, 1]
  
  model.RMSE["ols.linear", "RMSE.calibration"] <- sqrt(mean(residuals(ols.linear)^2))
  model.RMSE["ols.linear", "RMSE.alldata"] <- sqrt(
    mean((pred.ols.linear[match(yearly.cost$Year, rownames(pred.ols.linear)), "fit"] -
            yearly.cost$transf.cost)^2))
  
  # Calculation of heteroscedastic- and autocorrelation-robust variance covariance matrix of estimators for errors
  vcov.HAC.linear <- sandwich::vcovHAC(ols.linear)
  
  # Calculating 95% confidence intervals based on robust variance covariance matrix
  modelmatrix.linear <- stats::model.matrix(~ Year,
                                            data = prediction.years)
  
  # Variance of prediction years
  var.years.linear <- modelmatrix.linear %*% vcov.HAC.linear %*% t(modelmatrix.linear)
  # Standard errors
  se.years.linear <- sqrt(diag(var.years.linear))
  
  # Confidence intervals
  pred.ols.linear[, "lwr"] <-  pred.ols.linear[, "fit"] -
    se.years.linear * qt(confidence.interval + 
                           (1 - confidence.interval) / 2,
                         df = ols.linear$df.residual)
  pred.ols.linear[, "upr"] <-  pred.ols.linear[, "fit"] +
    se.years.linear * qt(confidence.interval + 
                           (1 - confidence.interval) / 2,
                         df = ols.linear$df.residual)
  
  
  
  # OLS - quadratic effect
  ols.quadratic <- lm(transf.cost ~ Year + I(Year^2), data = yearly.cost.calibration,
                      weights = incomplete.year.weights)
  pred.ols.quadratic <- predict(ols.quadratic, 
                                prediction.years,
                                interval = "confidence",
                                level = confidence.interval)
  rownames(pred.ols.quadratic) <- prediction.years[, 1]
  
  model.RMSE["ols.quadratic", "RMSE.calibration"] <- sqrt(mean(residuals(ols.quadratic)^2))
  model.RMSE["ols.quadratic", "RMSE.alldata"] <- sqrt(
    mean((pred.ols.quadratic[match(yearly.cost$Year, rownames(pred.ols.quadratic)), "fit"] -
            yearly.cost$transf.cost)^2))
  
  # Calculation of heteroscedastic- and autocorrelation-robust variance covariance matrix of estimators for errors
  vcov.HAC.quadratic <- sandwich::vcovHAC(ols.quadratic)
  
  # Calculating 95% confidence intervals based on robust variance covariance matrix
  modelmatrix.quadric <- stats::model.matrix(~ Year + I(Year^2),
                                             data = prediction.years)
  # Variance of prediction years
  var.years.quadratic <- modelmatrix.quadric %*% vcov.HAC.quadratic %*% t(modelmatrix.quadric)
  # Standard errors
  se.years.quadratic <- sqrt(diag(var.years.quadratic))
  
  
  # Confidence intervals
  pred.ols.quadratic[, "lwr"] <-  pred.ols.quadratic[, "fit"] -
    se.years.quadratic * qt(confidence.interval + 
                              (1 - confidence.interval) / 2,
                            df = ols.quadratic$df.residual)
  pred.ols.linear[, "upr"] <-  pred.ols.linear[, "fit"] +
    se.years.quadratic * qt(confidence.interval + 
                              (1 - confidence.interval) / 2,
                            df = ols.quadratic$df.residual)
  
  
  
  # Robust regression -------------------------------------------------------
  # Robust regression - Linear effect
  message("\n --- Computing robust regressions\n")
  robust.linear <- robustbase::lmrob(transf.cost ~ Year, data = yearly.cost.calibration, 
                                     weights = incomplete.year.weights)
  pred.robust.linear <- try(predict(robust.linear, 
                                    prediction.years,
                                    interval = "confidence", 
                                    level = confidence.interval),
                            silent = TRUE)
  if("try-error" %in% class(pred.robust.linear)) 
  {
    warning("Could not estimate confidence interval for robust linear regression")
    pred.robust.linear <- data.frame(
      fit = predict(robust.linear,
                    newdata = prediction.years),
      lwr = NA, upr = NA)
  }
  rownames(pred.robust.linear) <- prediction.years[, 1]
  
  
  
  model.RMSE["robust.linear", "RMSE.calibration"] <- sqrt(mean(residuals(robust.linear)^2))
  model.RMSE["robust.linear", "RMSE.alldata"] <- 
    sqrt(mean((pred.robust.linear[match(yearly.cost$Year, rownames(pred.robust.linear)), 
                                  "fit"] - yearly.cost$transf.cost)^2))
  
  
  
  # Robust regression - quadratic effect
  robust.quadratic <- robustbase::lmrob(transf.cost ~ Year + I(Year^2), data = yearly.cost.calibration, 
                                        weights = incomplete.year.weights,
                                        cov = ".vcov.w") # Covariance matrix estimated using asymptotic normality of the coefficients 
  # See ?lmrob and Koller & Stahel 2011 
  pred.robust.quadratic <- try(predict(robust.quadratic,
                                       prediction.years,
                                       interval = "confidence", 
                                       level = confidence.interval),
                               silent = TRUE)
  if("try-error" %in% class(pred.robust.quadratic)) 
  {
    warning("Could not estimate confidence interval for robust quadratic regression")
    pred.robust.quadratic <- data.frame(
      fit = predict(robust.quadratic,
                    newdata = prediction.years),
      lwr = NA, upr = NA)
  }
  
  rownames(pred.robust.quadratic) <- prediction.years[, 1]
  
  model.RMSE["robust.quadratic", "RMSE.calibration"] <- sqrt(mean(residuals(robust.quadratic)^2))
  model.RMSE["robust.quadratic", "RMSE.alldata"] <- sqrt(
    mean((pred.robust.quadratic[match(yearly.cost$Year, rownames(pred.robust.quadratic)), 
                                "fit"] - yearly.cost$transf.cost)^2))
  
  
  # Multiple Adapative Regression Splines -----------------------------------
  
  message("\n --- Computing MARS\n")
  mars <- earth::earth(transf.cost ~ Year, data = yearly.cost.calibration,
                       varmod.method = "lm",
                       # nk = mars.nk,
                       nprune = mars.nprune,
                       nfold = 5, 
                       ncross = 30, 
                       pmethod = "backward", # Would probably be better to use cross-validation but it does not work currently (I contacted the package author to fix this issue)
                       weights = incomplete.year.weights)
  
  pred.mars <- predict(mars,
                       prediction.years,
                       interval = "pint",
                       level = confidence.interval)
  rownames(pred.mars) <- prediction.years[, 1]
  
  model.RMSE["mars", "RMSE.calibration"] <- sqrt(mean(residuals(mars)^2))
  model.RMSE["mars", "RMSE.alldata"] <- sqrt(
    mean((pred.mars[match(yearly.cost$Year, rownames(pred.mars)), "fit"] -
            yearly.cost$transf.cost)^2))
  
  
  
  # Quantile regression -----------------------------------------------------
  
  message("\n --- Computing quantile regressions\n")
  qt0.1 <- quantreg::rq(transf.cost ~ Year, 
                        data = yearly.cost.calibration,
                        tau = 0.1,
                        weights = incomplete.year.weights)
  qt0.5 <- quantreg::rq(transf.cost ~ Year, 
                        data = yearly.cost.calibration,
                        tau = 0.5,
                        weights = incomplete.year.weights)
  qt0.9 <- quantreg::rq(transf.cost ~ Year, 
                        data = yearly.cost.calibration,
                        tau = 0.9,
                        weights = incomplete.year.weights)
  
  
  # quantreg sometimes throws errors in the prediction of confidence intervals
  # so we need to adapt the code
  
  pred.qt0.1 <- try(predict(qt0.1,
                            newdata = prediction.years,
                            interval = "confidence"),
                    silent = TRUE)
  if("try-error" %in% class(pred.qt0.1)) 
  {
    warning("Could not estimate confidence interval for quantile 0.1 regression")
    pred.qt0.1 <- data.frame(
      fit = predict(qt0.1,
                    newdata = prediction.years),
      lwr = NA, upr = NA)
  }
  pred.qt0.5 <- try(predict(qt0.5,
                            newdata = prediction.years,
                            interval = "confidence"),
                    silent = TRUE)
  if("try-error" %in% class(pred.qt0.5)) 
  {
    warning("Could not estimate confidence interval for quantile 0.5 regression")
    pred.qt0.5 <- data.frame(
      fit = predict(qt0.5,
                    newdata = prediction.years),
      lwr = NA, upr = NA)
  }
  pred.qt0.9 <- try(predict(qt0.9,
                            newdata = prediction.years,
                            interval = "confidence"),
                    silent = TRUE)
  if("try-error" %in% class(pred.qt0.9)) 
  {
    warning("Could not estimate confidence interval for quantile 0.9 regression")
    pred.qt0.9 <- data.frame(
      fit = predict(qt0.9,
                    newdata = prediction.years),
      lwr = NA, upr = NA)
  }
  colnames(pred.qt0.9) <- colnames(pred.qt0.5) <- colnames(pred.qt0.1) <- colnames(pred.ols.linear)
  rownames(pred.qt0.9) <- rownames(pred.qt0.5) <- rownames(pred.qt0.1) <- prediction.years[, 1]
  model.RMSE["qt0.1", "RMSE.calibration"] <- sqrt(mean(residuals(qt0.1)^2))
  model.RMSE["qt0.1", "RMSE.alldata"] <- sqrt(
    mean((pred.qt0.1[match(yearly.cost$Year, rownames(pred.qt0.1)), "fit"] -
            yearly.cost$transf.cost)^2))
  model.RMSE["qt0.5", "RMSE.calibration"] <- sqrt(mean(residuals(qt0.5)^2))
  model.RMSE["qt0.5", "RMSE.alldata"] <- sqrt(
    mean((pred.qt0.5[match(yearly.cost$Year, rownames(pred.qt0.5)), "fit"] -
            yearly.cost$transf.cost)^2))
  model.RMSE["qt0.9", "RMSE.calibration"] <- sqrt(mean(residuals(qt0.9)^2))
  model.RMSE["qt0.9", "RMSE.alldata"] <- sqrt(
    mean((pred.qt0.9[match(yearly.cost$Year, rownames(pred.qt0.9)), "fit"] -
            yearly.cost$transf.cost)^2))
  
  
  
  # Assembling predictions --------------------------------------------------
  
  
  message("\n --- Preparing the output objects\n")
  model.preds <- rbind.data.frame(data.frame(model = "OLS regression",
                                             Year = prediction.years$Year,
                                             Details = "Linear",
                                             pred.ols.linear),
                                  data.frame(model = "OLS regression",
                                             Year = prediction.years$Year,
                                             Details = "Quadratic",
                                             pred.ols.quadratic),
                                  data.frame(model = "Robust regression",
                                             Year = prediction.years$Year,
                                             Details = "Linear",
                                             pred.robust.linear),
                                  data.frame(model = "Robust regression",
                                             Year = prediction.years$Year,
                                             Details = "Quadratic",
                                             pred.robust.quadratic),
                                  data.frame(model = "MARS",
                                             Year = prediction.years$Year,
                                             Details = "",
                                             pred.mars),
                                  data.frame(model = "Quantile regression",
                                             Year = prediction.years$Year,
                                             Details = "Quantile 0.1",
                                             pred.qt0.1),
                                  data.frame(model = "Quantile regression",
                                             Year = prediction.years$Year,
                                             Details = "Quantile 0.5",
                                             pred.qt0.5),
                                  data.frame(model = "Quantile regression",
                                             Year = prediction.years$Year,
                                             Details = "Quantile 0.9",
                                             pred.qt0.9))
  
  
  # Model Summaries ---------------------------------------------------------
  
  
  # Creating the list containing the summary of model results
  testsummary <- list()
  # OLS
  testsummary$ols.linear$coeftest <- lmtest::coeftest(ols.linear, df = ols.linear$df.residual, vcov = vcov.HAC.linear)
  testsummary$ols.linear$r.squared <- summary(ols.linear)$r.squared
  testsummary$ols.linear$adjusted.r.squared <- summary(ols.linear)$adj.r.squared
  testsummary$ols.quadratic$coeftest <- lmtest::coeftest(ols.quadratic, df = ols.quadratic$df.residual, vcov = vcov.HAC.quadratic)
  testsummary$ols.quadratic$r.squared <- summary(ols.quadratic)$r.squared
  testsummary$ols.quadratic$adjusted.r.squared <- summary(ols.quadratic)$adj.r.squared
  # Robust
  testsummary$robust.linear <- summary(robust.linear)
  testsummary$robust.quadratic <- summary(robust.quadratic)
  # MARS
  testsummary$mars <- summary(mars)
  # Quantile
  testsummary$qt0.1 <- summary(qt0.1)
  testsummary$qt0.5 <- summary(qt0.5)
  testsummary$qt0.9 <- summary(qt0.9)
  
  class(testsummary) <- append("invacost.modelsummary", class(testsummary))
  
  
  
  
  # Preparing outputs -------------------------------------------------------
  
  # Formatting results for output object
  if(cost.transf == "log10")
  {
    # Transform log10 values back to actual US$
    model.preds[, c("fit", "lwr", "upr")] <- 
      apply(model.preds[, c("fit", "lwr", "upr")] ,
            2,
            function(x) 10^x)
    
    results <- list(cost.data = yearly.cost,
                    parameters = parameters, 
                    calibration.data = yearly.cost.calibration,
                    fitted.models = list(ols.linear = ols.linear,
                                         ols.quadratic = ols.quadratic,
                                         robust.linear = robust.linear,
                                         robust.quadratic = robust.quadratic,
                                         mars = mars,
                                         quantile = list(qt0.1 = qt0.1,
                                                         qt0.5 = qt0.5,
                                                         qt0.9 = qt0.9)),
                    estimated.annual.costs = model.preds,
                    model.summary = testsummary,
                    RMSE = model.RMSE,
                    final.year.cost = c(ols.linear = 
                                          unname(10^predict(ols.linear,
                                                            newdata = data.frame(Year = final.year))),
                                        ols.quadratic =  
                                          unname(10^predict(ols.quadratic,
                                                            newdata = data.frame(Year = final.year))),
                                        robust.linear = 
                                          unname(10^predict(robust.linear,
                                                            newdata = data.frame(Year = final.year))),
                                        robust.quadratic = 
                                          unname(10^predict(robust.quadratic,
                                                            newdata = data.frame(Year = final.year))),
                                        mars = 
                                          unname(10^predict(mars,
                                                            newdata = data.frame(Year = final.year))),
                                        quantile.0.1 = 
                                          unname(10^predict(qt0.1,
                                                            newdata = data.frame(Year = final.year))),
                                        quantile.0.5 = 
                                          unname(10^predict(qt0.5,
                                                            newdata = data.frame(Year = final.year))),
                                        quantile.0.9 = 
                                          unname(10^predict(qt0.9,
                                                            newdata = data.frame(Year = final.year)))))
  } else if(cost.transf == "none")
  {
    results <- list(cost.data = yearly.cost,
                    parameters = parameters, 
                    calibration.data = yearly.cost.calibration,
                    fitted.models = list(linear = ols.linear, 
                                         quadratic = ols.quadratic, 
                                         robust.linear = robust.linear,
                                         robust.quadratic = robust.quadratic,
                                         mars = mars,
                                         quantile = list(qt0.1 = qt0.1,
                                                         qt0.5 = qt0.5,
                                                         qt0.9 = qt0.9)),
                    estimated.annual.costs = model.preds,
                    model.summary = testsummary,
                    RMSE = model.RMSE,
                    final.year.cost = c(ols.linear =  
                                          unname(predict(ols.linear,
                                                         newdata = data.frame(Year = final.year))),
                                        ols.quadratic =  
                                          unname(predict(ols.quadratic,
                                                         newdata = data.frame(Year = final.year))),
                                        robust.linear = 
                                          unname(predict(robust.linear,
                                                         newdata = data.frame(Year = final.year))),
                                        robust.quadratic = 
                                          unname(predict(robust.quadratic,
                                                         newdata = data.frame(Year = final.year))),
                                        mars = 
                                          unname(predict(mars,
                                                         newdata = data.frame(Year = final.year))),
                                        quantile.0.1 = 
                                          unname(predict(qt0.1,
                                                         newdata = data.frame(Year = final.year))),
                                        quantile.0.5 = 
                                          unname(predict(qt0.5,
                                                         newdata = data.frame(Year = final.year))),
                                        quantile.0.9 = 
                                          unname(predict(qt0.9,
                                                         newdata = data.frame(Year = final.year)))))
  } else
  {
    results <- list(input.data = costdb,
                    cost.data = yearly.cost,
                    parameters = parameters, 
                    calibration.data = yearly.cost.calibration,
                    fitted.models = list(linear = ols.linear, # Inconsistent name, should be corrected (but need to check generic functions)
                                         quadratic = ols.quadratic, # Inconsistent name, should be corrected (but need to check generic functions)
                                         robust.linear = robust.linear,
                                         robust.quadratic = robust.quadratic,
                                         mars = mars,
                                         quantile = list(qt0.1 = qt0.1,
                                                         qt0.5 = qt0.5,
                                                         qt0.9 = qt0.9)),
                    model.summary = testsummary,
                    RMSE = model.RMSE)
  }
  class(results) <- append("invacost.costmodel", class(results))
  return(results)
}

####
# HELPER BREAKS FUNCTION

my_breaks <- function(x) {
  # x <- c(1,10000)
  rr <- log10(max(x)) - log10(min(x))
  
  #(10^(rr %/% 6))^(-15:15)
  
  if (rr <= 4) {
    (10)^(-15:15)
  } else if (rr <= 16) {
    (100)^(-15:15)
  } else {
    (100000000)^(-15:15)
  }
}

my_limits <- function(x) { c(1, ifelse(x[2]<=100000000,x[2],100000000)) }

################
# PLOT FUNCTION

plot.invacost.costmodel.custom <- function(y,
                                    plot.breaks = 10^(-15:15),
                                    plot.type = "facets",
                                    models = c("ols.linear", 
                                               "ols.quadratic", 
                                               "robust.linear",
                                               "robust.quadratic",
                                               "gam",
                                               "mars",
                                               "quantile"),
                                    evaluation.metric = FALSE,
                                    graphical.parameters = NULL,
                                    params = y[[1]]$models$parameters,
                                    pred.cutoff = 2010,
                                    nc=NULL,
                                    ...) {
  
  models <- c("ols.linear", "ols.quadratic")
  
  # 1. We create the ggplot here ---------------
  p <- ggplot()
  
  # Setting up graphical.parameters
  if(is.null(graphical.parameters)) {
    # 2a. If users do not specify graphical parameters we create them here ---------------
    p <- p + 
      ylab(paste0("Annual cost in $ US ", 
                  ifelse(params$in.millions, 
                         "millions",
                         ""), "(2023 Equiv.)\n")) +
      xlab("\nYear") +
      theme_bw()
    if(params$cost.transformation == "log10") {
      # 3a. We define axes here for log-transformed data ---------------
      p <- p +
        scale_y_log10(
          #breaks = plot.breaks,
          breaks = my_breaks,
          labels = scales::comma,
          limits = my_limits) +
        annotation_logticks(sides = "l")
    } else if(params$cost.transformation == "none") {
      # 3b. We define axes here for untransformed data ---------------
      p <- p +
        scale_y_continuous(labels = scales::comma)
    } else {
      stop("If you made a manual transformation (other than log10), then you
           will have to make the plot by yourself.")
    }
  } else {
    # Workaround for retrocompatibility
    if(graphical.parameters == "manual") {
      graphical.parameters <- NULL
    }
    # 4. Adding user-defined parameters to the plot ---------------
    p <- p + graphical.parameters
  }
  
  # y <- to.plot
  # i <- 1
  # evaluation.metric <- F
  
  mod.preds.long <- NULL
  
  cost.data.long <- NULL

  for (i in 1:length(y)) {
    
    x <- y[[i]]$models
    
    if(is.null(x)) {
      cost.data.long <- bind_rows(cost.data.long, bind_cols(Year=NA, transf.cost=NA, whichone=names(y)[i]))
      mod.preds.long <- bind_rows(mod.preds.long, bind_cols(model=NA,Year=NA,Details=NA,fit=NA,lwr=NA,upr=NA, whichone=names(y)[i]))
      next
    }
    
    # Changing order of factors for points 
    x$cost.data$Calibration <- factor(x$cost.data$Calibration, 
                                      levels = c("Included", "Excluded"))
    
    cost.data.long <- rbind(cost.data.long, cbind(x$cost.data, whichone=names(y)[i]))
    
    # Preparing model predictions for plots
    model.preds <- x$estimated.annual.costs
    model.preds$Model <- as.character(model.preds$model)
    model.preds$Model[model.preds$model == "OLS regression" & 
                        model.preds$Details == "Linear"] <- "OLS linear regression"
    model.preds$Model[model.preds$model == "OLS regression" 
                      & model.preds$Details == "Quadratic"] <- "OLS quadratic regression"
    model.preds$Model[model.preds$model == "Robust regression" & 
                        model.preds$Details == "Linear"] <- "Robust linear regression"
    model.preds$Model[model.preds$model == "Robust regression" & 
                        model.preds$Details == "Quadratic"] <- "Robust quadratic regression"
    model.preds$Model[model.preds$model == "Quantile regression"] <-
      paste0(model.preds$Details[model.preds$model == "Quantile regression"],
             " regression")
  
    # Ordering model names
    model.preds$Model <- factor(model.preds$Model,
                                levels = c("OLS linear regression", 
                                           "OLS quadratic regression",
                                           "Robust linear regression",
                                           "Robust quadratic regression",
                                           "MARS",
                                           "GAM",
                                           paste("Quantile", c(0.1, 0.5, 0.9), "regression")))
    model.preds$model <- factor(model.preds$model,
                                levels = c("OLS regression", "Robust regression",
                                           "GAM", "MARS", "Quantile regression"))
    
    # Limiting plots to user selected
    #Relabel models parameter to match plot labeling from above
    models[models=="ols.linear"] <- "OLS linear regression"
    models[models=="ols.quadratic"] <- "OLS quadratic regression"
    models[models=="gam"] <- "GAM"
    models[models=="mars"] <- "MARS"
    models <- rep(models,1+2*(models=="quantile"))
    models[models=="quantile"] <- c("Quantile 0.1 regression","Quantile 0.5 regression","Quantile 0.9 regression")
    models[models=="robust.linear"] <- "Robust linear regression"
    models[models=="robust.quadratic"] <- "Robust quadratic regression"
    model.preds <- model.preds[model.preds$Model %in% models,]
    
    mod.preds.long <- rbind(mod.preds.long, cbind(model.preds, whichone=names(y)[i]))

    if(evaluation.metric) {
      model.rmse <- x$RMSE[, 1]
      
      #Relabel models parameter to match plot labeling from above
      names(model.rmse)[names(model.rmse) == "ols.linear"] <- "OLS linear"
      names(model.rmse)[names(model.rmse) == "ols.quadratic"] <- "OLS quadratic"
      names(model.rmse)[names(model.rmse) == "gam"] <- "GAM"
      names(model.rmse)[names(model.rmse) == "mars"] <- "MARS"
      names(model.rmse)[names(model.rmse) == "robust.linear"] <- "Robust linear"
      names(model.rmse)[names(model.rmse) == "robust.quadratic"] <- "Robust quadratic"
      
      model.rmse <- model.rmse[-grep("qt", names(model.rmse))]
    }
  
  }
  # Creating a colourblind palette (Wong 2011)
  # to best distinguish models
  alpha <- round(.8 * 255)
  cols <- c(`OLS linear regression` = rgb(86, 180, 233, alpha = alpha,
                                          maxColorValue = 255), # Sky blue
            `OLS quadratic regression` = rgb(230, 159, 0, alpha = alpha,
                                             maxColorValue = 255), # Orange
            `Robust linear regression` = rgb(0, 114, 178, alpha = alpha,
                                             maxColorValue = 255), # Blue
            `Robust quadratic regression` = rgb(213, 94, 0, alpha = alpha,
                                                maxColorValue = 255), # Vermillion
            `MARS` = rgb(204, 121, 167, alpha = alpha,
                         maxColorValue = 255), # Reddish purple
            `GAM` = rgb(0, 158, 115, alpha = alpha,
                        maxColorValue = 255), # Bluish green
            `Quantile 0.5 regression` = grey(0.5, alpha = alpha / 255),
            `Quantile 0.1 regression` = grey(0.25, alpha = alpha / 255),
            `Quantile 0.9 regression` = grey(0, alpha = alpha / 255)
  )
  
  if(plot.type == "single") {
    # 5. Single plot --------------------
    p <-
      p +
      geom_point(data = x$cost.data, 
                 aes_string(x = "Year",
                            y = "Annual.cost",
                            shape = "Calibration"),
                 col = grey(.4)) +
      geom_line(data = model.preds, 
                aes_string(x = "Year",
                           y = "fit",
                           col = "Model"),
                size = 1) +
      scale_discrete_manual(aesthetics = "col",
                            values = cols)
    
    if(evaluation.metric) {
      p <-
        p + labs(tag = paste0("RMSE\n",
                              paste(names(model.rmse), "/", round(model.rmse, 3), collapse = "\n"))) +
        theme(plot.tag = element_text(hjust = 1, vjust = 0),
              plot.tag.position = c(1, 0.05))
    }
    
    
  } else if(plot.type == "facets") { 
    # 6. Facet plot --------------------
    
    #mod.preds.list
    #str(mod.preds.list)
    #str(model.preds)
    
    cost.data.long$whichone <- factor(cost.data.long$whichone, levels=names(y))
    mod.preds.long$whichone <- factor(mod.preds.long$whichone, levels=names(y))
    
    ymx <- 100000000#max(na.omit(c(cost.data.long$Annual.cost,mod.preds.long$fit)))
    
    cost.data.letters <- data.frame(whichone=unique(cost.data.long$whichone))%>%
      mutate(letter = LETTERS[as.integer(whichone)])
    
    if(!is.null(nc)) {
      cost.data.long <- cost.data.long %>% left_join(nc)
      mod.preds.long <- mod.preds.long %>% left_join(nc)
      cost.data.letters <- cost.data.letters %>% left_join(nc)
    }
    
    p <-
      p +
      geom_point(data = cost.data.long,                       # now faceted
                 aes_string(x = "Year",
                            y = "Annual.cost",
                            shape = "Calibration"),
                 col = grey(.4)) +
      geom_line(data = mod.preds.long,                        # faceted
                aes_string(x = "Year",
                           y = "fit",
                           col = "Model"),
                size = 1) +
      geom_ribbon(data = mod.preds.long,                      # faceted
#                  aes_string(x = "Year",
#                             ymin = "lwr",
#                             ymax = "upr",
#                             group = "Details"),
                  aes(x=Year, ymin=ifelse(lwr<1, 1, lwr), ymax=ifelse(upr>ymx, ymx, upr), group=Details),#, fill=ifelse(Year <=2010, grey(1), 'pink')),
                  alpha = .1) + 
      geom_text(data=cost.data.letters,aes(x=1969, y=20000000, label=letter))+
      scale_discrete_manual(aesthetics = "col",
                            values = cols)
    
    if (is.null(nc))
      p <- p +facet_wrap(~ whichone,ncol=3)
    else
      p <- p + facet_grid(where ~ what)

#      facet_wrap (~ whichone,
#                  ncol=3,
#                  labeller = label_both#,
                  #scales = "free_y"
#                  ) +

    
    
    
    message("Note that MARS error bands are prediction intervals and not confidence interval (see ?plot.invacost.costmodel)\n")
    
    
    if(evaluation.metric) {
      p <-
        p + labs(tag = paste0("RMSE\n",
                              paste(names(model.rmse), "/", round(model.rmse, 3), collapse = "\n"))) +
        theme(plot.tag = element_text(hjust = 1, vjust = 0),
              plot.tag.position = c(0.95, 0.05))
    }
  }
  return(p)
}

##############################
### generic filtering function
###
inva.filt <- function(db=invacost,
                      kingdoms="Fungi",
                      phyla="Arthropoda|Nematoda",
                      hab=c("1","2","3"),
                      damage_type="All",
                      continent="All",
                      min=1960,
                      max=2010,
                      expansion.cutoff=2010,
                      nomodel=F)
  {
    x <- db %>%
      filter(grepl(kingdoms,Kingdom) | grepl(phyla,Phylum)) %>%
      filter(Habitat %in% c("1","2","3"))
    
    if (damage_type != "All")
      x <- x  %>% filter(Type_of_cost_merged == damage_type)
    
    if (continent != "All")
      x <- x %>% filter(grepl(continent, Geographic_region))
    
    number_of_records <- x$Reference_ID %>% unique %>% length
    
    #trend <- costTrendOverTime(x, minimum.year = min, maximum.year = max)
    
    
    
    x2 <- expandYearlyCosts(
      x %>%
        filter(!(Probable_starting_year_adjusted %in% c(NA,""))) %>%
        filter(!(Probable_ending_year_adjusted %in% c(NA,""))),
      startcolumn='Probable_starting_year_adjusted',
      endcolumn='Probable_ending_year_adjusted') %>% filter(Impact_year <= expansion.cutoff)
   z1<-summarizeCosts(x2)
   if(nomodel) {
     z2 <- NULL
     trend <- NULL
   } else {
     z2<-modelCosts2(x2, minimum.year=min, maximum.year=max)
     
     a <- z2$model.summary$ols.linear$coeftest['(Intercept)','Estimate']
     b <- z2$model.summary$ols.linear$coeftest['Year','Estimate']
     e <- z2$model.summary$ols.linear$coeftest['Year','Std. Error']
     
     trend <- c(a=a,
                b=b,
                e=e,
                b10=10^b,
                bmin10=(10^b)/(1.96*10^e),
                bmax10=(10^b)*(1.96*10^e),
                t2025=10^(a+b*2025),
                t2025mse=10^(a+(b+e)*2025), 
                t2025pse=10^(a+(b-e)*2025),
                t2050=10^(a+b*2050),
                t2050mse=10^(a+(b+e)*2050), 
                t2050pse=10^(a+(b-e)*2050))
     
     bq <- z2$model.summary$ols.linear$coeftest['Year','Estimate']
     eq <- z2$model.summary$ols.linear$coeftest['Year','Std. Error']
     
     trendq <- c(b=10^b, bmin=(10^b)/(1.96*10^e), bmax=(10^b)*(1.96*10^e))
     
     
   }
   
   
   return(list(
     data = x,
     trend = trend,
     n = number_of_records,
     expanded = x2,
     summary = z1,
     models = z2
   ))
}



### fungi, arthropods, nematodes - world 

world.all <- inva.filt() ## - all
world.d <- inva.filt(damage_type='Damage', min =1980) ### - damage
world.m <- inva.filt(damage_type='Management') ### - damage



# fungi, arthropods, nematodes - North America

north.america.all <- inva.filt(continent="North America", min =1988) ## - all
north.america.d <- inva.filt(continent="North America", damage_type='Damage', min =1988) ### - damage
north.america.m <- inva.filt(continent="North America", damage_type='Management', min =2000) ### - damage

# fungi, arthropods, nematodes - EU

eu.all <- inva.filt(continent="Europe", min=1960, max=2020) ## - all
eu.d <- inva.filt(continent="Europe", damage_type='Damage', min=1960, max=2010, nomodel=T) ### - damage
eu.m <- inva.filt(continent="Europe", damage_type='Management', min=1960, max=2020) ### - damage

plot(eu.all$models)

# elsewhere

# rest of world
elsewhere<-invacost %>% filter(!grepl("North America|Europe", Geographic_region, perl=T))
elsewhere$Geographic_region %>% unique

rest.of.world.all <- inva.filt(elsewhere)
rest.of.world.d <- inva.filt(elsewhere, damage_type='Damage',min=1980)
rest.of.world.m <- inva.filt(elsewhere, damage_type='Management')

# trends
world.all$trend
world.all$models$model.summary$ols.linear
10^(-58.2903770+(0.0307258)*2050)
10^(-58.2903770+(0.0307258-.0028136)*2050)
10^(-58.2903770+(0.0307258+.0028136)*2050)

world.all$trend
north.america.all$trend
eu.all$trend
rest.of.world.all$trend

world.all$trend

# put them together

to.plot <- list(world.all, world.d, world.m, north.america.all, north.america.d, north.america.m)

names(to.plot)<-
  c(
    "WorldAll",
    "WorldDamage",
    "WorldManagement",
    "NorthAmericaAll",
    "NorthAmericaDamage",
    "NorthAmericaManagement"
  )

to.plot2 <- list(world.all, world.d, world.m, north.america.all, north.america.d, north.america.m, eu.all,
                 eu.d,
                 eu.m)

names(to.plot2)<-
  c(
    "WorldAll",
    "WorldDamage",
    "WorldManagement",
    "NorthAmericaAll",
    "NorthAmericaDamage",
    "NorthAmericaManagement",
    "EuropeAll",
    "EuropeDamage",
    "EuropeManagement"
  )
#PLOT

plot.invacost.costmodel.custom(to.plot, models = c("ols.linear",  "ols.quadratic"))+
  theme(legend.position = 'none', axis.text.x = element_text(angle=45, vjust=1, hjust=1))

plot.invacost.costmodel.custom(to.plot2, models = c("ols.linear",  "ols.quadratic"))+
  theme(legend.position = 'none', axis.text.x = element_text(angle=45, vjust=1, hjust=1))

### SAVE RESULTS

# annual costs

cbind(
  to.plot$WorldAll$summary$average.cost.per.period[,1:2],
  data.frame(
    WorldAll = to.plot$WorldAll$summary$average.cost.per.period$annual_cost*1.0213,
    WorldDamage = to.plot$WorldDamage$summary$average.cost.per.period$annual_cost*1.0213,
    WorldManagement = to.plot$WorldManagement$summary$average.cost.per.period$annual_cost*1.0213,
    NorthAmericaAll = to.plot$NorthAmericaAll$summary$average.cost.per.period$annual_cost*1.0213,
    NorthAmericaDamage = to.plot$NorthAmericaDamage$summary$average.cost.per.period$annual_cost*1.0213,
    NorthAmericaManagement = to.plot$NorthAmericaManagement$summary$average.cost.per.period$annual_cost*1.0213
  )
) %>% write.csv("Costs_Yearly_2024_Dollars_1960_2010.csv", row.names=F)

cbind(
  to.plot2$WorldAll$summary$average.cost.per.period[,1:2],
  data.frame(
    WorldAll = to.plot2$WorldAll$summary$average.cost.per.period$annual_cost*1.0213,
    WorldDamage = to.plot2$WorldDamage$summary$average.cost.per.period$annual_cost*1.0213,
    WorldManagement = to.plot2$WorldManagement$summary$average.cost.per.period$annual_cost*1.0213,
    NorthAmericaAll = to.plot2$NorthAmericaAll$summary$average.cost.per.period$annual_cost*1.0213,
    NorthAmericaDamage = to.plot2$NorthAmericaDamage$summary$average.cost.per.period$annual_cost*1.0213,
    NorthAmericaManagement = to.plot2$NorthAmericaManagement$summary$average.cost.per.period$annual_cost*1.0213,
    EuropeAll = to.plot2$EuropeAll$summary$average.cost.per.period$annual_cost*1.0213,
    EuropeDamage = to.plot2$EuropeDamage$summary$average.cost.per.period$annual_cost*1.0213,
    EuropeManagement = to.plot2$EuropeManagement$summary$average.cost.per.period$annual_cost*1.0213
  )
)  %>% write.csv("Costs_Total_2024_Dollars_includingEurope_1960_2010.csv", row.names=F)

# total cost

cbind(
  to.plot$WorldAll$summary$average.cost.per.period[,1:2],
  data.frame(
    WorldAll = to.plot$WorldAll$summary$average.cost.per.period$total_cost*1.0213,
    WorldDamage = to.plot$WorldDamage$summary$average.cost.per.period$total_cost*1.0213,
    WorldManagement = to.plot$WorldManagement$summary$average.cost.per.period$total_cost*1.0213,
    NorthAmericaAll = to.plot$NorthAmericaAll$summary$average.cost.per.period$total_cost*1.0213,
    NorthAmericaDamage = to.plot$NorthAmericaDamage$summary$average.cost.per.period$total_cost*1.0213,
    NorthAmericaManagement = to.plot$NorthAmericaManagement$summary$average.cost.per.period$total_cost*1.0213
  )
) %>% write.csv("Costs_Total_2024_Dollars_1960_2010.csv", row.names=F)



#####
### fungi, arthropods, nematodes - world 

world.all.2050 <- inva.filt(max=2050) ## - all
world.d.2050 <- inva.filt(damage_type='Damage', min =1980, max=2050) ### - damage
world.m.2050 <- inva.filt(damage_type='Management', max=2050) ### - damage

# fungi, arthropods, nematodes - North America

north.america.all.2050 <- inva.filt(continent="North America", min =1988, max=2050) ## - all
north.america.d.2050 <- inva.filt(continent="North America", damage_type='Damage', min =1988, max=2050) ### - damage
north.america.m.2050 <- inva.filt(continent="North America", damage_type='Management', min =2000, max=2050) ### - damage

# fungi, arthropods, nematodes - EU

eu.all.2050 <- inva.filt(continent="Europe", min=1960, max=2050) ## - all
eu.d.2050 <- inva.filt(continent="Europe", damage_type='Damage', min=1960, max=2050, nomodel=T) ### - damage
eu.m.2050 <- inva.filt(continent="Europe", damage_type='Management', min=1960, max=2050) ### - damage



# put them together

to.plot.2050 <- list(world.all.2050, world.d.2050, world.m.2050, north.america.all.2050, north.america.d.2050, north.america.m.2050)

names(to.plot.2050)<-
  c(
    "WorldAll",
    "WorldDamage",
    "WorldManagement",
    "NorthAmericaAll",
    "NorthAmericaDamage",
    "NorthAmericaManagement"
  )

to.plot.2050.2 <- list(world.all.2050, world.d.2050, world.m.2050)
names(to.plot.2050.2)<-
  c(
    "WorldAll",
    "WorldDamage",
    "WorldManagement"
  )

to.plot.2050.2$WorldDamage$models$fitted.models<-to.plot.2050.2$WorldDamage$models$fitted.models[[c(1,3)]]

to.plot.2050.3 <- list(world.all.2050,
                       world.d.2050,
                       world.m.2050,
                       north.america.all.2050, 
                       north.america.d.2050, 
                       north.america.m.2050,
                       eu.all.2050,
                       eu.d.2050,
                       eu.m.2050)

names(to.plot.2050.3)<-
  c(
    "WorldAll",
    "WorldDamage",
    "WorldManagement",
    "NorthAmericaAll",
    "NorthAmericaDamage",
    "NorthAmericaManagement",
    "EuropeAll",
    "EuropeDamage",
    "EuropeManagement"
  )

#PLOT

plot.invacost.costmodel.custom(to.plot.2050, models = c("ols.linear",  "ols.quadratic"))+
  theme(legend.position = 'none', axis.text.x = element_text(angle=45, vjust=1, hjust=1))


plot.invacost.costmodel.custom(to.plot.2050.3, models = c("ols.linear",  "ols.quadratic"))+
  theme(legend.position = 'none', axis.text.x = element_text(angle=45, vjust=1, hjust=1))


### SAVE RESULTS

cbind(
  to.plot.2050$WorldAll$summary$average.cost.per.period[,1:2],
  data.frame(
    WorldAll = to.plot.2050$WorldAll$summary$average.cost.per.period$annual_cost,
    WorldDamage = to.plot.2050$WorldDamage$summary$average.cost.per.period$annual_cost,
    WorldManagement = to.plot.2050$WorldManagement$summary$average.cost.per.period$annual_cost,
    NorthAmericaAll = to.plot.2050$NorthAmericaAll$summary$average.cost.per.period$annual_cost,
    NorthAmericaDamage = to.plot.2050$NorthAmericaDamage$summary$average.cost.per.period$annual_cost,
    NorthAmericaManagement = to.plot.2050$NorthAmericaManagement$summary$average.cost.per.period$annual_cost
  )
)


rest.of.world.all.2050 <- inva.filt(elsewhere,max=2050)
rest.of.world.d.2050 <- inva.filt(elsewhere, damage_type='Damage',min=1980,max=2050)
rest.of.world.m.2050 <- inva.filt(elsewhere, damage_type='Management',max=2050)

to.plot.2050.4 <- list(world.all.2050,
                       world.d.2050,
                       world.m.2050,
                       north.america.all.2050, 
                       north.america.d.2050, 
                       north.america.m.2050,
                       eu.all.2050,
                       eu.d.2050,
                       eu.m.2050,
                       rest.of.world.all.2050,
                       rest.of.world.d.2050,
                       rest.of.world.m.2050)

names(to.plot.2050.4)<-
  c(
    "WorldAll",
    "WorldDamage",
    "WorldManagement",
    "NorthAmericaAll",
    "NorthAmericaDamage",
    "NorthAmericaManagement",
    "EuropeAll",
    "EuropeDamage",
    "EuropeManagement",
    "AsiaGlobalSAll",
    "AsiaGlobalSDamage",
    "AsiaGlobalSManagement"
  )

name.combos <- data.frame(whichone=names(to.plot.2050.4),
                          where=ordered(c(rep("World",3),rep("N. Amer.", 3),rep("Europe",3),rep("Asia & G. South",3))),
                          what=ordered(rep(c("All","Damage","Management"),4)))
name.combos$where <- factor(name.combos$where, levels=c("World","N. Amer.","Europe","Asia & G. South"))
#name.combos$what <- factor(name.combos$what, levels=)
levels(name.combos$where)
plot.invacost.costmodel.custom(to.plot.2050.4, models = c("ols.linear",  "ols.quadratic"), nc=name.combos)+
  theme(legend.position = 'none', axis.text.x = element_text(angle=45, vjust=1, hjust=1))
  


world.all.2050$n
world.d.2050$n
world.m.2050$n
north.america.all.2050$n
north.america.d.2050$n
north.america.m.2050$n
eu.all.2050$n
eu.d.2050$n
eu.m.2050$n
rest.of.world.all$n
rest.of.world.d$n
rest.of.world.m$n

##########
rest.of.world.all$expanded$Genus %>% table
rest.of.world.all$expanded$Genus %>% table %>% (function (x) x/sum(x))
north.america.all$expanded$Genus %>% table
north.america.all$expanded$Genus %>% table %>% (function (x) x/sum(x))

eu.all$expanded$Genus %>% table
eu.all$expanded$Genus %>% table%>% (function (x) x/sum(x))

world.all$expanded$Genus %>% table

world.all$expanded$Reference_title[which(world.all$expanded$Genus == "Diverse/Unspecified")] %>% unique
rest.of.world.all$expanded$Reference_title[which(world.all$expanded$Genus == "Diverse/Unspecified")] %>% unique
north.america.all$expanded$Reference_title[which(world.all$expanded$Genus == "Diverse/Unspecified")] %>% unique
eu.all$expanded$Reference_title[which(world.all$expanded$Genus == "Diverse/Unspecified")] %>% unique


plot(rest.of.world$models)
