MatthewsCC <- function(data, lev = NULL, model = NULL)
{
  out <- c(mltools::mcc(preds = data$pred, actuals = data$obs))
  names(out) <- c("MCC")
  out
}