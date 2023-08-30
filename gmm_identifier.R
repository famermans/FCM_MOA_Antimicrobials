gmm_identifier <- function(flow_set, PhenoGMM_model, use.exprs = TRUE)
{
  fcs_PGMM_param <- flow_set[, attributes(PhenoGMM_model[[1]])$param]
  data_table <- flowCore::fsApply(fcs_PGMM_param, FUN = function(x) {
    data.table(gmm_id = mclust::predict.Mclust(PhenoGMM_model[[2]], x)$classification, x)
  }, use.exprs = use.exprs, simplify = FALSE)
}