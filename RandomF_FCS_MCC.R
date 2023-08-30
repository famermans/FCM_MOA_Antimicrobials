RandomF_FCS_MCC <- function (x, sample_info, sample_col = "name", target_label, 
                             downsample = 0, classification_type = "sample", param = c("FL1-H", 
                                                                                       "FL3-H", "FSC-H", "SSC-H"), p_train = 0.75, seed = 777, 
                             cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", 
                             plot_fig = FALSE, method = "rf", metric) 
{
  set.seed(seed)
  if (cleanFCS == TRUE) {
    cat(paste0("-------------------------------------------------------------------------------------------------", 
               "\n"))
    cat(date(), paste0("--- Using the following parameters for removing errant collection events\n in samples:\n \n"))
    cat(paste0(param), "\n \n")
    cat(date(), paste0("--- Scatter parameters will be automatically excluded", 
                       "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------"))
    cat("\n", paste0("Please cite:", "\n"))
    cat("\n", paste0("Monaco et al., flowAI: automatic and interactive anomaly discerning tools for flow cytometry data,\n Bioinformatics, Volume 32, Issue 16, 15 August 2016, Pages 2473-2480, \n https://doi.org/10.1093/bioinformatics/btw191", 
                     "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------", 
               "\n \n"))
    x <- x[sample_info[, sample_col]]
    sam_names <- flowCore::sampleNames(x)
    param_f <- BiocGenerics::unique(gsub(param, pattern = "-H|-A|-W", 
                                         replacement = ""))
    filter_param <- BiocGenerics::colnames(x)
    filter_param <- BiocGenerics::unique(gsub(filter_param, 
                                              pattern = "-H|-A|-W", replacement = ""))
    filter_param <- filter_param[!filter_param %in% param_f & 
                                   filter_param != TimeChannel]
    filter_param <- c(filter_param, "FSC", "SSC")
    add_measuredparam <- base::unique(gsub(".*-([A-Z])$", 
                                           "\\1", param))[1]
    x <- flowAI::flow_auto_qc(x, alphaFR = 0.01, folder_results = "QC_flowAI", 
                              fcs_highQ = "HighQ", output = 1, timeCh = TimeChannel, 
                              ChExcludeFM = paste0(param_f[param_f %in% c("FSC", 
                                                                          "SSC")], "-", add_measuredparam), ChExcludeFS = filter_param, 
                              second_fractionFR = timesplit)
    x <- x[, param]
    x@phenoData@data$name <- sam_names
    flowCore::sampleNames(x) <- sam_names
    for (i in 1:length(x)) {
      x[[i]]@parameters@data[, 3] <- as.numeric(x[[i]]@parameters@data[, 
                                                                       3])
      x[[i]]@parameters@data[, 4] <- as.numeric(x[[i]]@parameters@data[, 
                                                                       4])
      x[[i]]@parameters@data[, 5] <- as.numeric(x[[i]]@parameters@data[, 
                                                                       5])
    }
  }
  x <- x[as.character(sample_info[, sample_col])]
  Biobase::pData(x) <- base::cbind(Biobase::pData(x), sample_info[base::order(base::match(as.character(sample_info[, 
                                                                                                                   sample_col]), as.character(Biobase::pData(x)[, "name"]))), 
                                                                  target_label])
  base::colnames(Biobase::pData(x))[2] <- target_label
  x <- x[, param]
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  x <- Phenoflow::FCS_resample(x, sample = downsample)
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  for (i in 1:length(x)) {
    tmp <- data.table::data.table(label = Biobase::pData(x[i])[, 
                                                               target_label], flowCore::exprs(x[[i]]))
    if (i == 1) {
      full_data <- tmp
    }
    else full_data <- base::rbind(full_data, tmp)
  }
  full_data <- base::droplevels(full_data)
  trainIndex <- caret::createDataPartition(full_data$label, 
                                           p = p_train)
  train_data <- full_data[trainIndex$Resample1, ]
  test_data <- full_data[-trainIndex$Resample1, ]
  if (metric == "Accuracy") {
    fitControl <- caret::trainControl(method = "repeatedcv", 
                                      number = 10, repeats = 3)
    metric <- "Accuracy"
    mtry <- round(base::sqrt(ncol(train_data)), 0)
    tunegrid <- base::expand.grid(.mtry = mtry)
    cat(date(), paste0("--- Training Random Forest classifier on ", 
                       100 * p_train, "% training set with options: \n", "\t-Performance metric: ", 
                       metric, "\n", "\t-Number of trees: ", 500, "\n", "\t-mtry: ", 
                       round(mtry, 2), "\n", "\t-method: ", fitControl$number, 
                       "x ", fitControl$method, "\n", "\t-repeats: ", fitControl$repeats, 
                       "\n"))
    cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
    RF_train <- caret::train(label ~ ., data = train_data, method = method, 
                             metric = metric, tuneGrid = tunegrid, trControl = fitControl, 
                             ntree = 500)
    print(RF_train)
    cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
    performance_metrics <- data.frame(metric = 1, n_cells = c(table(test_data$label)), 
                                      label = levels(as.factor(test_data$label)))
    for (n_label in 1:length(unique(test_data$label))) {
      tmp <- test_data[test_data$label == unique(test_data$label)[n_label], 
      ]
      tmp_pred <- stats::predict(RF_train, newdata = tmp)
      index <- performance_metrics$label == unique(test_data$label)[n_label]
      performance_metrics$metric[index] <- round(100 * base::sum(tmp_pred == 
                                                                   tmp$label)/performance_metrics$n_cells[index], 2)
    }
    colnames(performance_metrics)[1] <- "Accuracy"
    cat(date(), paste0("--- Performance on ", 100 * (1 - p_train), 
                       "% test set: \n"))
    print(performance_metrics)
    cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
    RF_pred <- stats::predict(RF_train, newdata = test_data)
    results_list <- list()
    results_list[[1]] <- RF_train
    results_list[[2]] <- caret::confusionMatrix(data = RF_pred, 
                                                as.factor(test_data$label))
  }
  
  if (metric == "MCC") {
    MatthewsCC <- function(data, lev = NULL, model = NULL)
    {
      out <- c(mltools::mcc(preds = data$pred, actuals = data$obs))
      names(out) <- c("MCC")
      out
    }
    
    fitControl <- caret::trainControl(method = "repeatedcv", 
                                      number = 5, repeats = 3, summaryFunction = MatthewsCC)
    metric <- "MCC"
    mtry <- round(base::sqrt(ncol(train_data)), 0)
    tunegrid <- base::expand.grid(.mtry = mtry)
    cat(date(), paste0("--- Training Random Forest classifier on ", 
                       100 * p_train, "% training set with options: \n", "\t-Performance metric: ", 
                       metric, "\n", "\t-Number of trees: ", 500, "\n", "\t-mtry: ", 
                       round(mtry, 2), "\n", "\t-method: ", fitControl$number, 
                       "x ", fitControl$method, "\n", "\t-repeats: ", fitControl$repeats, 
                       "\n"))
    cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
    RF_train <- caret::train(label ~ ., data = train_data, method = method, 
                             metric = metric, tuneGrid = tunegrid, trControl = fitControl, 
                             ntree = 500)
    print(RF_train)
    cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
    performance_metrics <- data.frame(metric = 1, n_cells = c(table(test_data$label)), 
                                      label = levels(as.factor(test_data$label)))
    for (n_label in 1:length(unique(test_data$label))) {
      tmp <- test_data[test_data$label == unique(test_data$label)[n_label], 
      ]
      tmp_pred <- stats::predict(RF_train, newdata = tmp)
      index <- performance_metrics$label == unique(test_data$label)[n_label]
      performance_metrics$metric[index] <- round(100 * base::sum(tmp_pred == 
                                                                   tmp$label)/performance_metrics$n_cells[index], 2)
    }
    colnames(performance_metrics)[1] <- "Accuracy"
    cat(date(), paste0("--- Performance on ", 100 * (1 - p_train), 
                       "% test set: \n"))
    print(performance_metrics)
    cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
    RF_pred <- stats::predict(RF_train, newdata = test_data)
    results_list <- list()
    results_list[[1]] <- RF_train
    results_list[[2]] <- caret::confusionMatrix(data = RF_pred, 
                                                as.factor(test_data$label))
  }
  
  if (plot_fig == TRUE) {
    mytable <- base::round(data.frame(performance = results_list[[2]]$overall), 
                           2)
    p_conf <- ggplot2::ggplot(data.frame(results_list[[2]]$table), 
                              ggplot2::aes(x = Prediction, y = Reference, fill = 100 * 
                                             Freq/sum(Freq))) + ggplot2::geom_raster() + 
      ggplot2::geom_text(ggplot2::aes(label = round(100 * 
                                                      Freq/sum(Freq), 0)), size = 6) + ggplot2::scale_fill_distiller(name = "% of total cells\n classified\n") + 
      ggplot2::theme_bw() + ggplot2::scale_x_discrete(position = "top") + 
      ggplot2::theme(axis.title = ggplot2::element_text(size = 16), 
                     strip.text.x = ggplot2::element_text(size = 14), 
                     legend.title = ggplot2::element_text(size = 14), 
                     legend.text = ggplot2::element_text(size = 14), 
                     axis.text.y = ggplot2::element_text(size = 13), 
                     axis.text.x = ggplot2::element_text(size = 13, 
                                                         angle = 55, hjust = 0), title = ggplot2::element_text(size = 20), 
                     plot.margin = ggplot2::unit(c(1.1, 1.1, 1.1, 
                                                   1.1), "cm"), panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(), 
                     panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "transparent", 
                                                                                                       colour = NA), plot.background = ggplot2::element_rect(fill = "transparent", 
                                                                                                                                                             colour = NA))
    p_conf_table <- ggplot2::ggplot(data.frame(results_list[[2]]$table), 
                                    ggplot2::aes(x = Prediction, y = Reference, fill = 100 * 
                                                   Freq/sum(Freq))) + ggplot2::theme(axis.title = ggplot2::element_text(size = 16), 
                                                                                     strip.text.x = ggplot2::element_text(size = 14), 
                                                                                     legend.title = ggplot2::element_text(size = 14), 
                                                                                     legend.text = ggplot2::element_text(size = 14), 
                                                                                     axis.text.y = ggplot2::element_text(size = 13), 
                                                                                     axis.text.x = ggplot2::element_text(size = 13, angle = 55, 
                                                                                                                         hjust = 0), title = ggplot2::element_text(size = 20), 
                                                                                     plot.margin = ggplot2::unit(c(1.1, 1.1, 1.1, 1.1), 
                                                                                                                 "cm"), panel.grid.major = ggplot2::element_blank(), 
                                                                                     panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_blank(), 
                                                                                     panel.background = ggplot2::element_rect(fill = "transparent", 
                                                                                                                              colour = NA), plot.background = ggplot2::element_rect(fill = "transparent", 
                                                                                                                                                                                    colour = NA)) + ggplot2::ylab("") + ggplot2::xlab("") + 
      ggplot2::annotation_custom(gridExtra::tableGrob(mytable))
    print(cowplot::plot_grid(p_conf, p_conf_table, align = "h", 
                             ncol = 2, rel_widths = c(1/2, 1/5)))
  }
  return(results_list)
}