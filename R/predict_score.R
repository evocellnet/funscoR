train_funscore <- function(preprocessed_phos, residue, gold, parameters = NULL, seed = NULL, ncores = 1) {
    cpu <- list()
    if (is.null(parameters) & residue == "ST") {
        ## Constants
        param <- list()
        param$outfoldtimes <- 10
        param$outfoldk <- 3
        param$alogrithm <- "gbm"
        param$n_trees <- 500
        param$interaction.depth <- 9
        param$shrinkage <- 0.0405
        param$n_minobsinnode <- 10
    } else if (is.null(parameters) & residue == "Y") {
        ## Constants
        param <- list()
        param$outfoldtimes <- 10
        param$outfoldk <- 3
        param$alogrithm <- "gbm"
        param$n_trees <- 500
        param$interaction.depth <- 9
        param$shrinkage <- 0.0105
        param$n_minobsinnode <- 10
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (ncores > 1) {
        mycluster <- makePSOCKcluster(ncores)
        registerDoParallel(mycluster)
        cpu$parworker <- getDoParWorkers()
        cpu$parname <- getDoParName()
        cpu$parversion <- getDoParVersion()
    }
    feat <- preprocessed_phos
    goldreg <- gold
    ## Preparing features
    feat <- feat[, !names(feat) == "residue"]
    funclasses <- c("nonfunctional", "functional")
    feat$is_regulatory <- funclasses[unclass(as.factor(row.names(feat) %in% paste(goldreg$V1, goldreg$V2, sep = "_")))]
    feat$is_regulatory <- as.factor(feat$is_regulatory)

    folds <- lapply(1:param$outfoldtimes, function(x) createFolds(feat$is_regulatory, k = param$outfoldk, returnTrain = TRUE))

    fit_control <- trainControl(method = "none", classProbs = TRUE, allowParallel = TRUE)

    grid_control <- data.frame(interaction.depth = param$interaction.depth, n.trees = param$n_trees, shrinkage = param$shrinkage,
        n.minobsinnode = param$n_minobsinnode)
    time <- system.time({
        foldmodels <- mclapply(folds, function(iteration) {
            mclapply(iteration, function(x) {
                train_data <- feat[x, ]
                model <- train(is_regulatory ~ ., data = train_data, method = param$alogrithm, metric = "ROC", trControl = fit_control,
                  tuneGrid = grid_control, verbose = FALSE)
                return(model)
            })
        })
    })
    if (ncores > 1) {
        stopCluster(mycluster)
    }
    training <- list()
    training$folds <- folds
    training$models <- foldmodels
    training$variables <- list(parameters = param, cpu = cpu, time = time)

    return(training)
}

predict_funscore <- function(preprocessed_phos, training, ncores = 1) {
    if (ncores > 1) {
        mycluster <- makePSOCKcluster(ncores)
        registerDoParallel(mycluster)
    }

    foldmodels <- training$models
    feat <- preprocessed_phos
    folds <- training$folds

    predictions <- mclapply(seq_len(length(foldmodels)), function(y) {
        do.call("rbind", mclapply(seq_len(length(foldmodels[[y]])), function(x) {
            model <- foldmodels[[y]][[x]]
            test_data <- feat[-folds[[y]][[x]], ]
            this_predictions <- predict(model, newdata = test_data, type = "prob")
            pred_model <- data.frame(sites = row.names(test_data), probabilities = as.vector(this_predictions[, "functional"]))
            return(pred_model)
        }))
    })
    if (ncores > 1) {
        stopCluster(mycluster)
    }
    ## merge all
    predictions <- suppressWarnings(Reduce(predictions, function(dtf1, dtf2) full_join(dtf1, dtf2, by = "sites"), .))
    median_predictions <- data.frame(sites = predictions[, 1], probabilities = apply(predictions[, -1], 1, median))
    return(median_predictions)
}
