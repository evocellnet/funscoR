#' Training the phosphosite functional score
#'
#' Function to train a phosphosite functional score based on a phosphoproteome and a provided gold standard.
#'
#' @param preprocessed_phos data frame of preprocessed features
#' @param residue residue set to train. Options are "ST" or "Y"
#' @param gold data frame containing a gold standard set of known functional phosphosites. D
#' ata frame should contain at least 2 columns: "acc" and "position".
#' @param parameters list of parameters with training settings
#' @param seed seed to reproduce the training
#' @param ncores number of cores to use in parallelized training. Defaults to 1.
#'
#' @return object with the trained model
#'
#' @import caret dplyr doParallel gbm e1071 pROC parallel foreach
#' @export
train_funscore <- function(preprocessed_phos, residue, gold, parameters = NULL, seed = NULL, ncores = 1) {
    cpu <- list()
    if (is.null(parameters) & residue == "ST") {
        ## Constants
        param <- list()
        param$outfoldtimes <- 5
        param$outfoldk <- 3
        param$alogrithm <- "gbm"
        param$n_trees <- 500
        param$interaction.depth <- 9
        param$shrinkage <- 0.0405
        param$n_minobsinnode <- 10
    } else if (is.null(parameters) & residue == "Y") {
        ## Constants
        param <- list()
        param$outfoldtimes <- 5
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
        mycluster <- parallel::makePSOCKcluster(ncores)
        doParallel::registerDoParallel(mycluster)
        cpu$parworker <- foreach::getDoParWorkers()
        cpu$parname <- foreach::getDoParName()
        cpu$parversion <- foreach::getDoParVersion()
    }
    feat <- preprocessed_phos
    goldreg <- gold
    ## Preparing features
    feat <- feat[, !names(feat) == "residue"]
    funclasses <- c("nonfunctional", "functional")
    feat$is_regulatory <- funclasses[unclass(as.factor(row.names(feat) %in% paste(goldreg$acc, goldreg$position, sep = "_")))]
    feat$is_regulatory <- as.factor(feat$is_regulatory)

    folds <- lapply(1:param$outfoldtimes, function(x)
        createFolds(feat$is_regulatory, k = param$outfoldk, returnTrain = TRUE))

    fit_control <- trainControl(method = "none",
                                classProbs = TRUE,
                                allowParallel = TRUE)

    grid_control <- data.frame(interaction.depth = param$interaction.depth,
                               n.trees = param$n_trees,
                               shrinkage = param$shrinkage,
                               n.minobsinnode = param$n_minobsinnode)
    time <- system.time({
        foldmodels <- parallel::mclapply(folds, function(iteration) {
            parallel::mclapply(iteration, function(x) {
                train_data <- feat[x, ]
                model <- train(is_regulatory ~ .,
                               data = train_data,
                               method = param$alogrithm,
                               metric = "ROC",
                               trControl = fit_control,
                               keep.data=FALSE,
                               tuneGrid = grid_control,
                               verbose = FALSE)
                return(model)
            })
        })
    })
    if (ncores > 1) {
        parallel::stopCluster(mycluster)
    }
    training <- list()
    training$folds <- folds
    training$models <- foldmodels
    training$variables <- list(parameters = param, cpu = cpu, time = time)

    return(training)
}

#' Predict functional score using trained model
#'
#' Apply model to predict functional score to a data frame of preprocessed features.
#'
#' @param preprocessed_phos preprocessed data frame of features
#' @param model trained model
#' @param ncores number of cores for parallelized computation. Defauls to 1.
#'
#' @return A data frame of scored phosphosites
#' @import caret dplyr doParallel gbm e1071 pROC parallel
#' @export
predict_funscore <- function(preprocessed_phos, model, ncores = 1) {
    if (ncores > 1) {
        mycluster <- parallel::makePSOCKcluster(ncores)
        doParallel::registerDoParallel(mycluster)
    }

    foldmodels <- model$models
    feat <- preprocessed_phos
    folds <- model$folds

    predictions <- parallel::mclapply(seq_len(length(foldmodels)), function(y) {
        do.call("rbind", parallel::mclapply(seq_len(length(foldmodels[[y]])), function(x) {
            model <- foldmodels[[y]][[x]]
            test_data <- feat[-folds[[y]][[x]], ]
            this_predictions <- predict(model, newdata = test_data, type = "prob")
            pred_model <- data.frame(sites = row.names(test_data),
                                     probabilities = as.vector(this_predictions[, "functional"]),
                                     stringsAsFactors=FALSE)
            return(pred_model)
        }))
    })
    if (ncores > 1) {
        parallel::stopCluster(mycluster)
    }
    ## merge all
    predictions <- predictions %>%
        Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "sites"), .)
    median_predictions <- data.frame(sites = predictions[, 1],
                                     probabilities = apply(predictions[, -1], 1, median))
    return(median_predictions)
}
