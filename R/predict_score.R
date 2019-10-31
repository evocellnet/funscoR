predict_score <- function(annotated_phos, residue, gold, parameters = NULL, seed = NULL, ncores = 1) {
    if (is.null(parameters) & residue == "ST") {
        ## Constants
        param <- list()
        param$outfoldtimes <- 30
        param$outfoldk <- 3
        param$alogrithm <- "gbm"
        param$n.trees <- 500
        param$interaction.depth <- 9
        param$shrinkage <- 0.0405
        param$n.minobsinnode <- 10
    } else if (is.null(parameters) & residue == "Y") {
        ## Constants
        param <- list()
        param$outfoldtimes <- 30
        param$outfoldk <- 3
        param$alogrithm <- "gbm"
        param$n.trees <- 500
        param$interaction.depth <- 9
        param$shrinkage <- 0.0105
        param$n.minobsinnode <- 10
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (ncores > 1) {
        registerDoParallel(cores = ncores)
        cpu$parworker <- getDoParWorkers()
        cpu$parname <- getDoParName()
        cpu$parversion <- getDoParVersion()
    }
    feat <- annotated_phos
    gold <- goldreg
    
    ## Preparing features
    feat <- feat[, !names(feat) == "residue"]
    funclasses <- c("nonfunctional", "functional")
    feat$isRegulatory <- funclasses[unclass(as.factor(row.names(feat) %in% paste(goldreg$V1, goldreg$V2, sep = "_")))]
    feat$isRegulatory <- as.factor(feat$isRegulatory)
    
    folds <- lapply(1:param$outfoldtimes, function(x) createFolds(feat$isRegulatory, k = param$outfoldk, returnTrain = TRUE))
    
    fitControl <- trainControl(method = "none", classProbs = TRUE)
    
    gridControl <- data.frame(interaction.depth = param$interaction.depth, n.trees = param$n.trees, shrinkage = param$shrinkage, 
        n.minobsinnode = param$n.minobsinnode)
    time <- system.time({
        foldmodels <- lapply(folds, function(iteration) {
            lapply(iteration, function(x) {
                trainData <- feat[x, ]
                model <- train(isRegulatory ~ ., data = trainData, method = param$alogrithm, metric = "ROC", trControl = fitControl, 
                  tuneGrid = gridControl, verbose = FALSE)
                return(model)
            })
        })
    })
    
    predictions <- lapply(1:length(foldmodels), function(y) {
        do.call("rbind", lapply(1:length(foldmodels[[y]]), function(x) {
            model <- foldmodels[[y]][[x]]
            testData <- feat[-folds[[y]][[x]], ]
            pred.model <- data.frame(sites = row.names(testData), probabilities = as.vector(predict(model, newdata = testData, 
                type = "prob")[, "functional"]))
            return(pred.model)
        }))
    })
    
    ## merge all
    predictions <- predictions %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "sites"), .)
    
    median_predictions <- data.frame(sites = predictions[, 1], probabilities = apply(predictions[, -1], 1, median))
    
    variables <- list(parameters = param, models = foldmodels, cpu = cpu, time = time)
    
    out <- list()
    out$result <- median_predictions
    out$variables <- variables
    
    return(out)
}
