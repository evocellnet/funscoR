
#' Annotate phosphoproteome with functional features
#'
#' This function takes annotation files contained within the package to annotate a
#' phosphoproteome with extra functional information.
#'
#' @param phosphoproteome A data frame with the columns `acc`, `position` and `residue`
#' referring to a list of phosphorylation sites
#'
#' @return A data frame expanding the input data frame with a list of functional features.
#'
#' @importFrom utils data
#' @export
annotate_sites <- function(phosphoproteome) {
    feature_list <- grep("feature", data(package = "funscoR")$results[, "Item"], value = TRUE)
    for (feat in feature_list) {
        df <- get(feat)
        df <- df[, !names(df) == "window"]
        bycols <- intersect(names(phosphoproteome), names(df))
        phosphoproteome <- merge(phosphoproteome, df, all.x = TRUE, by = bycols)
    }
    phosphoproteome <- unique(phosphoproteome)
    row.names(phosphoproteome) <- paste(phosphoproteome$acc, phosphoproteome$position, sep = "_")
    return(phosphoproteome)
}

#' Preprocess features
#'
#' @param annotated_phos A data frame containing an annotated phosphoproteome.
#' @param residue The residue set to preprocess. Options are "ST" or "Y"
#' @param methods Methods used to preprocessed the features using the `preProcess`
#' function included in the `caret` package.
#' @param features_to_exclude A list of features to exclude from the preprocessing step.
#' Defaults to "acc" and "pos" to be excluded.
#'
#' @return A data frame with the preprocessed features ready for training
#' @import caret dplyr tibble
#' @importFrom tidyr replace_na
#' @importFrom stats predict
#' @importFrom stats median
#' @export
preprocess_features <- function(annotated_phos, residue, methods = c("center", "scale", "YeoJohnson", "nzv", "knnImpute"), features_to_exclude = c("acc", "position")) {
    phos <- annotated_phos
    ## Exclude correlated and unnecessary features
    phos <- phos[, !names(phos) %in% features_to_exclude]

    if (residue == "ST") {
        phos <- phos[phos$residue %in% c("S", "T"), ]
    } else {
        phos <- phos[phos$residue == "Y", ]
    }

    ## Complete variables
    phos <- phos %>%
        tibble::rownames_to_column(var = "rowname") %>%
        mutate(PEP = replace(PEP, PEP == 0, min(PEP[PEP != 0])),
               PEP = -log10(PEP),
               isProteinDomain = domain_uniprot != "" & !is.na(domain_uniprot),
               isProteinKinaseDomain = grepl("inase", domain_uniprot),
               isUniprotRegion = region_uniprot != "" & !is.na(region_uniprot),
               isUniprotCompBias = comp_bias_uniprot != "" & !is.na(comp_bias_uniprot),
               isUniprotRepeat = repeat_uniprot != "" & !is.na(repeat_uniprot),
               isUniprotZnFinger = zn_finger_uniprot != "" & !is.na(zn_finger_uniprot),
               isCytoplasmic = topology_uniprot == "Cytoplasmic" & !is.na(topology_uniprot),
               isMotif = motif_uniprot != "" & !is.na(motif_uniprot),
               isELMLinearMotif = elm_CLV | elm_DEG | elm_DOC | elm_LIG | elm_MOD | elm_TRG | isELMkinaseMotif,
               isEV_ala_prediction_epistatic5 = EV_ala_prediction_epistatic < -5 & !is.na(EV_ala_prediction_epistatic),
               isInterface = isInterface == TRUE & !is.na(isInterface),
               isKinaseCoreg = !is.na(ptmdb_maxcoreg_kinase),
               SSpro = as.factor(replace_na(SSpro, "-")),
               SSpro8 = as.factor(replace_na(SSpro8, "-")),
               ACCpro = as.factor(replace_na(ACCpro, "-")),
               exp3d_ala_ddG_effect = as.factor(replace_na(exp3d_ala_ddG_effect, "unknown")),
               exp3d_acid_ddG_effect = as.factor(replace_na(exp3d_acid_ddG_effect, "unknown")),
               isHotspot = replace_na(isHotspot, FALSE),
               isELMLinearMotif = replace_na(isELMLinearMotif, FALSE),
               is_disopred = replace_na(is_disopred, TRUE),
               quant_top1 = replace_na(quant_top1, 0),
               quant_top5 = replace_na(quant_top5, 0),
               w0_mya = replace_na(w0_mya, 0),
               w3_mya = replace_na(w3_mya, 0),
               adj_ptms_w21 = replace_na(adj_ptms_w21, 0),
               log10_hotspot_pval_min = replace_na(log10_hotspot_pval_min, 0),
               disopred_score = replace_na(disopred_score, median(disopred_score, na.rm = TRUE)),
               netpho_max_all = replace_na(netpho_max_all, median(netpho_max_all, na.rm = TRUE)),
               netpho_max_KIN = replace_na(netpho_max_KIN, median(netpho_max_KIN, na.rm = TRUE)),
               paxdb_abundance_log10 = replace_na(paxdb_abundance_log10, median(paxdb_abundance_log10, na.rm = TRUE)),
               sift_ala_score = replace_na(sift_ala_score, median(sift_ala_score, na.rm = TRUE)),
               pubmed_counts = replace_na(pubmed_counts, median(pubmed_counts, na.rm = TRUE)),
               PWM_max_mss = replace_na(PWM_max_mss, median(PWM_max_mss, na.rm = TRUE)))
    row.names(phos) <- phos$rowname
    phos <- phos %>% select(-rowname)
    ## only complete columns or sift are included in preprocessing
    columns_to_process <- apply(phos, 2, function(x) sum(is.na(x))) == 0 | grepl("sift", names(phos))
    phos <- phos[, columns_to_process]

    pre_proc_values <- preProcess(phos, method = methods)

    transformed <- predict(pre_proc_values, newdata = phos)
    return(transformed)
}
