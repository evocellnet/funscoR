
annotate_sites <- function(phosphoproteome) {
    feature_list <- grep("feature", data(package = "funscoR")$results[, "Item"], value = TRUE)
    for (feat in feature_list) {
        df <- get(feat)
        df <- suppressWarnings(df[, -"window"])
        bycols <- intersect(names(phosphoproteome), names(dt))
        phosphoproteome <- merge(phosphoproteome, df, all.x = TRUE, by = bycols)
    }
    return(phosphoproteome)
}


preprocessFeatures <- function(annotated_phos, residue, features_to_exclude = c("acc", "position", "mq_siteid", "PWM_nkinTop01", 
    "PWM_nkinTop02", "sift_mean_score", "sift_acid_score", "sift_min_score", "quant_top10", "spectralcounts", "biological_samples")) {
    phos <- annotated_phos
    ## Exclude correlated and unnecessary features
    features_to_exclude <- c("acc", "position", "mq_siteid", "PWM_nkinTop01", "PWM_nkinTop02", "sift_mean_score", "sift_acid_score", 
        "sift_min_score", "quant_top10", "spectralcounts", "biological_samples")
    phos <- phos[, !names(phos) %in% toExcludeFeatures]
    
    if (residue == "ST") {
        phos <- phos[phos$residue %in% c("S", "T"), ]
    } else {
        phos <- phos[phos$residue == "Y", ]
    }
    
    ## Complete variables
    phos <- within(phos, {
        PEP[PEP == 0] <- min(PEP[PEP != 0])
        PEP <- -log10(PEP)
        is_disopred[is.na(is_disopred)] <- TRU
        disopred_score[is.na(disopred_score)] <- median(disopred_score, na.rm = TRUE)
        isProteinDomain <- domain_uniprot != "" & !is.na(domain_uniprot)
        isProteinKinaseDomain <- grepl("inase", domain_uniprot)
        isUniprotRegion <- region_uniprot != "" & !is.na(region_uniprot)
        isUniprotCompBias <- comp_bias_uniprot != "" & !is.na(comp_bias_uniprot)
        isUniprotRepeat <- repeat_uniprot != "" & !is.na(repeat_uniprot)
        isUniprotZnFinger <- zn_finger_uniprot != "" & !is.na(zn_finger_uniprot)
        isCytoplasmic <- topology_uniprot == "Cytoplasmic" & !is.na(topology_uniprot)
        isMotif <- motif_uniprot != "" & !is.na(motif_uniprot)
        isELMLinearMotif <- elm_CLV == TRUE | elm_DEG == TRUE | elm_DOC == TRUE | elm_LIG == TRUE | elm_MOD == TRUE | 
            elm_TRG == TRUE | isELMkinaseMotif == TRUE
        isELMLinearMotif[is.na(isELMLinearMotif)] <- FALSE
        isEV_ala_prediction_epistatic5 <- EV_ala_prediction_epistatic < -5 & !is.na(EV_ala_prediction_epistatic)
        exp3d_acid_ddG_effect <- factor(exp3d_acid_ddG_effect, levels = c(levels(exp3d_acid_ddG_effect), "unknown"))
        exp3d_ala_ddG_effect <- factor(exp3d_ala_ddG_effect, levels = c(levels(exp3d_ala_ddG_effect), "unknown"))
        exp3d_acid_ddG_effect[is.na(exp3d_acid_ddG_effect)] <- "unknown"
        exp3d_ala_ddG_effect[is.na(exp3d_ala_ddG_effect)] <- "unknown"
        isHotspot[is.na(isHotspot)] <- FALSE
        isHotspot[is.na(isHotspot)] <- FALSE
        log10_hotspot_pval_min[is.na(log10_hotspot_pval_min)] <- 0
        isInterface <- isInterface == TRUE & !is.na(isInterface)
        adj_ptms_w21[is.na(adj_ptms_w21)] <- 0
        netpho_max_all[is.na(netpho_max_all)] <- median(netpho_max_all, na.rm = TRUE)
        netpho_max_KIN[is.na(netpho_max_KIN)] <- median(netpho_max_KIN, na.rm = TRUE)
        paxdb_abundance_log10[is.na(paxdb_abundance_log10)] <- median(paxdb_abundance_log10, na.rm = TRUE)
        w0_mya[is.na(w0_mya)] <- 0
        w3_mya[is.na(w3_mya)] <- 0
        isKinaseCoreg <- !is.na(ptmdb_maxcoreg_kinase)
        medianPubmedCounts <- median(pubmed_counts, na.rm = TRUE)
        pubmed_counts[is.na(pubmed_counts)] <- medianPubmedCounts
        quant_top1[is.na(quant_top1)] <- 0
        quant_top5[is.na(quant_top5)] <- 0
        PWM_max_mss[is.na(PWM_max_mss)] <- median(PWM_max_mss, na.rm = TRUE)
        SSpro <- factor(SSpro, levels = c(levels(SSpro), "-"))
        SSpro8 <- factor(SSpro8, levels = c(levels(SSpro8), "-"))
        ACCpro[is.na(ACCpro)] <- "-"
        SSpro[is.na(SSpro)] <- "-"
        SSpro8[is.na(SSpro8)] <- "-"
        sift_ala_score[is.na(sift_ala_score)] <- median(sift_ala_score, na.rm = TRUE)
    })
    
    ## only complete columns or sift are included in preprocessing
    columnsToProcess <- apply(phos, 2, function(x) sum(is.na(x))) == 0 | grepl("sift", names(phos))
    phos <- phos[, columnsToProcess]
    
    preProcValues <- preProcess(phos, method = c("center", "scale", "YeoJohnson", "nzv", "knnImpute"))
    
    transformed <- predict(preProcValues, newdata = phos)
    return(transformed)
}
