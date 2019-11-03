#' Disorder predictions from disopred
#'
#' A dataset containing disorder predictions for every S, T or Y.
#'
#' @format Data frame with 1.87M rows
#' \describe{
#'   \item{acc}{uniprot accession for protein}
#'   \item{position}{position of the acceptor residue}
#'   \item{residue}{residue of the acceptor residue}
#'   \item{is_disopred}{is the residue in a disordered region}
#'   \item{disopred_score}{score for disopred prediction}
#' }
#' @source Ochoa et al. The functional landscape of the human phosphoproteome
"feature_disopred"

#' Disorder annotations in Disprot
#'
#' A dataset containing disorder annotations for S, T or Y.
#'
#' @format Data frame with 5,255 rows
#' \describe{
#'   \item{acc}{uniprot accession for protein}
#'   \item{position}{position of the acceptor residue}
#'   \item{residue}{residue of the acceptor residue}
#'   \item{disprot_disorder}{is annotated as a disordered region}
#' }
#' @source Ochoa et al. The functional landscape of the human phosphoproteome
"feature_disprot"

#' Disorder annotations in Disprot
#'
#' A dataset containing region annotations from Uniprot.
#'
#' @format Data frame with 5,255 rows
#' \describe{
#'   \item{acc}{uniprot accession for protein}
#'   \item{position}{position of the acceptor residue}
#'   \item{residue}{residue of the acceptor residue}
#'   \item{Coiled_coil_uniprot}{is annotated as a coiled coil in uniprot}
#'   \item{comp_bias_uniprot}{is annotated as a compositional bias region in uniprot}
#'   \item{domain_uniprot}{is annotated as a domain in uniprot}
#'   \item{motif_uniprot}{is annotated as a motif in uniprot}
#'   \item{region_uniprot}{is annotated as a region in uniprot}
#'   \item{repeat_uniprot}{is annotated as a repeat region in uniprot}
#'   \item{zn_finger_uniprot}{is annotated as a Zn finger region in uniprot}
#' }
#' @source Ochoa et al. The functional landscape of the human phosphoproteome
"feature_domains"

"feature_elm"
"feature_evmut"
"feature_exac"
"feature_foldx"
"feature_hotspots"
"feature_interfaces"
"feature_ms_pride"
"feature_neighPTMs"
"feature_netphorest"
"feature_paxdb"
"feature_proteinlength"
"feature_ptmdb_age"
"feature_ptmdb_coregulation"
"feature_ptmdb_counts"
"feature_ptmdb_regulation"
"feature_pwm_match"
"feature_scratch1Dfeatures"
"feature_sift_scores"
"feature_spectral_counts"
"feature_topology"
"feature_transitpeptides"

#' Reference human phosphoproteome
#'
#' List of human phosphorylation sites gathered from a big search of mass spectrometry
#' data contained in the PRIDE database.
#'
#' @format Data frame with 116,258 rows
#' \describe{
#'   \item{acc}{uniprot accession for protein}
#'   \item{position}{position of the acceptor residue}
#'   \item{residue}{residue of the acceptor residue}
#' }
#' @source Ochoa et al. The functional landscape of the human phosphoproteome
"phosphoproteome"

#' List of human phosphosites with known regulatory functions
#'
#' This list is the result of downloading and parsing the sites reported in the PhosphositePlus database (see license).
#'
#' @format Data frame with 6,249 rows
#' \describe{
#'   \item{acc}{uniprot accession for protein}
#'   \item{position}{position of the acceptor residue}
#' }
#' @source PhosphositePlus \url{https://www.phosphosite.org/}
"psp"
