# Functions intended to aid in analysis of multiple objects
# (such as a data set with all cell types and another data set
# at higher resolution of a single cell type)
# TODO Write tests for these functions (will require test objects)
require(Matrix)
require(tidyverse)

#' Get incoming ligands to domino object
#'
#' This function returns a list of ligands that are incoming to a domino object.
#' 
#' @param dom A domino object
#' @param clusters A vector of cluster names to restrict the search to
#' @return A list with two elements: lig_names and complex_names
#' @keywords internal
get_incoming_ligands <- function(dom) {
    all_lig <- unlist(dom@linkages$clust_incoming_lig)
    all_lig <- unique(all_lig)
    all_lig <- all_lig[!all_lig == ""]
    all_lig_names_resolved <- resolve_names(dom, all_lig)
    
    if (length(dom@linkages$complexes) > 0) {
        all_lig_complexes_resolved_list <- resolve_complexes(dom, all_lig_names_resolved)
        all_lig_names_resolved <- unlist(all_lig_complexes_resolved_list)
    }
    
    return(list("lig_names" = unique(all_lig_names_resolved), "complex_names" = all_lig_complexes_resolved_list))
    }

#' Get network information from sending domino object to receiving domino object
#' 
#' Given a sending domino object and desired clusters as well as a receiving domino object and desired clusters,
#' this function returns a data frame of ligands expressed in the sending domino object by cluster
#' and links that information to matched receptors and transcription factors in the receiving domino object
#' 
#' @param send_dom A domino object (expressing ligands)
#' @param rec_dom A domino object (expressing receptors to receive ligands)
#' @param send_clusters A vector of cluster names to restrict the search to for the sending object
#' @param rec_clusters A vector of cluster names to restrict the search to for the receiving object
#' @param exp_type A character string specifying the type of expression data to use (either "counts" or "z_scores")
#' @return A data frame with columns for ligand, receptor, transcription factor, ligand expression,
#' receptor expression, and transcription factor AUC
#' @keywords internal
two_dom_network <- function(send_dom, rec_dom, send_clusters = NULL, rec_clusters = NULL, exp_type = c("counts", "z_scores")) {
    all_lig_names_resolved_lists <- get_incoming_ligands(rec_dom)
    all_lig_names_resolved <- all_lig_names_resolved_lists$lig_names
    all_lig_complexes_resolved <- all_lig_names_resolved_lists$complex_names
    
    if (exp_type == "counts") {
        lig_genes <- intersect(all_lig_names_resolved, rownames(send_dom@counts))
    } else if (exp_type == "z_scores") {
        lig_genes <- intersect(all_lig_names_resolved, rownames(send_dom@z_scores))
    }

    if (is.null(send_clusters)) {
        send_clusters <- levels(send_dom@clusters)
    }
    
    cl_ligands <- get_ligand_expression(send_dom, send_clusters, lig_genes, all_lig_complexes_resolved, exp_type)
    
    cl_ligands_sub <- reshape2::melt(cl_ligands)
    colnames(cl_ligands_sub) <- c("ligand", "cluster", "mean_counts")
    
    if (is.null(rec_clusters)) {
        rec_clusters <- levels(rec_dom@clusters)
    }
    
    df <- get_signaling_info(rec_dom, rec_clusters, cl_ligands_sub, exp_type)
    
    return(df)
}

#' Get network information between two domino objects
#' 
#' Given two domino objects and desired clusters, this function returns a data frame of ligands
#' expressed in each object that are matched to receptors and transcription factors in the other object
#' 
#' @param dom1 A domino object
#' @param dom2 A domino object
#' @param dom1_clust A vector of cluster names to restrict the search to for dom1
#' @param dom2_clust A vector of cluster names to restrict the search to for dom2
#' @param exp_type A character string specifying the type of expression data to use (either "counts" or "z_scores")
#' @return A data frame with columns for ligand, receptor, transcription factor, ligand expression,
#' receptor expression, and transcription factor AUC
#' @export
#' @examples
#' network_between_doms(pbmc_dom_built_tiny, pbmc_dom_built_tiny1, exp_type = "counts")
network_between_doms <- function(dom1, dom2, dom1_clust = NULL, dom2_clust = NULL, exp_type = c("counts", "z_scores")) {
    dom1_to_dom2 = two_dom_network(dom1, dom2, dom1_clust, dom2_clust, exp_type)
    dom2_to_dom1 = two_dom_network(dom2, dom1, dom2_clust, dom1_clust, exp_type)
    df = rbind(dom1_to_dom2, dom2_to_dom1)
    return(df)
}
