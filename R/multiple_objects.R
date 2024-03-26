# Functions intended to aid in analysis of multiple objects
# (such as a data set with all cell types and another data set
# at higher resolution of a single cell type)

# Functions structured such that ligands are coming from dom1
# and receptors and tf activation are from dom2

# Get all incoming ligands from dom2
get_incoming_ligands <- function(dom, clusters = NULL) {
    ligs = dom@linkages$clust_incoming_lig
    if (is.null(clusters)) {clusters = levels(dom@clusters)}
    ligs = ligs[clusters]
    # List name to receiving cluster, lig to ligand in data frame:
    cl_lig = lapply(1:length(ligs), function(i) {
        data.frame(rec_clust = names(ligs)[i], ligand = ligs[[i]])
    })
    cl_lig_df = purrr::list_rbind(cl_lig)
}

# Get ligand expression results for all ligands from dom1
# Use existing get_ligand_expression function with
# dom as dom1, send_clusters from dom1, lig_genes from get_incoming_ligands,
# complexes from dom1, and exp_type as zscore or counts
# Result is sending_cluster by lig_gene matrix
# Melt to sending_cluster, lig_gene, and mean_exp
# Match to receiving_cluster from output of get_incoming_ligands
add_receiving_clusters = function(lig_df, lig_exp) {
    combined_df = full_join(lig_df, lig_exp, by = "ligand", relationship = "many-to-many")
}