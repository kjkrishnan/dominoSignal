# Write tests for functions and component functions turning network
# data into a data frame.

library(testthat)

# Start with resolve_names
test_that("resolve_names replaces L.name with gene if L.name and gene are different, and the gene name has no comma.", {
    expect_equal(resolve_names(pbmc_dom_built_tiny, "integrin_aMb2_complex"), "integrin_aMb2_complex")
    expect_equal(resolve_names(pbmc_dom_built_tiny, "IGF1"), "IGF1")
})

# Test resolve_complexes
test_that("resolve_complexes replaces complex names with their components.", {
    expect_equal(resolve_complexes(pbmc_dom_built_tiny, "integrin_aMb2_complex"), c("ITGB2", "ITGAM"))
    expect_equal(resolve_complexes(pbmc_dom_built_tiny, "IL7_receptor"), c("IL7R", "IL2RG"))
    expect_equal(resolve_complexes(pbmc_dom_built_tiny, "IGF1"), "IGF1")
})

# Test avg_exp_for_complexes
test_that("function returns list with average expression of component genes in complexes.", {
    exp_mat = pbmc_dom_built_tiny@counts
    complexes = pbmc_dom_built_tiny@linkages$complexes
    expect_equal(length(avg_exp_for_complexes(exp_mat, complexes), length(complexes)))
    expect_equal(length(avg_exp_for_complexes(exp_mat, complexes)[[1]]), ncol(exp_mat))
    expect_equal(avg_exp_for_complexes(exp_mat, complexes)[[1]][1], mean(exp_mat[c("ITGB2", "ITGAM"), ]))
})

# Test get_resolved_ligands
test_that("function returns list of resolved ligands.", {
    expect_equal(get_resolved_ligands(pbmc_dom_built_tiny), c("ITGB2", "ITGAM", "IGF1", "TNF", "IL7"))
})

# Test get_ligand_expression
test_that("function returns a matrix of expression", {
    expect_equal(get_ligand_expression(pbmc_dom_built_tiny, send_clusters = levels(pbmc_dom_built_tiny@clusters),
        lig_genes = "ITGB2", exp_type = "counts"), pbmc_dom_built_tiny@counts["ITGB2", ])
    expect_equal(get_ligand_expression(pbmc_dom_built_tiny, send_clusters = levels(pbmc_dom_built_tiny@clusters), 
        outgoing_cluster = "B_cell", lig_genes = c("IGF1", "IL7"), exp_type = "counts"),
        pbmc_dom_built_tiny@counts["IGF1", which(pbmc_dom_built_tiny@clusters == "B_cell")])
})

# Test get_signaling_info


# Test network_to_df