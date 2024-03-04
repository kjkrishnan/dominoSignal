#' check inputs to create_domino function
#' 
#' Issues a stop and error message if any of the data inputs to create_domino are incorrectly formatted or conflict with one another
#' 
#' @param rl_map Data frame where each row describes a receptor-ligand interaction with required columns gene_A & gene_B including the gene names for the receptor and ligand and type_A & type_B annotating if genes A and B are a ligand (L) or receptor (R)
#' @param features Either a path to a csv containing cell level features of interest (ie. the auc matrix from pySCENIC) or named matrix with cells as columns and features as rows.
#' @param ser Seurat object containing scaled RNA expression data in the RNA assay slot and cluster identity. Either a ser object OR z_scores and clusters must be provided. If ser is present z_scores and clusters will be ignored.
#' @param counts Counts matrix for the data. If a Seurat object is provided this will be ignored. This is only used to threshold receptors on dropout.
#' @param z_scores Matrix containing z-scored expression data for all cells with cells as columns and features as rows. Either z_scores and clusters must be provided OR a ser object. If ser is present z_scores and clusters will be ignored.
#' @param clusters Named factor containing cell cluster with names as cells. Either clusters and z_scores OR ser must be provided. If ser is present z_scores and clusters will be ignored.
#' @param tf_targets A list where names are transcription factors and the stored values are character vectors of genes in the transcription factor's regulon.
#' @return executes a function stop and error if checks do not pass or provides no output
#' @keywords internal
#' 

check_create_domino <- function(rl_map, features, ser, counts, clusters, tf_targets) {
  # format of rl_map
  stopifnot(
    `rl_map must be a data.frame with column names gene_A, gene_B, type_A, and type_B` = 
      (is(rl_map,"data.frame") & 
      c("gene_A", "gene_B", "type_A", "type_B") %in% colnames(rl_map))
    )
  # matrix or data.frame format of features
  stopifnot(
    `features must be either a file path or a named matrix with cells as columns and features as rows` = 
      (
        (is(features,"character") & length(features) == 1) | 
          (is(features, "matrix") & !is.null(rownames(features)) & !is.null(colnames(features))) |
          (is(features, "data.frame") & !is.null(rownames(features)) & !is.null(colnames(features)))
      )
  )
  # tf_targets is a list
  if(!is.null(tf_targets)) {
    stopifnot(
      `tf_targets should be a list containing character vectors of genes within a transcription factor's regulon named by the transcription factor` = 
        is(tf_targets, "list")
    )
  }
  # Either a Seurat object or expression matrices are provided as arguments
  stopifnot(
    `Either a Seurat object OR counts, z scores, and clusters must be provided` =
      (!is.null(ser) & (!is.null(counts) | !is.null(z_scores) | !is.null(clusters)))
  ) 
  # checks on expression matrix inputs if a Seurat object is not provided
  if(is.null(ser)) {
    # clusters provided are a named factor
    stopifnot(
      `clusters should be a named factor where value names are cell barcodes` = 
        (is(clusters, "factor") & !is.null(names(clusters)))
    )
    # counts matrix has rownames and colnames
    stopifnot(
      `counts should have feature rownames and sample colnames` = 
        (!is.null(rownames(counts)) & !is.null(colnames(counts)))
    )
    # z_scores matrix has rownames and colnames
    stopifnot(
      `z_scores should have feature rownames and sample colnames` = 
        (!is.null(z_scores) & !is.null(rownames(z_scores)) & !is.null(colnames(z_scores)))
    )
  }
}
