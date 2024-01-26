#' @import stats
#' @importFrom utils read.csv
#' @importFrom Matrix rowSums
#' @import methods
#'
NULL

#' Create a list of genes in regulons inferred by SCENIC
#'
#' Generates a list of transcription factors and the genes targeted by the transcription factor as part of their regulon inferred by pySCENIC
#'
#' @param regulons Dataframe or file path to the table of the output of the grn (gene regulatory network) function from pySCENIC
#' @return A list where names are transcription factors and the stored values are character vectors of genes in the inferred regulons
#' @export create_regulon_list_scenic
#' @examples
#' regulon_list_tiny <- create_regulon_list_scenic(regulons = domino2:::regulons_tiny)
#'
create_regulon_list_scenic <- function(regulons) {
  if (is(regulons, "character")) {
    regulons <- read.csv(regulons)
  }
  TFS <- unique(regulons[["TF"]])
  TF_targets <- lapply(TFS, function(tf) {
    regulons_small <- regulons[regulons[["TF"]] == tf, ]
    targets <- regulons_small[["TargetGenes"]]
    target_genes <- lapply(targets, function(x) {
      split_targs <- unlist(strsplit(x, ""))
      split_targs <- split_targs[seq(2, length(split_targs))]
      split_targs <- split_targs[seq(1, length(split_targs) - 1)]
      split_targs <- paste(split_targs, collapse = "")
      split_targs <- unlist(strsplit(split_targs, "['), (']"))
      split_targs_pre <- split_targs[split_targs != ""]
      split_targs_post <- split_targs_pre[seq(1, length(split_targs_pre), 2)]
      return(split_targs_post)
    })
    return(unique(unlist(target_genes)))
  })
  names(TF_targets) <- TFS
  return(TF_targets)
}

#' Create a domino object and prepare it for network construction
#'
#' This function reads in a receptor ligand signaling database, cell level
#' features of some kind (ie. output from pySCENIC), z-scored single cell data,
#' and cluster id for single cell data, calculates a correlation matrix between
#' receptors and other features (this is transcription factor module scores if
#' using pySCENIC), and finds features enriched by cluster. It will return a
#' domino object prepared for [build_domino()], which will calculate a signaling
#' network.
#'
#' @param rl_map Data frame where each row describes a receptor-ligand interaction with required columns gene_A & gene_B including the gene names for the receptor and ligand and type_A & type_B annotating if genes A and B are a ligand (L) or receptor (R)
#' @param features Either a path to a csv containing cell level features of interest (ie. the auc matrix from pySCENIC) or named matrix with cells as columns and features as rows.
#' @param ser Seurat object containing scaled RNA expression data in the RNA assay slot and cluster identity. Either a ser object OR z_scores and clusters must be provided. If ser is present z_scores and clusters will be ignored.
#' @param counts Counts matrix for the data. If a Seurat object is provided this will be ignored. This is only used to threshold receptors on dropout.
#' @param z_scores Matrix containing z-scored expression data for all cells with cells as columns and features as rows. Either z_scores and clusters must be provided OR a ser object. If ser is present z_scores and clusters will be ignored.
#' @param clusters Named factor containing cell cluster with names as cells. Either clusters and z_scores OR ser must be provided. If ser is present z_scores and clusters will be ignored.
#' @param use_clusters Boolean indicating whether to use the clusters from a Seurat object. If a Seurat object is not provided then this parameter is ignored.
#' @param tf_targets Optional. A list where names are transcription factors and the stored values are character vectors of genes in the transcription factor's regulon.
#' @param verbose Boolean indicating whether or not to print progress during computation.
#' @param use_complexes Boolean indicating whether you wish to use receptor/ligand complexes in the receptor ligand signaling database. If FALSE, receptor/ligand pairs where either functions as a protein complex will not be considered when constructing the signaling network.
#' @param rec_min_thresh Minimum expression level of receptors by cell. Default is 0.025 or 2.5 percent of all cells in the data set. This is important when calculating correlation to connect receptors to transcription activation. If this threshold is too low then correlation calculations will proceed with very few cells with non-zero expression.
#' @param remove_rec_dropout Whether to remove receptors with 0 expression counts when calculating correlations. This can reduce false positive correlation calculations when receptors have high dropout rates.
#' @param tf_selection_method Selection of which method to target transcription factors. If 'clusters' then differential expression for clusters will be calculated. If 'variable' then the most variable transcription factors will be selected. If 'all' then all transcription factors in the feature matrix will be used. Default is 'clusters'. Note that if you wish to use clusters for intercellular signaling downstream to MUST choose clusters.
#' @param tf_variance_quantile What proportion of variable features to take if using variance to threshold features. Default is 0.5. Higher numbers will keep more features. Ignored if tf_selection_method is not 'variable'
#' @return A domino object
#' @export create_domino
#' @examples 
#' pbmc_dom_tiny_all <- pbmc_dom_tiny <- create_domino(
#'  rl_map = domino2:::rl_map_tiny, features = domino2:::auc_tiny, 
#'  counts = domino2:::RNA_count_tiny, z_scores = domino2:::RNA_zscore_tiny,
#'  clusters = domino2:::clusters_tiny, tf_targets = domino2:::regulon_list_tiny, 
#'  use_clusters = FALSE, use_complexes = FALSE, 
#'  rec_min_thresh = 0.1, remove_rec_dropout = TRUE,
#'  tf_selection_method = "all")
#' 
#' pbmc_dom_tiny_clustered <- create_domino(
#'  rl_map = domino2:::rl_map_tiny, features = domino2:::auc_tiny, 
#'  counts = domino2:::RNA_count_tiny, z_scores = domino2:::RNA_zscore_tiny,
#'  clusters = domino2:::clusters_tiny, tf_targets = domino2:::regulon_list_tiny,
#'  use_clusters = TRUE, use_complexes = TRUE, remove_rec_dropout = FALSE)
#' 
create_domino <- function(
    rl_map, features, ser = NULL, counts = NULL, z_scores = NULL, clusters = NULL,
    use_clusters = TRUE, tf_targets = NULL, verbose = TRUE, use_complexes = TRUE, rec_min_thresh = 0.025,
    remove_rec_dropout = TRUE, tf_selection_method = "clusters", tf_variance_quantile = 0.5) {
  # Check inputs:
  stopifnot(`rl_map must be a data.frame with column names gene_A, gene_B, type_A, and type_B` = (is(
    rl_map,
    "data.frame"
  ) & c("gene_A", "gene_B", "type_A", "type_B") %in% colnames(rl_map)))
  stopifnot(`features must be either a file path or a named matrix with cells as columns and features as rows` = ((is(
    features,
    "character"
  ) & length(features) == 1) | (is(features, "matrix") & !is.null(rownames(features)) &
    !is.null(colnames(features))) | (is(features, "data.frame") & !is.null(rownames(features)) &
    !is.null(colnames(features)))))
  stopifnot(`Either a Seurat object OR counts, z scores, and clusters must be provided` = (is(ser, "Seurat") |
    (!is.null(counts) & !is.null(rownames(counts)) & !is.null(colnames(counts)) &
      is(z_scores, "matrix") & !is.null(rownames(z_scores)) & !is.null(colnames(z_scores)) & is(
      clusters,
      "factor"
    ) & !is.null(names(clusters)))))
  stopifnot(`rec_min_thresh must be a number between 0 and 1` = (is(rec_min_thresh, "numeric") &
    rec_min_thresh <= 1 & rec_min_thresh >= 0))
  # Create object
  dom <- domino()
  dom@misc[["create"]] <- TRUE
  dom@misc[["build"]] <- FALSE
  dom@misc[["build_vars"]] <- NULL
  if (!is.null(ser) & (!is.null(clusters) | !is.null(z_scores) | !is.null(counts))) {
    warning("Ser and z_score, clusters, or counts provided. Defaulting to ser.")
  }
  if (is.null(ser) & (is.null(clusters) | is.null(z_scores) | is.null(counts))) {
    stop("Either ser or clusters and z_scores must be provided")
  }
  if (!(tf_selection_method %in% c("all", "clusters", "variable"))) {
    stop("tf_selection_method must be one of all, clusters, or variable")
  }
  # Read in lr db info
  if (verbose) {
    message("Reading in and processing signaling database")
  }
  if ("database_name" %in% colnames(rl_map)) {
    dom@db_info <- rl_map
    if (verbose) {
      message(paste0("Database provided from source: ", unique(rl_map[["database_name"]])))
    }
  } else {
    dom@db_info <- rl_map
  }
  # check for receptors that match receptor complex syntax of comma seperated genes
  non_complex_index <- which(!grepl("\\,", rl_map[["gene_A"]]) & !grepl("\\,", rl_map[["gene_B"]]))
  # discard interactions including complexes if requested
  if (use_complexes == FALSE) {
    rl_map <- rl_map[non_complex_index, ]
  }
  # Get genes for receptors
  rl_reading <- NULL
  for (i in 1:nrow(rl_map)) {
    rl <- list()
    inter <- rl_map[i, ]
    p <- ifelse(inter[["type_A"]] == "R", "A", "B")
    q <- ifelse(p == "A", "B", "A")
    R.gene <- inter[[paste0("gene_", p)]]
    L.gene <- inter[[paste0("gene_", q)]]
    rl[["R.gene"]] <- R.gene
    rl[["L.gene"]] <- L.gene
    if (paste0("uniprot_", p) %in% names(inter)) {
      rl[["R.uniprot"]] <- inter[[paste0("uniprot_", p)]]
    }
    if (paste0("uniprot_", q) %in% names(inter)) {
      rl[["L.uniprot"]] <- inter[[paste0("uniprot_", q)]]
    }
    if (paste0("name_", p) %in% names(inter)) {
      rl[["R.name"]] <- inter[[paste0("name_", p)]]
    }
    if (paste0("name_", q) %in% names(inter)) {
      rl[["L.name"]] <- inter[[paste0("name_", q)]]
    }
    rl <- as.data.frame(rl)
    rl_reading <- rbind(rl_reading, rl)
  }
  # save a list of complexes and their components
  dom@linkages$complexes <- NULL
  if (use_complexes) {
    complex_list <- list()
    for (i in 1:nrow(rl_reading)) {
      inter <- rl_reading[i, ]
      if (grepl("\\,", inter[["L.gene"]])) {
        complex_list[[inter[["L.name"]]]] <- unlist(strsplit(inter[["L.gene"]], split = "\\,"))
      }
      if (grepl("\\,", inter[["R.gene"]])) {
        complex_list[[inter[["R.name"]]]] <- unlist(strsplit(inter[["R.gene"]], split = "\\,"))
      }
    }
    dom@linkages$complexes <- complex_list
  }
  rec_genes <- unique(unlist(strsplit(rl_reading[["R.gene"]], split = "\\,")))
  rec_names <- rl_reading[["R.name"]]
  lig_genes <- unique(unlist(strsplit(rl_reading[["L.gene"]], split = "\\,")))
  lig_names <- rl_reading[["L.name"]]
  # building RL linkages
  rec_lig_linkage <- list()
  for (rec in rec_names) {
    inter <- rl_reading[rl_reading[["R.name"]] == rec, ]
    ligs <- inter[["L.name"]]
    rec_lig_linkage[[rec]] <- ligs
  }
  dom@linkages[["rec_lig"]] <- rec_lig_linkage
  dom@misc[["rl_map"]] <- rl_reading
  # Get z-score and cluster info
  if (verbose) {
    message("Getting z_scores, clusters, and counts")
  }
  if (!is.null(ser)) {
    z_scores <- ser@assays$RNA@scale.data
    if (use_clusters) {
      clusters <- ser@active.ident
    }
    counts <- ser@assays$RNA@counts
  }
  dom@z_scores <- z_scores
  if (!is.null(clusters)) {
    dom@clusters <- clusters
  }
  # Read in features matrix and calculate differential expression by cluster.
  if (is(features, "character")) {
    features <- read.csv(features, row.names = 1, check.names = FALSE)
  }
  features <- features[, colnames(dom@z_scores)]
  dom@features <- as.matrix(features)
  if (tf_selection_method == "clusters") {
    p_vals <- matrix(1, nrow = nrow(features), ncol = length(levels(dom@clusters)))
    rownames(p_vals) <- rownames(features)
    colnames(p_vals) <- levels(dom@clusters)
    if (verbose) {
      message("Calculating feature enrichment by cluster")
      clust_n <- length(levels(dom@clusters))
    }
    for (clust in levels(dom@clusters)) {
      if (verbose) {
        cur <- which(levels(dom@clusters) == clust)
        message(paste0(cur, " of ", clust_n))
      }
      cells <- which(dom@clusters == clust)
      for (feat in rownames(dom@features)) {
        p_vals[feat, clust] <- stats::wilcox.test(
          dom@features[feat, cells], dom@features[feat, -cells],
          alternative = "g"
        )$p.value
      }
    }
    dom@clust_de <- p_vals
  }
  if (tf_selection_method == "all") {
    dom@clusters <- factor()
  }
  if (tf_selection_method == "variable") {
    dom@clusters <- factor()
    variances <- apply(dom@features, 1, function(x) {
      sd(x) / mean(x)
    })
    keep_n <- length(variances) * tf_variance_quantile
    keep_id <- which(rank(variances) > keep_n)
    dom@features <- dom@features[names(keep_id), ]
  }
  # store tf_targets in linkages if they are provided as a list
  if (!is(tf_targets, "list")) {
    dom@linkages[["tf_targets"]] <- NULL
    warning("tf_targets is not a list. No regulons stored")
  } else {
    dom@linkages[["tf_targets"]] <- tf_targets
  }
  # Calculate correlation matrix between features and receptors.
  dom@counts <- counts
  zero_sum <- rowSums(counts == 0)
  keeps <- which(zero_sum < (1 - rec_min_thresh) * ncol(counts))
  ser_receptors <- intersect(names(keeps), rec_genes)
  rho <- matrix(0, nrow = length(ser_receptors), ncol = nrow(dom@features))
  rownames(rho) <- ser_receptors
  colnames(rho) <- rownames(dom@features)
  if (verbose) {
    message("Calculating correlations")
    n_tf <- nrow(dom@features)
  }
  for (module in rownames(dom@features)) {
    # If df is provided then check if receptors are targets of TF. If they are then set
    # correlation equal to 0.
    if (verbose) {
      cur <- which(rownames(dom@features) == module)
      message(paste0(cur, " of ", n_tf))
    }
    if (!is.null(dom@linkages$tf_targets)) {
      tf <- gsub(pattern = "\\.\\.\\.", replacement = "", module) # correction for AUC values from pySCENIC that append an elipses to TF names due to (+) characters in the orignial python output
      module_targets <- tf_targets[[tf]]
      module_rec_targets <- intersect(module_targets, ser_receptors)
    } else {
      module_rec_targets <- NULL
    }
    scores <- dom@features[module, ]
    rhorow <- rep(0, length(ser_receptors))
    names(rhorow) <- ser_receptors
    for (rec in ser_receptors) {
      if (remove_rec_dropout) {
        keep_id <- which(dom@counts[rec, ] > 0)
        rec_z_scores <- dom@z_scores[rec, keep_id]
        tar_tf_scores <- scores[keep_id]
      } else {
        rec_z_scores <- dom@z_scores[rec, ]
        tar_tf_scores <- scores
      }
      # There are some cases where all the tfs are zero for the cells left after trimming
      # dropout for receptors. Skip those and set cor to zero manually.
      if (sum(tar_tf_scores) == 0) {
        rhorow[rec] <- 0
        next
      }
      cor <- stats::cor.test(rec_z_scores, tar_tf_scores, method = "spearman", alternative = "greater")
      rhorow[rec] <- cor$estimate
    }
    if (length(module_rec_targets > 0)) {
      rhorow[module_rec_targets] <- 0
    }
    rho[, module] <- rhorow
  }
  colnames(rho) <- rownames(dom@features)
  dom@misc$rec_cor <- rho
  # assess correlation among genes in the same receptor complex
  cor_list <- list()
  for (i in 1:length(names(dom@linkages$rec_lig))) {
    r <- names(dom@linkages$rec_lig)[i]
    if (r %in% names(dom@linkages$complexes)) {
      r_genes <- dom@linkages$complexes[[r]]
    } else {
      r_genes <- r
    }
    if (sum(rownames(rho) %in% r_genes) != length(r_genes)) {
      message(paste0(r, " has component genes that did not pass testing parameters"))
      cor_list[[r]] <- rep(0, ncol(rho))
      next
    }
    if (length(r_genes) > 1) {
      gene_cor <- rho[rownames(rho) %in% r_genes, ]
      cor_med <- apply(gene_cor, 2, function(x) {
        median(x)
      })
      cor_list[[r]] <- cor_med
    } else {
      cor_list[[r]] <- rho[rownames(rho) == r_genes, ]
    }
  }
  c_cor <- t(as.data.frame(cor_list))
  dom@cor <- c_cor
  # If cluster methods are used, calculate percentage of non-zero expression of receptor genes
  # in clusters
  if (tf_selection_method == "clusters") {
    cl_rec_percent <- NULL
    for (rec in ser_receptors) {
      rec_percent <- sapply(X = levels(dom@clusters), FUN = function(x) {
        # percentage of cells in cluster with non-zero expression of receptor gene
        sum(dom@counts[rec, dom@clusters == x] > 0) / length(dom@counts[rec, dom@clusters ==
          x])
      })
      cl_rec_percent <- rbind(cl_rec_percent, rec_percent)
    }
    rownames(cl_rec_percent) <- ser_receptors
    dom@misc$cl_rec_percent <- cl_rec_percent
  }
  return(dom)
}

#' Adds a column to the RL signaling data frame.
#'
#' This function adds a column to the internal rl 'map' used to map all
#' receptor and receptor complexes to all ligand and ligand complexes.
#'
#' @param map RL signaling data frame.
#' @param map_ref Name of column to match new data to
#' @param conv Data frame matching current data in map to new data.
#' @param new_name Name of new column to be created in RL map
#' @return An updated RL signaling data frame
#' @export
#' @examples 
#' lr_name <- data.frame("abbrev" = c("L", "R"), "full" = c("Ligand", "Receptor"))
#' rl_map_expanded <- add_rl_column(map = domino2:::rl_map_tiny, map_ref = "type_A",
#'  conv = lr_name, new_name = "type_A_full")
#' 
add_rl_column <- function(map, map_ref, conv, new_name) {
  map_in_ref <- match(map[[map_ref]], conv[, 1])
  not_in_ref <- which(is.na(map_in_ref))
  if (length(not_in_ref > 0)) {
    not_in_ref_map <- cbind.data.frame(map[not_in_ref, ], as.character(NA), stringsAsFactors = FALSE)
    colnames(not_in_ref_map)[ncol(not_in_ref_map)] <- new_name
    rownames(not_in_ref_map) <- c()
  } else {
    not_in_ref_map <- c()
  }
  new_map <- c()
  for (r_id in 1:nrow(map)) {
    row <- map[r_id, ]
    conv_ids <- which(conv[, 1] == row[[map_ref]])
    for (id in conv_ids) {
      new_row <- c(as.matrix(row), conv[id, 2])
      new_map <- rbind(new_map, new_row)
    }
  }
  rownames(new_map) <- c()
  colnames(new_map) <- c(colnames(map), new_name)
  new_map <- rbind.data.frame(new_map, not_in_ref_map, stringsAsFactors = FALSE)
  new_map <- data.frame(new_map, stringsAsFactors = FALSE)
}

#' Calculate mean ligand expression as a data.frame for plotting in circos plot
#'
#' Creates a data frame of mean ligand expression for use in plotting a circos
#' plot of ligand expression and saving tables of mean expression.
#'
#' @param x Gene by cell expression matrix
#' @param ligands Character vector of ligand genes to be quantified
#' @param cell_ident Vector of cell type (identity) names for which to calculate mean ligand gene expression
#' @param cell_barcodes Vector of cell barcodes (colnames of x) belonging to cell_ident to calculate mean expression across
#' @param destination Name of the receptor with which each ligand interacts
#' @return A data frame of ligand expression targeting the specified receptor
#' @export
#' @examples
#' counts <- dom_counts(domino2:::pbmc_dom_built_tiny)
#' mean_exp <- mean_ligand_expression(counts,
#'  ligands = c("PTPRC", "FASLG"), cell_ident = "CD14_monocyte",
#'  cell_barcodes = colnames(counts), destination = "FAS")
#' 
mean_ligand_expression <- function(x, ligands, cell_ident, cell_barcodes, destination){
  # initiate data frame to store results
  df <- NULL
  for(feat in ligands){
    # index of ligand row
    lig_index <- grep(paste0("^", feat, "$"), rownames(x))
    # column indecies of cells belonging to cell_ident
    cell_index <- colnames(x) %in% cell_barcodes
    cell_df <- data.frame(
      origin = paste0(cell_ident, "_", feat),
      destination = destination,
      mean.expression = mean(x[lig_index, cell_index])
    )
    df <- rbind(df, cell_df)
  }
  return(df)
}
