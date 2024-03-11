# internal scripts for the create_domino() function

neo_create_domino <- function(
    rl_map, features, ser = NULL, counts = NULL, z_scores = NULL,
    clusters = NULL, use_clusters = TRUE, tf_targets = NULL, verbose = TRUE,
    use_complexes = TRUE, rec_min_thresh = 0.025, remove_rec_dropout = TRUE,
    tf_selection_method = "clusters", tf_variance_quantile = 0.5){
  # check function inputs
  check_arg(rl_map, allow_class = "data.frame",
            need_vars = c("gene_A", "gene_B", "type_A", "type_B"))
  check_arg(features, allow_class = c("data.frame", "character", "matrix"))
  if (any(class(features) %in% c("data.frame", "matrix"))) {
    check_arg(features, need_rownames = TRUE, need_colnames = TRUE)
  }
  if (!is.null(ser)) {
    check_arg(ser, allow_class = "Seurat")
  } else {
    check_arg(counts, allow_class = c("matrix", "data.frame", "Matrix", "dgCMatrix"),
              need_rownames = TRUE, need_colnames = TRUE)
    check_arg(z_scores, allow_class = "matrix", need_rownames = TRUE,
              need_colnames = TRUE)
    check_arg(clusters, allow_class = "factor", need_names = TRUE)
  }
  check_arg(rec_min_thresh, allow_class = c("numeric"), allow_range = c(0, 1))
  check_arg(tf_selection_method,
            allow_values = c("clusters", "variable", "all"))
  if (!is.null(ser) & (!is.null(clusters) | !is.null(z_scores) | !is.null(counts))) {
    warning("Ser and z_score, clusters, or counts provided. Defaulting to ser.")
  }
  # rl_map processing
  rl_parse <- read_rl_map_genes(rl_map)
  linkage_rec_lig <- write_rec_lig_linkages(rl_parse = rl_parse)
  linkage_complexes <- read_rl_map_complexes(rl_parse = rl_parse, use_complexes = use_complexes)
  receptor_names <- rl_parse[["R.name"]]
  receptor_genes <- unique(unlist(strsplit(rl_parse[["R.gene"]], split = "\\,")))
  
  # Establish TF-rec linkages
  cluster_de <- select_cluster_tf(
    features = features, clusters = clusters, 
    method = "one.sided.wilcox", verbose = verbose
  )
  tfs_rec_cor <- test_tfs_rec_linkage(
    features = features, z_scores = z_scores, counts = counts,
    feature_de = cluster_de, receptors = receptor_names, 
    method = "spearman.correlation", verbose = verbose
  )
  tfs_rec_cor_filter <- filter_tf_regulon_receptors(cor_mat = tfs_rec_cor, tf_targets = tf_targets)
  if(use_complexes) {
    tfs_rec_cor_filter <- assess_complex_receptor_cor(receptors = receptor_names, complexes_list = linkage_complexes, cor_mat = tfs_rec_cor)
  }
  cl_rec_precent <- calc_rec_percentage(counts = counts, clusters = clusters, receptor_genes = receptor_genes)
  
  # create domino result
  dom <- new(
    "domino",
    db_info = rl_map,
    z_scores = z_scores,
    counts = counts,
    clusters = clusters,
    features = features,
    cor = tfs_rec_cor_filter,
    linkages = list(
      "complexes" = linkage_complexes,
      "rec_lig" = linkage_rec_lig,
      "tf_targets" = tf_targets
    ),
    clust_de = cluster_de,
    misc = list(
      "create" = TRUE,
      "build" = FALSE,
      "build_vars" = NULL,
      "rec_cor" = tfs_rec_cor,
      "cl_rec_percent" = cl_rec_precent
    )
  )
  return(dom)
}

calc_rec_percentage <- function(counts, clusters, receptor_genes) {
  cluster_lvls <- levels(clusters)
  names(receptor_genes) <- receptor_genes
  cluster_percent_ls <- lapply(
    receptor_genes, FUN = function(r) {
      r_row_logic <- rownames(counts) == r
      rec_percent <- vapply(
        cluster_lvls, FUN = function(cl) {
          cell_in_cl <- clusters == cl
          expressing_cells <- sum(counts[r_row_logic, cell_in_cl] > 0)
          rec_prop <- expressing_cells / length(cell_in_cl) 
          return(rec_prop)
      }, 
      FUN.VALUE = numeric(length(1))
      )
    }
  )
  cluster_percent <- do.call(rbind, cluster_percent_ls)
  return(cluster_percent)
 }

assess_complex_receptor_cor <- function(receptors, complexes_list, cor_mat, method = "median") {
  valid_methods <- c("median")
  if(!(method %in% valid_methods)){
    stop("Invalid method supplied")
  }
  names(receptors) <- receptors
  complexes <- names(complexes_list)
  complex_cor_ls <- lapply(
    receptors, FUN = function(r) {
      if(r %in% complexes) {
        r_comp <- complexes_list[[r]]
        # determine if all components are present
        r_comp_logic <- r_comp %in% rownames(cor_mat)
        if(sum(r_comp_logic) == length(r_comp_logic)) {
          gene_cor <- cor_mat[rownames(cor_mat) %in% r_comp, ]
          if(method == "median") {
            cor_aggr <- apply(gene_cor, 2, function(x) {
              median(x)
            })
          }
          return(cor_aggr)
        } else {
          return()
        }
      } else {
        # obtain rows of correlation matrix for non-complex receptors
        if(r %in% rownames(cor_mat)) {
          r_cor <- cor_mat[r,]
          return(r_cor)
        } else {
          return()
        }
      }
    }
  )
  complex_cor <- do.call(rbind, complex_cor_ls)
  colnames(complex_cor) <- colnames(cor_mat)
  return(complex_cor)
}

test_tfs_rec_linkage <- function(
    features, z_scores, counts, 
    feature_de,
    receptors, 
    method = "spearman.correlation", verbose = TRUE) {
  tfs <- rownames(feature_de)
  names(tfs) <- tfs
  n_tfs <- length(tfs)
  rec_z_scores <- z_scores[rownames(z_scores) %in% receptors,]
  rec_counts <- counts[rownames(counts) %in% receptors,]
  confirmed_recs <- rownames(rec_z_scores)
  if(verbose) {message("Calculating feature-receptor linkages")}
  cl_link_list <- lapply(
    tfs,
    function(tf) {
      a <- which(tf == tfs)
      if(verbose) {message(paste0(a, " of ", n_tfs))}
      tf_scores <- features[tf,]
      linkage_score_ls <- lapply(
        confirmed_recs,
        function(r) {
          r_exp <- rec_z_scores[r,]
          if(method == "spearman.correlation") {
            # warnings supressed for cases where exact p-values cannot be calculated with ties
            test_res <- suppressWarnings(stats::cor.test(r_exp, tf_scores, method = "spearman", alternative = "greater"))
            est <- test_res[["estimate"]]
            names(est) <- r
            return(est)
          }
        }
      )
      linkage_score <- unlist(linkage_score_ls)
      return(linkage_score)
    }
  )
  linkage_score_mat <- do.call(cbind, cl_link_list)
  colnames(linkage_score_mat) <- tfs
  return(linkage_score_mat)
}

filter_tf_regulon_receptors <- function(cor_mat, tf_targets) {
  tfs <- names(tf_targets)
  recs <- rownames(cor_mat)
  cor_prune_ls <- lapply(
    tfs, FUN = function(tf) {
      regulon <- tf_targets[[tf]]
      tf_col <- cor_mat[,tf]
      recs <- names(tf_col)
      regulon_logic <- recs %in% regulon
      tf_col[regulon_logic] <- 0
      return(tf_col)
    }
  )
  cor_prune <- do.call(cbind, cor_prune_ls)
  colnames(cor_prune) <- tfs
  return(cor_prune)
}

select_cluster_tf <- function(features, clusters, verbose = TRUE, method = "one.sided.wilcox") {
  cluster_lvls <- levels(clusters)
  n_lvls <- length(cluster_lvls)
  feat_names <- rownames(features)
  names(feat_names) <- feat_names
  if(verbose) {message("Calculating feature enrichment by cluster")}
  cl_pval_list <- lapply(
    cluster_lvls,
    function(cl) {
      a <- which(cluster_lvls == cl)
      if(verbose) {message(paste0(a, " of ", n_lvls))}
      cl_cells <- which(clusters == cl)
      out_cells <- which(clusters != cl)
      p_val_ls <- lapply(
        feat_names,
        function(f) {
          cl_value <- features[f, cl_cells]
          out_value <- features[f, out_cells]
          if(method == "one.sided.wilcox") {
            test_res <- stats::wilcox.test(
              x = cl_value, y = out_value, alternative = "greater"
            )
            return(test_res[["p.value"]])
          }
        }
      )
      p_vals <- unlist(p_val_ls)
      return(p_vals)
    }
  )
  pval_mat <- do.call(cbind, cl_pval_list)
  colnames(pval_mat) <- cluster_lvls
  return(pval_mat)
}

write_rec_lig_linkages <- function(rl_parse) {
  rec_names <- unique(rl_parse[["R.name"]])
  names(rec_names) <- rec_names
  rec_lig_ls <- lapply(
    rec_names,
    function(r) {
      rl_int_row <- rl_parse[rl_parse[["R.name"]] == r,]
      ligs <- rl_int_row[["L.name"]]
      return(ligs)
    }
  )
  return(rec_lig_ls)
}

read_rl_map_complexes <- function(rl_parse, use_complexes) {
  if(use_complexes) {
    complex_lig <- rl_parse[grepl("\\,", rl_parse$L.gene), "L.name"]
    complex_rec <- rl_parse[grepl("\\,", rl_parse$R.gene), "R.name"]
    complex_names <- unique(c(complex_lig, complex_rec))
    names(complex_names) <- complex_names
    gene_complexes <- lapply(
      complex_names, FUN = function(x) {
        comp_cs <- c(
          rl_parse[rl_parse$L.name == x, "L.gene"],
          rl_parse[rl_parse$R.name == x, "R.gene"]
        )
        if(length(unique(comp_cs)) > 1) {
          comp_first <- comp_cs[1]
          warning(
            paste0(
              "protein complex ", x, 
              " has multiple, differing descriptions of component genes. \n",
              "Defaulting to first encountered component description: \n",
              comp_first
            )
          )
          comp <- unique(unlist(strsplit(comp_first, split = "\\,")))
          return(comp)
        } else {
          comp <- unique(unlist(strsplit(comp_cs, split = "\\,")))
          return(comp)
        }
      }
    )
    return(gene_complexes)
  } else {
    return(NULL)
  }
}

read_rl_map_genes <- function(rl_map) {
  rl_list <- apply(
    rl_map, MARGIN = 1,
    FUN = function(x) {
      rl <- c()
      p <- ifelse(x[["type_A"]] == "R", "A", "B")
      q <- ifelse(p == "A", "B", "A")
      R.gene <- x[[paste0("gene_", p)]]
      L.gene <- x[[paste0("gene_", q)]]
      rl[["R.gene"]] <- R.gene
      rl[["L.gene"]] <- L.gene
      if (paste0("uniprot_", p) %in% names(x)) {
        rl[["R.uniprot"]] <- x[[paste0("uniprot_", p)]]
      }
      if (paste0("uniprot_", q) %in% names(x)) {
        rl[["L.uniprot"]] <- x[[paste0("uniprot_", q)]]
      }
      if (paste0("name_", p) %in% names(x)) {
        rl[["R.name"]] <- x[[paste0("name_", p)]]
      }
      if (paste0("name_", q) %in% names(x)) {
        rl[["L.name"]] <- x[[paste0("name_", q)]]
      }
      return(as.data.frame(rl))
    }
  )
  rl_df <- as.data.frame(do.call(rbind, rl_list))
  return(rl_df)
}
