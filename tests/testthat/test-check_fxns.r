# unit tests for check_warning_fxn.R

test_that("incorrect inputs to create_domino are caught", {
  # expected formats
  good_barcodes <- c("AAA-1", "AAT-1", "AAC-1", "AAG-1", "ATA-1")
  good_rl_map <- data.frame(
    "gene_A" = c("HUMA", "HUMB"),
    "gene_B" = c("HUMZ", "HUMY"),
    "type_A" = c("L", "L"),
    "type_B" = c("R", "R")
  )
  good_auc <- matrix(
    rbinom(n = 25, size = 5, prob = 0.5),
    ncol = 5, nrow = 5,
    dimnames = list(
      LETTERS[1:5],
      good_barcodes
    )
  )
  good_counts <- Matrix::Matrix(
    data = rbinom(n = 25, size = 5, prob = 0.5),
    ncol = 5, nrow = 5, sparse = TRUE,
    dimnames = list(
      LETTERS[6:10],
      good_barcodes
    )
  )
  good_clusters <- factor(
    c(rep("X", 2), rep("Y", 3)), levels = c("X", "Y")
  )
  names(good_clusters) <- good_barcodes
  good_regulon_ls <- list(
    "A" = c("F"), "B" = c("G", "H"), "C" = c(""), "D" = c("I"), "E" = c("J")
  )
  good_ser <- matrix(
    1, nrow = 5, ncol = 5,
    dimnames = list(
      LETTERS[6:10],
      good_barcodes
    )
  )
  # missing data in rl_map
  rl_map_missing_col <- good_rl_map[,-1]
  rl_map_mat <- matrix(0, ncol = 4, nrow = 2)
  
  expect_error(
    check_create_domino(
      rl_map = rl_map_missing_col,
      ser = NULL, features = good_auc, counts = good_counts, clusters = good_clusters, tf_targets = good_regulon_ls
    )
  )
  expect_error(
    check_create_domino(
      rl_map = rl_map_mat,
      ser = NULL, features = good_auc, counts = good_counts, clusters = good_clusters, tf_targets = good_regulon_ls
    )
  )
  
  test_rl_map_df <- data.frame()
  test_rl_map_mat <- matrix(0, ncol = 4, nrow = 4)
})