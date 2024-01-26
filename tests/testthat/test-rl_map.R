# read_cellphonedb_file
test_that("read_cellphonedb_file: replace empty cells with FALSE", {
  tdir <- tempdir()
  df <- data.frame(
    "uniprot" = c("A", "B", "C"), 
    "transmembrane" = c("True", "False", ""), 
    "peripheral" = c("True", "False", ""), 
    "secreted" = c("True", "False", ""), 
    "secreted_highlight" = c("True", "False", ""), 
    "receptor" = c("True", "False", ""), 
    "integrin" = c("True", "False", ""), 
    "other" = c("True", "False", "")
  )
  df_check <- data.frame(
    "uniprot" = c("A", "B", "C"), 
    "transmembrane" = c("TRUE", "FALSE", "FALSE"), 
    "peripheral" = c("TRUE", "FALSE", "FALSE"), 
    "secreted" = c("TRUE", "FALSE", "FALSE"), 
    "secreted_highlight" = c("TRUE", "FALSE", "FALSE"), 
    "receptor" = c("TRUE", "FALSE", "FALSE"), 
    "integrin" = c("TRUE", "FALSE", "FALSE"), 
    "other" = c("TRUE", "FALSE", "FALSE")
  )
  test_file <- paste0(tdir, "/test_prot.csv")
  write.csv(df, file = test_file, row.names = FALSE)
  df2 <- read_cellphonedb_file(filename = test_file, set_unannotated = TRUE)
  expect_equal(df2, df_check)
})

# TF_syntax_to_R
test_that("TF_syntax_to_R: Convert python True/False syntax to TRUE/FALSE", {
  x <- c("True", "True", "False", "True")
  y <- TF_syntax_to_R(x = x)
  expect_equal(
    y, 
    c(TRUE, TRUE, FALSE, TRUE)
  )
})

# check_filepath
test_that("check_filepath: read file and ignore non character vectors", {
  tdir <- tempdir()
  df <- data.frame(
    "uniprot" = c("A", "B", "C"), 
    "transmembrane" = c("True", "False", "")
  )
  df_check <- data.frame(
    "uniprot" = c("A", "B", "C"), 
    "transmembrane" = c("TRUE", "FALSE", "FALSE")
  )
  test_file <- paste0(tdir, "/test_prot.csv")
  write.csv(df, file = test_file, row.names = FALSE)
  df2 <- check_filepath(test_file)
  expect_equal(df2, df_check)
  non_df <- check_filepath(1)
  expect_equal(non_df, 1)
})

# parse_interaction
test_that("parse_interaction: format simple or complpex partners as list", {
  test_int_AB <- c(
    "partner_a" = "simpleA",
    "partner_b" = "complexB",
    "protein_name_a" = "PROT_A",
    "protein_name_b" = "",
    "annotation_strategy" = "test",
    "source" = "test"
  )
  test_int_AC <- c(
    "partner_a" = "simpleA",
    "partner_b" = "simpleC",
    "protein_name_a" = "PROT_A",
    "protein_name_b" = "PROT_C",
    "annotation_strategy" = "test",
    "source" = "test"
  )
  test_genes <- data.frame(
    "gene_name" = c("GENEA", "GENEB1", "GENEB2", "GENEC"),
    "uniprot" = c("simpleA", "complexB1", "complexB2", "simpleC"),
    "hgnc_symbol" = c("GA", "GB1", "GB2", "GC"),
    "ensembl" = c("ENSGA", "ENSGB1", "ENSGB2", "ENSGC")
  )
  test_complexes <- data.frame(
    "complex_name" = "complexB",
    "uniprot_1" = "complexB1",
    "uniprot_2" = "complexB2",
    "uniprot_3" = "",
    "uniprot_4" = "",
    "receptor" = TRUE
  )
  test_proteins <- data.frame(
    "uniprot" = c("simpleA", "complexB1", "complexB2", "simpleC"),
    "protein_name" = c("PROT_A", "PROT_B1", "PROT_B2", "PROT_C"),
    "receptor" = c(FALSE, TRUE, TRUE, TRUE)
  )
  lig_int <- parse_interaction(
    int = test_intAB, partner = "A",
    genes = test_genes, proteins = test_proteins,
    complexes = test_complexes
  )
  rec_int <- parse_interaction(
    int = test_intAB, partner = "B",
    genes = test_genes, proteins = test_proteins,
    complexes = test_complexes
  )
  rec_int_noComplex <- parse_interaction(
    int = test_intAB, partner = "B",
    genes = test_genes, proteins = test_proteins,
    complexes = NULL
  )
  expect_equal(
    lig_int,
    list("gene_A" = "GENEA", "type_A" = "L", "name_A" = "GENEA")
  )
  expect_equal(
    rec_int,
    list("gene_B" = "GENEB1,GENEB2", "type_B" = "R", "name_B" = "complexB")
  )
  expect_equal(
    rec_int_noComplex,
    list("gene_B" = "", "type_B" = "", "name_B" = "complexB")
  )
})

# create_rl_map_cellphonedb
test_that("create_rl_map_cellphonedb: using interactions with and without complexes", {
  test_interactions <- data.frame(
    "partner_a" = c("simpleA", "simpleA"),
    "partner_b" = c("complexB", "simpleC"),
    "protein_name_a" = c("PROT_A", "PROT_A"),
    "protein_name_b" = c("", "PROT_C"),
    "annotation_strategy" = c("test", "test"),
    "source" = c("test", "test")
  )
  test_genes <- data.frame(
    "gene_name" = c("GENEA", "GENEB1", "GENEB2", "GENEC"),
    "uniprot" = c("simpleA", "complexB1", "complexB2", "simpleC"),
    "hgnc_symbol" = c("GA", "GB1", "GB2", "GC"),
    "ensembl" = c("ENSGA", "ENSGB1", "ENSGB2", "ENSGC")
  )
  test_complexes <- data.frame(
    "complex_name" = "complexB",
    "uniprot_1" = "complexB1",
    "uniprot_2" = "complexB2",
    "uniprot_3" = "",
    "uniprot_4" = "",
    "receptor" = TRUE
  )
  test_proteins <- data.frame(
    "uniprot" = c("simpleA", "complexB1", "complexB2", "simpleC"),
    "protein_name" = c("PROT_A", "PROT_B1", "PROT_B2", "PROT_C"),
    "receptor" = c(FALSE, TRUE, TRUE, TRUE)
  )
  
  rl_map_complexes <- create_rl_map_cellphonedb(
    genes = test_genes, proteins = test_proteins, interactions = test_interactions,
    complexes = test_complexes
  )
  rl_map_NOcomplexes <- create_rl_map_cellphonedb(
    genes = test_genes, proteins = test_proteins, interactions = test_interactions,
    complexes = NULL
  )
  
  rl_map_check <- data.frame(
    gene_A = c("GENEA", "GENEA"),
    type_A = c("L", "L"),
    name_A = c("GENA", "GENEA"),
    gene_B = c("GENEB1,GENEB2", "GENEC"),
    type_B = c("R", "R"),
    name_B = c("complexB", "GENEC"),
    int_pair = c("GENEA & complexB", "GENEA & GENEC"),
    annotation_strategy = c("test", "test"),
    source = c("test", "test"),
    database_name = c("CellPhoneDB", "CellPhoneDB")
  )
  expect_equal(rl_map_complexes, rl_map_check)
  expect_equal(rl_map_NOcomplexes, rl_map_check[2,])
})

# ortholog_mapping 
test_that("ortholog_mapping", {
  test_from <- c("HUMA1", "HUMB1", "HUMB1", "HUMC1", "HUMD1", "HUMD2")
  test_to <- c("Musa1", "Musb1", "Musb2", "", "Musd1", "Musd1")
  
  map_OneOne <- ortholog_mapping(gene = "HUMA1", from = test_from, to = test_to)
  map_OneMany <- ortholog_mapping(gene = "HUMB1", from = test_from, to = test_to)
  map_ManyOne <- ortholog_mapping(gene = c("HUMD1", "HUMD2"), from = test_from, to = test_to)
  map_OneNone <- ortholog_mapping(gene = "HUMC1", from = test_from, to = test_to)
  
  expect_equal(map_OneOne, list(c("Musa1")))
  expect_equal(map_OneMany, list(c("Musb1", "Musb2")))
  expect_equal(map_ManyOne, list(c("Musd1"), c("Musd1")))
  expect_equal(map_OneNone, list(c("")))
  
  # trying to find an ortholog for a gene not included in 'from'
  map_missing <- ortholog_mapping(gene = "HUME", from = test_from, to = test_to)
  expect_equal(map_missing, list(c(NA)))
  
  # return error when 'from' and 'to' are unequal lengths
  expect_error(ortholog_mapping(gene = "HUMA1", from = test_from[1:2], to = test_to[1:3]))
})

# rl_map_ortholog_conversion
test_that("rl_map_ortholog_conversion: cases for simple interactions, complex interactions, and running without complexes", {
  test_rl_map <- data.frame(
    gene_A = c("HUMA1", "HUMB1", "HUMC1", "HUMD1", "HUMD2", "HUMA1,HUMA2", "HUMB1,HUMA2", "HUMC1,HUMA2", "HUMD1,HUMA2", "HUMD2,HUMA2"),
    type_A = rep("R", 10),
    name_A = c("HUMA1", "HUMB1", "HUMC1", "HUMD1", "HUMD2", "complexAA", "complexBA", "complexCA", "complexD1A", "complexD2A"),
    gene_B = rep("HUME", 10),
    type_B = rep("L", 10),
    name_B = rep("HUME", 10),
    int_pair = c("HUMA1 & HUME", "HUMB1 & HUME", "HUMC1 & HUME", "HUMD1 & HUME", "HUMD2 & HUME", "complexAA & HUME", "complexBA & HUME", "complexCA & HUME", "complexD1A & HUME", "complexD2A & HUME"),
    annotation = rep("test", 10),
    source = c("simple int, one-one ortho", "simple int, one-many ortho", "simple int, one-none ortho", "simple int, many-one ortho", "simple int, many-one ortho", "complex int, one-one ortho", "complex int, one-many ortho", "complex int, one-none ortho", "complex int, many-one ortho", "complex int, many-one ortho"),
    database_name = rep("test", 10)
  )
  test_from <- c("HUMA1", "HUMA2", "HUMB1", "HUMB1", "HUMC1", "HUMD1", "HUMD2", "HUME")
  test_to <- c("Musa1", "Musa2", "Musb1", "Musb2", "", "Musd1", "Musd1", "Muse")
  
  ortho_rl_map <- rl_map_ortholog_conversion(rl_map = test_rl_map, to = test_to, from = test_from, use_complexes = TRUE)
  ortho_rl_map_noComplex <- rl_map_ortholog_conversion(rl_map = test_rl_map, to = test_to, from = test_from, use_complexes = FALSE)
  
  # simple int, one-one ortho
  expect_equal(
    ortho_rl_map[ortho_rl_map$name_A == "HUMA1",],
    c(
      gene_A = "Musa1",
      type_A = "R",
      name_A = "HUMA1",
      gene_B = "Muse",
      type_B = "L",
      name_B = "HUME",
      int_pair = "HUMA1 & HUME",
      annotation = "test",
      source = "simple int, one-one ortho",
      database_name = "test_conversion_"
    )
  )
  # simple int, one-many ortho
  # expectation is to include an interaction for each ortholog
  expect_equal(
    ortho_rl_map[ortho_rl_map$name_A == "HUMB1",],
    c(
      gene_A = c("Musb1", "Musb2"),
      type_A = c("R", "R"),
      name_A = c("HUMB1", "HUMB1"),
      gene_B = c("Muse", "Muse"),
      type_B = c("L", "L"),
      name_B = c("HUME", "HUME"),
      int_pair = c("HUMB1 & HUME", "HUMB1 & HUME"),
      annotation = c("test", "test"),
      source = c("simple int, one-many ortho", "simple int, one-many ortho"),
      database_name = c("test_conversion_", "test_conversion_")
    )
  )
  # simple int, one-none ortho
  expect_equal(
    dim(ortho_rl_map[ortho_rl_map$name_A == "HUMC1",]),
    c(0, 10)
  )
  # simple int, many-one ortho
  expect_equal(
    ortho_rl_map[ortho_rl_map$name_A %in% c("HUMD1", "HUMD2"),],
    c(
      gene_A = c("Musd1", "Musd1"),
      type_A = c("R", "R"),
      name_A = c("HUMD1", "HUMD2"),
      gene_B = c("Muse", "Muse"),
      type_B = c("L", "L"),
      name_B = c("HUME", "HUME"),
      int_pair = c("HUMD1 & HUME", "HUMD2 & HUME"),
      annotation = c("test", "test"),
      source = c("simple int, many-one ortho", "simple int, many-one ortho"),
      database_name = c("test_conversion_", "test_conversion_")
    )
  )
  # complex int, one-one ortho
  expect_equal(
    ortho_rl_map[ortho_rl_map$name_A == "complexAA",],
    c(
      gene_A = "Musa1,Musa2",
      type_A = "R",
      name_A = "complexAA",
      gene_B = "Muse",
      type_B = "L",
      name_B = "HUME",
      int_pair = "complexAA & HUME",
      annotation = "test",
      source = "complex int, one-one ortho",
      database_name = "test_conversion_"
    )
  )
  # complex int, one-many ortho
  expect_equal(
    ortho_rl_map[ortho_rl_map$name_A == "complexBA",],
    c(
      gene_A = c("Musb1,Musa2", "Musb2,Musa2"),
      type_A = c("R", "R"),
      name_A = c("complexBA", "complexBA"),
      gene_B = c("Muse", "Muse"),
      type_B = c("L", "L"),
      name_B = c("HUME", "HUME"),
      int_pair = c("complexBA & HUME", "complexBA & HUME"),
      annotation = c("test", "test"),
      source = c("complex int, one-many ortho", "complex int, one-many ortho"),
      database_name = c("test_conversion_", "test_conversion_")
    )
  )
  # complex int, one-none ortho
  expect_equal(
    dim(ortho_rl_map[ortho_rl_map$name_A == "complexCA",]),
    c(0, 10)
  )
  # complex int, many-one ortho
  expect_equal(
    ortho_rl_map[ortho_rl_map$name_A %in% c("complexD1A", "complexD2A"),],
    c(
      gene_A = c("Musd1,Musa2", "Musd1,Musa2"),
      type_A = c("R", "R"),
      name_A = c("complexD1A", "complexD2A"),
      gene_B = c("Muse", "Muse"),
      type_B = c("L", "L"),
      name_B = c("HUME", "HUME"),
      int_pair = c("complexD1A & HUME", "complexD2A & HUME"),
      annotation = c("test", "test"),
      source = c("complex int, many-one ortho", "complex int, many-one ortho"),
      database_name = c("test_conversion_", "test_conversion_")
    )
  )
  
  # running without complexes returns only simple interactions
  expect_equal(
    dim(ortho_rl_map_noComplex),
    c(5, 10)
  )
})

