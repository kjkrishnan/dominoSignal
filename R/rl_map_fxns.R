
#' Read CellPhoneDB file
#' 
#' Read CellPhoneDB file into R session as a data frame.
#' 
#' @param filename character. file path to CellPhoneDB data base file
#' @param set_unnanotated logical. Replace empty values in logical vector columns with FALSE
#' @export
#' @examples
#' tdir <- tempdir()
#' df <- data.frame(
#'   "uniprot" = c("A", "B", "C"), 
#'   "transmembrane" = c("True", "False", ""), 
#'   "peripheral" = c("True", "False", ""), 
#'   "secreted" = c("True", "False", ""), 
#'   "secreted_highlight" = c("True", "False", ""), 
#'   "receptor" = c("True", "False", ""), 
#'   "integrin" = c("True", "False", ""), 
#'   "other" = c("True", "False", "")
#' )
#' test_file <- paste0(tdir, "/test_prot.csv")
#' write.csv(df, file = test_file, row.names = FALSE)
#' read_cellphonedb_file(filename = test_file, set_unannotated = TRUE)
#' 

read_cellphonedb_file <- function(filename, set_unannotated = FALSE) {
  my_file <- read.csv(filename, stringsAsFactors = FALSE)
  if (set_unannotated) {
    # replace empty data frame cells with "False"
    gene_features <- c("transmembrane", "peripheral", "secreted", "secreted_highlight", "receptor", "integrin", "other")
    my_file[my_file$receptor == "", colnames(my_file) %in% gene_features] <- "False"
  }
  # change cases of True/False syntax from Python to TRUE/FALSE R syntax
  for (x in colnames(my_file)) {
    my_file[[x]] <- TF_syntax_to_R(my_file[[x]])
  }
  return(my_file)
}

#' Parse interactions
#' 
#' Use CellPhoneDB database files to parse each interactions and the genes comprising interacting ligand-receptor partners into lists genes comprising the partner, whether the partner is a ligand or recptor, and the name of the partner
#' 
#' @param int character vector describing the interaction. Values have names "partner_a", "partner_b", "protein_name_a", "protein_name_b", "annotation_strategy", "source".
#' @param genes data.frame. Genes table from CellPhoneDB.
#' @param proteins data.frame. Proteins table from CellPhoneDB.
#' @param partner character. Whether the function should return a list describing partner A or partner B in the interaction described in "int".
#' @param complexes data.frame. OPTIONAL. Complexes table from CellPhoneDB.
#' @return list
#' @export
#' @examples 
#' test_intAB <- c(
#'   "partner_a" = "simpleA", 
#'   "partner_b" = "complexB", 
#'   "protein_name_a" = "PROT_A", 
#'   "protein_name_b" = "", 
#'   "annotation_strategy" = "test",
#'   "source" = "test"
#'  )
#' test_genes <- data.frame(
#'   "gene_name" = c("GENEA", "GENEB1", "GENEB2", "GENEC"),
#'   "uniprot" = c("simpleA", "complexB1", "complexB2", "simpleC"),
#'   "hgnc_symbol" = c("GA", "GB1", "GB2", "GC"),
#'   "ensembl" = c("ENSGA", "ENSGB1", "ENSGB2", "ENSGC")
#'  )
#' test_complexes <- data.frame(
#'   "complex_name" = "complexB",
#'   "uniprot_1" = "complexB1",
#'   "uniprot_2" = "complexB2",
#'   "uniprot_3" = "",
#'   "uniprot_4" = "",
#'   "receptor" = TRUE
#'  )
#' test_proteins <- data.frame(
#'   "uniprot" = c("simpleA", "complexB1", "complexB2", "simpleC"),
#'   "protein_name" = c("PROT_A", "PROT_B1", "PROT_B2", "PROT_C"),
#'   "receptor" = c(FALSE, TRUE, TRUE, TRUE)
#'  )
#' parse_interaction(
#'   int = test_intAB, partner = "A",
#'   genes = test_genes, proteins = test_proteins,
#'   complexes = test_complexes
#' )
#' # $gene_A
#' # [1] "GENEA"
#' # 
#' # $type_A
#' # [1] "L"
#' # 
#' # $name_A
#' # [1] "GENEA"
#' 

parse_interaction <- function(int, genes, proteins, partner = c("A", "B"), complexes = NULL) {
  res <- list()
  p <- paste0("partner_", tolower(partner))
  pn <- paste0("protein_name_", tolower(partner))
  feat <- int[[p]]
  prot_name <- int[[pn]]
  if(feat %in% complexes[["complex_name"]]) {
    cpx <- complexes[complexes[["complex_name"]] == feat, ]
    comp_prot <- cpx[, c("uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4")]
    comp_prot <- comp_prot[comp_prot != ""]
    comp_gene <- vapply(
      comp_prot,
      FUN = function(x) {
        g <- unique(genes[genes[["uniprot"]] == x, c("gene_name")])
        if(length(g) > 1) {
          warning(paste0(
            "Multiple encoding genes for protein '", x, "' found in genes table.\n",
            "Defaulting to first encountered gene [", g[1], "]"
          ))
          g <- g[1]
        }
        return(g)
      },
      FUN.VALUE = character(1)
    )
    res[[paste0("gene_", partner)]] <- paste(comp_gene, collapse = ",")
    res[[paste0("protein_", partner)]] <- paste(comp_prot, collapse = ",")
    res[[paste0("type_", partner)]] <- ifelse(cpx[["receptor"]], "R", "L")
    res[[paste0("name_", partner)]] <- gsub(" ", "_", feat)
  } else if (feat %in% proteins[["uniprot"]]) {
    prot <- proteins[proteins[["uniprot"]] == feat, ]
    comp_prot <- prot[["uniprot"]]
    comp_gene <- unique(genes[genes[["uniprot"]] == comp_prot, c("gene_name")])
    if(length(comp_gene) > 1) {
      warning(paste0(
        "Multiple encoding genes for protein '", comp_prot, "' found in genes table.\n",
        "Defaulting to first encountered gene [", comp_gene[1], "]"
      ))
      comp_gene <- comp_gene[1]
    }
    res[[paste0("gene_", partner)]] <- comp_gene
    res[[paste0("protein_", partner)]] <- comp_prot
    res[[paste0("type_", partner)]] <- ifelse(prot[["receptor"]], "R", "L")
    res[[paste0("name_", partner)]] <- comp_gene
  } else {
    res[[paste0("gene_", partner)]] <- ""
    res[[paste0("protein_", partner)]] <- ""
    res[[paste0("type_", partner)]] <- ""
    res[[paste0("name_", partner)]] <- feat
  }
  return(res)
}

#' create receptor-ligand map (CellPhoneDB)
#' 
#' Re-formats the CellPhoneDB data base into a single data frame where each row 
#' describes a possible receptor-lignad interaction, the genes encoding the 
#' partners in the interaction, whether the partner is a liand or recptor, the 
#' name of the partner, how the interaction was annotated and the source for the
#' annotation, and the name of the database the rl_map was constructed from. 
#' Ligands and receptors that function as protein complexes may be included. 
#' Interactions involving these complexes are included as comma-seperated 
#' character strings of genes that comprise these complexes.
#' 
#' @param genes data.frame or file path to table of gene names in uniprot, hgnc_symbol, or ensembl format in cellphonedb database format
#' @param proteins data.frame or file path to table of protein features in cellphonedb format
#' @param interactions data.frame or file path to table of protein-protein interactions in cellphonedb format
#' @param complexes optional: data.frame or file path to table of protein complexes in cellphonedb format
#' @param database_name name of the database being used, stored in output. Default is "CellPhoneDB"
#' @return data.frame
#' @export
#' @examples
#' rl_map_tiny <- create_rl_map_cellphonedb(genes = domino2:::genes_tiny, 
#'  proteins = domino2:::proteins_tiny, interactions = domino2:::interactions_tiny, 
#'  complexes = domino2:::complexes_tiny)
#' 

create_rl_map_cellphonedb <- function(
    genes, proteins, interactions, 
    complexes = NULL, 
    database_name = "CellPhoneDB") {
  genes <- check_filepath(genes)
  proteins <- check_filepath(proteins, set_unannotated = TRUE)
  interactions <- check_filepath(interactions)
  if(!is.null(complexes)) {complexes <- check_filepath(complexes)}
  rl_map_ls <- list()
  for(i in seq(nrow(interactions))){
    int <- interactions[i,]
    partner_a <- parse_interaction(
      int=int, genes = genes, proteins = proteins, complexes = complexes,
      partner = "A"
    )
    partner_b <- parse_interaction(
      int=int, genes = genes, proteins = proteins, complexes = complexes,
      partner = "B"
    )
    if("" %in% partner_a | "" %in% partner_b) next
    if(
      (partner_a$type_A == "R" & partner_b$type_B == "R") |
      (partner_a$type_A == "L" & partner_b$type_B == "L")
    ) next
    i_features <- c(partner_a, partner_b)
    i_features[["int_pair"]] <- paste(i_features[["name_A"]], i_features[["name_B"]], sep = " & ")
    i_features[["annotation_strategy"]] <- int[["annotation_strategy"]]
    i_features[["source"]] <- int[["source"]]
    i_features[["database_name"]] <- database_name
    rl_map_ls[[i]] <- i_features
  }
  rl_mat <- t(do.call(cbind, rl_map_ls))
  rl_map <- data.frame(
    "gene_A" = as.character(rl_mat[,c("gene_A")]),
    "protein_A" = as.character(rl_mat[,c("protein_A")]),
    "type_A" = as.character(rl_mat[,c("type_A")]),
    "name_A" = as.character(rl_mat[,c("name_A")]),
    "gene_B" = as.character(rl_mat[,c("gene_B")]),
    "protein_B" = as.character(rl_mat[,c("protein_B")]),
    "type_B" = as.character(rl_mat[,c("type_B")]),
    "name_B" = as.character(rl_mat[,c("name_B")]),
    "int_pair" = as.character(rl_mat[,c("int_pair")]),
    "annotation_strategy" = as.character(rl_mat[,c("annotation_strategy")]),
    "source" = as.character(rl_mat[,c("source")]),
    "database_name" = as.character(rl_mat[,c("database_name")])
  )
  rownames(rl_map) <- seq(nrow(rl_map))
  return(rl_map)
}

#' Gene ortholog mapping
#' 
#' Given equal length vectors of genes from a reference genome and orthologous 
#' genes from another genome, converts a query gene to its ortholog in the
#' orthologous genome. Can handle cases of one-to-many mapping, one-to-none 
#' mapping, and many-to-one mapping. Returns a list where each value contains
#' a character vector containing all orthologs that map to the query gene in
#' the same position in the character vector query. "from" and "to" must be
#' equal length vectors where the reference genes share the same index in the
#' "from" vector as their orthologs in the "to" vector. Repeated values are 
#' permitted in both vectors
#' 
#' @param gene character. Vector of genes to be converted to gene orthologs
#' @param from character. Vector of genes from a reference database that contains the query genes
#' @param to character. Vector of gene orthologs from another reference genome to genes in from
#' @export
#' @return list
#' @examples
#' test_from <- c("HUMA1", "HUMB1", "HUMB1", "HUMC1", "HUMD1", "HUMD2")
#' test_to <- c("Musa1", "Musb1", "Musb2", "", "Musd1", "Musd1")
#' 
#' map_OneOne <- ortholog_mapping(gene = "HUMA1", from = test_from, to = test_to)
#' map_OneMany <- ortholog_mapping(gene = "HUMB1", from = test_from, to = test_to)
#' map_ManyOne <- ortholog_mapping(gene = c("HUMD1", "HUMD2"), from = test_from, to = test_to)
#' map_OneNone <- ortholog_mapping(gene = "HUMC1", from = test_from, to = test_to)
#' 
ortholog_mapping <- function(gene, from, to) {
  if(length(from) != length(to)) stop("'from' and 'to' vectors are of differing lengths. Ensure that all original genes in 'from' have an orthologous gene in 'to' at the same index of each vector.")
  ortho_ls <- lapply(gene, function(g) {
    if(g %in% from) return(to[from == g])
    else return(NA)
  })
  ortho_n <- vapply(ortho_ls, length, numeric(1))
  ortho <- list()
  if(0 %in% ortho_n) {ortho <- ""}
  if(max(ortho_n) > 1) {
    for(i in seq_along(gene)) {
      cap_g <- toupper(gene[i])
      cap_o <- toupper(ortho_ls[[i]])
      if(cap_g %in% cap_o) {
        ortho[[i]] <- ortho_ls[[i]][cap_o == cap_g]
      } else {
        ortho[[i]] <- ortho_ls[[i]]
      }
    }
  } else {
    ortho <- ortho_ls
  }
  return(ortho)
}

#' Receptor-Ligand Map Ortholog Conversion
#' 
#' Converts genes in an rl_map to orthologous genes in an alternative reference
#' genome. Partner names are maintained from the original rl_map passed to this
#' function. In cases where there are many orthologs to a single reference gene,
#' a seperate interaction is included for the new rl_map for each possible 
#' ortholog of the query gene. This can lead to interactions involving ligand
#' and receptor complexes to be expanded to many seperate rows. An option to
#' disregard interactions involving complexes is provided if users would like to
#' completely avoid spurious annotations of all gene orthologs participating in
#' the same protein complex as the converted reference gene.
#' 
#' @param rl_map data.frame. A receptor-ligand map. Must include columns for "gene_A" and "gene_B" at minimum.
#' @param from character. Vector of genes from a reference database that contains the query genes
#' @param to character. Vector of gene orthologs from another reference genome to genes in from
#' @param use_complexes logical. Whether to consider ligand-receptor interactions for gene conversions. If FALSE, the resulting rl_map will only include interactions of ligands and receptors each encoded by a single gene.
#' @param conversion_name character. A string to append to values in the "source" column of the rl_map to keep notes of how ortholog conversion was conducted.
#' @export
#' @return data.frame
#' @examples
#' test_rl_map <- data.frame(
#' gene_A = c(
#'   "HUMA1", "HUMB1", "HUMC1", "HUMD1", "HUMD2", 
#'   "HUMA1,HUMA2", "HUMB1,HUMA2", "HUMC1,HUMA2", "HUMD1,HUMA2", "HUMD2,HUMA2"
#'  ),
#' type_A = rep("R", 10),
#' name_A = c(
#'   "HUMA1", "HUMB1", "HUMC1", "HUMD1", "HUMD2", 
#'   "complexAA", "complexBA", "complexCA", "complexD1A", "complexD2A"
#'  ),
#' gene_B = rep("HUME", 10),
#' type_B = rep("L", 10),
#' name_B = rep("HUME", 10),
#' int_pair = c(
#'   "HUMA1 & HUME", "HUMB1 & HUME", "HUMC1 & HUME", "HUMD1 & HUME", "HUMD2 & HUME", 
#'   "complexAA & HUME", "complexBA & HUME", "complexCA & HUME", "complexD1A & HUME", 
#'   "complexD2A & HUME"
#'  ),
#' annotation = rep("test", 10),
#' source = c(
#'   "simple int, one-one ortho", "simple int, one-many ortho", "simple int, one-none ortho",
#'   "simple int, many-one ortho", "simple int, many-one ortho", "complex int, one-one ortho",
#'   "complex int, one-many ortho", "complex int, one-none ortho", "complex int, many-one ortho",
#'   "complex int, many-one ortho"
#'  ),
#' database_name = rep("test", 10)
#' )
#' test_from <- c("HUMA1", "HUMA2", "HUMB1", "HUMB1", "HUMC1", "HUMD1", "HUMD2", "HUME")
#' test_to <- c("Musa1", "Musa2", "Musb1", "Musb2", "", "Musd1", "Musd1", "Muse")
#' 
#' ortho_rl_map <- rl_map_ortholog_conversion(
#'   rl_map = test_rl_map, to = test_to, from = test_from, use_complexes = TRUE
#'  )
#' ortho_rl_map_noComplex <- rl_map_ortholog_conversion(
#'   rl_map = test_rl_map, to = test_to, from = test_from, use_complexes = FALSE
#'  )
#' 
rl_map_ortholog_conversion <- function(rl_map, from, to, use_complexes = FALSE, conversion_name = NULL) {
  ortho_map_ls <- list()
  ittr <- 1
  for(i in seq(nrow(rl_map))) {
    rl <- rl_map[i,]
    orig_A <- rl[["gene_A"]]
    orig_B <- rl[["gene_B"]]
    
    if( (grepl(",", orig_A) | grepl(",", orig_B)) & use_complexes == FALSE) next
    
    gene_A <- unlist(strsplit(orig_A, split = ","))
    ortho_A <- ortholog_mapping(gene = gene_A, from = from, to = to)
    if(sum(is.na(ortho_A)) | sum(ortho_A == "")) next
    gene_B <- unlist(strsplit(orig_B, split = ","))
    ortho_B <- ortholog_mapping(gene = gene_B, from = from, to = to)
    if(sum(is.na(ortho_B)) | sum(ortho_B == "")) next
    df_A <- expand.grid(ortho_A, stringsAsFactors = FALSE)
    df_B <- expand.grid(ortho_B, stringsAsFactors = FALSE)
    
    for(a in seq(nrow(df_A))) {
      replace_A <- paste(df_A[a,], collapse = ",")
      for(b in seq(nrow(df_B))) {
        replace_B <- paste(df_B[b,], collapse = ",")
        
        rl_ortho <- rl
        rl_ortho[["gene_A"]] <- replace_A
        rl_ortho[["gene_B"]] <- replace_B
        ortho_map_ls[[ittr]] <- rl_ortho
        ittr <- ittr + 1
      }
    }
  }
  rl_map_ortho <- do.call(rbind, ortho_map_ls)
  rl_map_ortho$database_name <- paste0(
    rl_map_ortho$database_name, "_conversion_", conversion_name
  )
  return(rl_map_ortho)
}
