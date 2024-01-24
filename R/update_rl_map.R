library(dplyr)
library(purrr)

gene_map_vec <- readRDS("../gene_ortholog_conversion/references/HumanMouseMap.rds")
gene_map_df <- data.frame(
  "HGNC" = gene_map_vec,
  "MGI" = names(gene_map_vec)
) |>
  dplyr::arrange(desc(MGI))

# occurances of one to many mapping
human_genes <- unique(gene_map_df$HGNC)
table_df <- gene_map_df$HGNC |>
  table() |>
  as.data.frame()
multi_ortho <- table_df[(table_df$Freq > 1) & !grepl("^$", table_df$Var1),]  



# path to cellphoneDB
cellphone_path <- "../CellPhoneDB_data/cellphonedb-data-4.0.0/data"
cpdb_files <- list.files(cellphone_path, full.names = TRUE)
cpdb_files <- cpdb_files[grep("_input.csv", cpdb_files)]
names(cpdb_files) <- sapply(cpdb_files, function(x) strsplit(x, "/") %>% unlist %>% rev %>% head(1) %>% strsplit("_") %>% unlist %>% head(1))

genes <- read_cellphoneDB_file(cpdb_files["gene"])
proteins <- read_cellphoneDB_file(cpdb_files["protein"], set_unannotated = TRUE)
interactions <- read_cellphoneDB_file(cpdb_files["interaction"])
complexes <- read_cellphoneDB_file(cpdb_files["complex"])

# construction of the rl_map will occur first, followed by gene conversion
# - this allows for rl_maps to be written from other references but they all
# can utilize this helper function.

# component functions for creation of map
read_cellphoneDB_file <- function(filename, set_unannotated = FALSE) {
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

# converts a character vector of "True" and "False" values to a logical vector of corresponding values.
# x: character vector. If all values are "True" or "False", the values will be converted to TRUE and FALSE, respectively.
# return: logical vector 
TF_syntax_to_R <- function(x) {
  if (identical(unique(x), c("True", "False")) | identical(unique(x), c("False","True"))) {
    y <- ifelse(x == "True", TRUE, FALSE)
  } else {
    y <- x
  }
  return(y)
}

check_filepath <- function(x, ...) {
  if(is.character(x)) {
    y <- read_cellphoneDB_file(x, ...)
  } else {
    y <- x
  }
  return(y)
}

parse_interaction <- function(int, genes, proteins, partner = c("A", "B"), complexes = NULL) {
  res <- list()
  p <- paste0("partner_", tolower(partner))
  pn <- paste0("protein_name_", tolower(partner))
  feat <- int[[p]]
  prot_name <- int[[pn]]
  if(feat %in% complexes[["complex_name"]]) {
    # annotation of protein complex
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
    res[[paste0("type_", partner)]] <- ifelse(prot[["receptor"]], "R", "L")
    res[[paste0("name_", partner)]] <- comp_gene
  } else {
    res[[paste0("gene_", partner)]] <- ""
    res[[paste0("type_", partner)]] <- ""
    res[[paste0("name_", partner)]] <- feat
  }
  return(res)
}


# create rl map

create_rl_map_cellphoneDB <- function(
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
    "type_A" = as.character(rl_mat[,c("type_A")]),
    "name_A" = as.character(rl_mat[,c("name_A")]),
    "gene_B" = as.character(rl_mat[,c("gene_B")]),
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

df <- create_rl_map_cellphoneDB(
  genes = genes, proteins = proteins, interactions = interactions,
  complexes = complexes
)
View(df)

# protein standin for non-protein ligand
parse_interaction(int = interactions[1,], genes = genes, proteins = proteins, partner = "A", complexes = complexes)
parse_interaction(int = interactions[1,], genes = genes, proteins = proteins, partner = "B", complexes = complexes)

# individual ligand receptor pair
parse_interaction(int = interactions[14,], genes = genes, proteins = proteins, partner = "A", complexes = NULL)
parse_interaction(int = interactions[14,], genes = genes, proteins = proteins, partner = "B", complexes = NULL)

# complex receptor
parse_interaction(int = interactions[17,], genes = genes, proteins = proteins, partner = "A", complexes = complexes)
parse_interaction(int = interactions[17,], genes = genes, proteins = proteins, partner = "B", complexes = complexes)

# complex ligand
parse_interaction(int = interactions[50,], genes = genes, proteins = proteins, partner = "A", complexes = complexes)
parse_interaction(int = interactions[50,], genes = genes, proteins = proteins, partner = "B", complexes = complexes)

# ortholog conversion will be based on an approach similar to map values
# - from vector (original gene names for conversion)
# - to vector (orthologous genes to replace the prior value)
# in a one to many match, exact gene name match is used
# - if there are no exact matches, each possible ortholog is used in the conversion

# component functions
ortholog_mapping <- function(gene, from, to) {
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

# rl map ortholog conversion

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



df2_c <- rl_map_ortholog_conversion(rl_map = df, from = gene_map_df$HGNC, to = gene_map_df$MGI, use_complexes = TRUE)
df2_o <- rl_map_ortholog_conversion(rl_map = df, from = gene_map_df$HGNC, to = gene_map_df$MGI, use_complexes = FALSE)

table_convert <- function(gene, from, to){
  res <- vapply(
    gene,
    FUN = function(g){
      i <- which(from == g)
      ortho <- to[i]
      if(length(ortho) > 1) {
        ortho <- ortho[1]
      }
      return(ortho)
    },
    FUN.VALUE = character(1)
  )
}



map_interaction <- function(interacting_partner, map_ortholog=NULL) {
  features <- NULL
  if (interacting_partner %in% proteins[["uniprot"]]) {
    protein <- proteins[proteins[["uniprot"]] == interacting_partner, ]
    gene <- unique(genes[genes[["uniprot"]] == protein[["uniprot"]], c("gene_name")])  # Sometimes we get multiple genenames for one uniprot: Q7Z5A7. This code will create a separate entry for each gene 
    if (!is.null(map_ortholog)) {
      orthologs <- mapOrtholog(gene, map_ortholog)
      if (is.null(orthologs)) {
        features <- NULL
        return(features)
      } else {
        #Handle one-many human-mouse mapping
        #Check if one of the mouse orthologs has the same name - keep only that. Else keep all
        matching_name <- str_to_title(gene)
        if (matching_name %in% orthologs) {
          gene <- matching_name
        } else {
          gene <- orthologs
        }
      }
    }
    features <- data.frame(uniprot = protein[["uniprot"]], gene = gene, type = ifelse(protein[["receptor"]], "R", "L"), name = gene)
  }
  return(features)
}

mapOrtholog <- function(gene, orthoMap) {
  if (sum(gene %in% orthoMap[,1]) < length(gene)) {
    ortho <- NULL
  } else {
    ortho <- unique(orthoMap[orthoMap[, 1] %in% gene, 2])
  }
  return(ortho)
}


