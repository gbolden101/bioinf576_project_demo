#' @import biomaRt
#' @import dplyr
#' @import tidyverse
# common example datasets: "hsapiens_gene_ensembl", "mmusculus_gene_ensembl"
#connects to the ensembledb, grabs and returns
#arguments: an entrez id vector, and the name of ensembldataset
#It would be wise to load in your chosen OrgDb set before running this code

#' @export
entrezid_to_genename <-function(entrez_vector,
                                ensembl_dataset = "hsapiens_gene_ensembl"){
  #connect to the Ensembl database for Mmu
  ensembl = useMart("ensembl", dataset = ensembl_dataset)
  #getting the corresponding Entrez Gene IDs
  entrezid_genename_df <- getBM(attributes = c("entrezgene_id", "external_gene_name"), filters = "entrezgene_id", values = entrez_vector, mart = ensembl)
  return(entrezid_genename_df)
}

#' enrichgo_unpack description
#' unpacks the genes based on the specified Description term from the enrichGO results an return all the genes under the the term
#' along side their deseq2 results. This will prepare them for volcano plot visualization
#' arguments
#' deseq2_table: processed data with the following columns: entrezgene_id,log2FoldChange, pvalue, padjust (high reccommend these four columns be present)
#' enrichGO reslts from ClusterProfiler package, the the chosen description (chosen_desc), and the ensembl_dataset that will passed to the above fxn
#' This will prep your data for the Volcano plot visual function and you can export and report the genes
#' BIG SIDE NOTE: this will only pre your data for the volcano plot fxn and not thed dotplot, visual function

#' example inputs
#' enrichgo_result <- key@result
#' chosen_desc <- "stress response to metal ion"
#' ensemble_dataset <- "hsapiens_gene_ensembl"
#' deseq2_table <- object1
#' enrichgo_unpack(deseq2_table, enrichgo_result, chosen_desc, ensemble_dataset)

#' @export
enrichgo_unpack <- function(deseq2_table, enrichgo_result, chosen_desc, ensembl_dataset){
  ###Extracting a specific gene into a cluster

  cluster <- tibble(key@result) %>%
    filter(str_detect(Description, chosen_desc)) %>%
    mutate(Description = as.character(Description))
  cluster

  ###creates a helper function to unlist and strsplits the enriched Gene IDs per term

  convert_to_integer <- function(str) {
    bruh <- as.integer(unlist(strsplit(str, "/")))
    return(bruh)
  }

  #apply this function to the geneID column
  extract_result <- lapply(cluster$geneID, convert_to_integer)
  #Next check the length of the new vector to see if it matches the Count column in the GSEA DF.
  length(extract_result[[1]])

  #create an empty vector
  path_vector <- c()
  #for each index in the first column of numbers
  for (i in seq_along(cluster$Description)) {
    #replicate the name of the Description for the length of the GO term in column 1
    path_vector <- c(path_vector, rep(cluster$Description[i], length(extract_result[[i]])))
  }


  # make new tibble that pairs the name of the Pathway from the GOTerm and enriched Entrez ID from the "post-processed" GeneID column
  pathway_tibble <- tibble(gene_desc = path_vector, entrezgene_id = unlist(extract_result))
  pathway_tibble

  pathway_genenames <- entrezid_to_genename(pathway_tibble$entrezgene_id, ensembl_dataset = ensemble_dataset)



  final_path_tibble <- inner_join(pathway_tibble, pathway_genenames, by = "entrezgene_id")
  final_path_tibble <- final_path_tibble %>% mutate(entrezgene_id = as.character(entrezgene_id))
  final_path_tibble

  if (!"external_gene_name" %in% colnames(deseq2_table)) {
    colnames(deseq2_table)[1] <- "external_gene_name"
  }

  final_path_fc_tibble <- inner_join(final_path_tibble, deseq2_table, by = "external_gene_name")
  final_path_fc_tibble
  return(final_path_fc_tibble)
}
