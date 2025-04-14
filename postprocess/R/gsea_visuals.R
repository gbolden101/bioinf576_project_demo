#' arguments
#' Required columns: GSEA pvalue, GSEA Term, GSEA Count column
#' organize by Count or pvalue GSEA results (org = 'Count', org = 'pvalue')
#' num_disp
#' GSEA pvalue filter
#' graph title (title)
#' file_name write out
#' output visual dotp


#' @export
#' @import tidyverse
#' @import dplyr
#' @import RColorBrewer
#' @import ggrepel

gsea_dotplot <- function(gsea_results, gspval_cutoff, org_by, num_disp, graph_title, file_name){
  if(org_by %in% 'Count'){
    #inside if
    dotplot <- gsea_results%>%
      arrange(desc(Count))%>%
      head(num_disp)%>%
      mutate(as.factor(Description))%>%
      ggplot() +
      aes(y = Count, x = fct_reorder(Description, Count)) +
      scale_colour_gradient(low = "red", high = "blue") +
      geom_point(aes(size = GeneRatio, color = pvalue)) +
      labs(title = graph_title, x = "Enriched Term", y = "Gene Count") +
      coord_flip() +
      theme(text = element_text(face = "bold") , plot.title = element_text(hjust = 1))

  }
  if(org_by %in% 'pvalue'){
    #inside if
    dotplot <- gsea_results %>%
      filter(pvalue < gspval_cutoff) %>%
      head(num_disp) %>%
      mutate(as.factor(Description)) %>%
      mutate(pvalue = -1 * log(pvalue)) %>%
      ggplot() +
      aes(x = pvalue, y = fct_reorder(Description, pvalue))+
      scale_colour_gradient(low = "red", high = "blue") +
      geom_point(aes(size = GeneRatio, color = Count)) +
      labs(title = graph_title, x = "-log10(P)", y = "Enriched Term") +
      theme(text = element_text(face = "bold"), plot.title = element_text(hjust = 1))


  }
  return(list(dotplot, ggsave(file_name, device = "png", width = 8, height = 6, units = "in")))
}


#'arguments: gsea_results,
#'gsea_result desc: should have following columns: a gene_id, gsea_pvalue, gsea_term,gsea_count


#' @export
#' @import tidyverse
#' @import dplyr
#' @import RColorBrewer
#' @import ggrepel


gsea_volplot <- function(gsea_results, unique = FALSE, lower_xlim = 10, upper_xlim = 10, step = 1){
  #setting the color

  gsea_results$color <- ifelse(gsea_results$log2FoldChange > 0, "Upregulated",  "Downregulated")

  #arrange by padj, may give that option to the user
  #the idea is that we are using DAVID so that means we some strict
  #It looks like this function is under the assumption that annotation
  gsea_results <- gsea_results %>%
    arrange(padj) %>%
    mutate(delabel = ifelse(row_number() <= 10, as.character(external_gene_name), ""))

  if(unique) {
    gsea_results <- gsea_results %>% distinct(external_gene_name, .keep_all = TRUE)
  }

  volplot <- ggplot(gsea_results) +
    aes(x = log2FoldChange, y = -log10(pvalue), color = color, label = delabel) +


    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +

    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +

    theme_set(theme_classic(base_size = 12) +

                theme(
                  axis.title.y = element_text(margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),

                  axis.title.x = element_text(hjust = 0.5, margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),

                  plot.title = element_text(hjust = 0.5)
                )) +

    geom_point(size = 2) +

    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +

    geom_text(nudge_y = 0.5, check_overlap = TRUE) +

    geom_text_repel(max.overlaps = Inf) +
    coord_cartesian(xlim = c(lower_xlim, upper_xlim)) +
    scale_x_continuous(breaks = seq(lower_xlim, upper_xlim, step))

  return(list(volplot, gsea_results))
}
