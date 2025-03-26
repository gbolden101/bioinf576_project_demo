#' deseq filtration uses common filtration methods applied to post-processed RNA-Seq data
#' The goal is make my life easier by streamlining this process
#' Added bones each step is written out to excel spread sheet so that you can see your changes happen and show others
#' deseq_filtration will return a data.fram of the results and the name of the file you've chosen


#' @export
#' @import tidyverse
#' @import dplyr
#' @import openxlsx

deseq_filteration <- function(tibble, 
                              file_name, 
                              pval_cutoff, 
                              rev = FALSE){
  #added argument: pval_cutoff, allows user to set their pvalue cutoff
  # This function takes into consideration that the Deseq2 Results where flipped
  # I'll have to make a conditional on that front
  #
  wbfxn <- createWorkbook()
  
  addWorksheet(wbfxn, paste("Initial Data", nrow(tibble)))
  writeData(wbfxn, sheet = paste("Initial Data", nrow(tibble)) , tibble)
  
  tibble_exp <- tibble %>%
    # filter in value <=0.05
    
    #In the past I have encountered Deseq sets wherelog2FCs where the results are flipped
    #This conditional takes that into consideration
    {
      if (rev) {
        . <- mutate(., INUN_log2FC = log2FoldChange * -1)
        . <- relocate(., INUN_log2FC, .after = 2)
        addWorksheet(wbfxn, paste("Corrected Log2FC", nrow(.)))
        writeData(wbfxn, sheet = paste("Corrected Log2FC", nrow(.)), .)
      } else {
        . <- mutate(., log2FoldChange = log2FoldChange)
      }
      
      .}%>%
    filter(pvalue < pval_cutoff) %>%
    {
      addWorksheet(wbfxn, paste("pvalue <=", pval_cutoff, nrow(.)))
      writeData(wbfxn, sheet = paste("pvalue <=", pval_cutoff, nrow(.)), .)
      
      .}%>%
    
    {
      if ("INUN_log2FC" %in% colnames(.)) {
        mutate(., plusFC = 2^INUN_log2FC)
      } else {
        mutate(., plusFC = 2^log2FoldChange)
      }
    } %>%
    
    relocate(plusFC, .before = 4) %>%
    {
      
      addWorksheet(wbfxn, "add plusFC") 
      writeData(wbfxn,sheet = "add plusFC", .) 
      
      .}%>%
    
    arrange(desc(plusFC)) %>%
    {
      addWorksheet(wbfxn, "Sort by plusFC desc.") 
      writeData(wbfxn, sheet = "Sort by plusFC desc.", .)
      
      .}%>%
    
    #add a column called negFC: if a value in plusFC is less than 0.05 then convert it to
    #negFC (-1/FC) else leave the value alone by multiplying it by 1
    mutate(negFC = if_else(plusFC < 0.5, -1/plusFC, plusFC * 1)) %>%
    relocate(negFC, .before = 5) %>%
    {
      addWorksheet(wbfxn, "add negFC") 
      writeData(wbfxn, sheet = "add negFC", .)
      .}%>%
    
    #In a new column called FC, 
    #If the values in plusFC and negativeFC match (plusFC == negFC) then copy the values that do match into the new column and convert them to the character type 
    #Else: paste both values from negFC and plusFC into the same point and seperate them by a "/" character.
    mutate(FC = if_else(plusFC == negFC, as.character(plusFC), paste(plusFC, negFC, sep = "/"))) %>%
    #In the new column FinalFC, If an entry has a slash charcter in the column FC, then get rid of everycharacter before the slash in each entry including the slash.
    mutate(FinalFC = gsub(".*?/","", FC)) %>%
    #convert the column FinalFC into a numeric type variable
    #The goal of the FinalFC column is only keep the strongest directional fold change
    mutate(FinalFC = as.numeric(FinalFC)) %>%
    {
      addWorksheet(wbfxn, "add FinalFC") 
      writeData(wbfxn,sheet = "add FinalFC", .)
      
      .}%>%
    # exlude the following columns: plusFC, negFC, and FC
    dplyr::select(-plusFC, -negFC, -FC) %>%
    relocate(FinalFC, .before = 4) %>%
    {
      addWorksheet(wbfxn, "exclude plusFC, negFC, intFC") 
      writeData(wbfxn, sheet = "exclude plusFC, negFC, intFC", .)
      
      
      .}%>%
    #keep genes if the FinalFC is greater than +2 or less than -2
    filter(FinalFC > 2 | FinalFC < -2) %>%
    #keep genes if the log2foldchange is greater than +1 or less than -1
    filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
    {
      addWorksheet(wbfxn, "FinalFC > 2, <-2") 
      writeData(wbfxn, sheet = "FinalFC > 2, <-2", .)
      
      .}
  tibble_exp
  file_path <- file.path(getwd(), paste0(file_name, ".xlsx"))
  
  worksheetOrder(wbfxn) <- c(1:length(wbfxn$worksheets))
  saveWorkbook(wbfxn, file_path, overwrite = TRUE)
  return(list(tibble_exp = tibble_exp, file_name = file_name))
}

#example use
#results frame being the
#object <- deseq_filteration(tibble =  results_frame, file_name = "test1_2", pval_cutoff = 0.01, rev = FALSE)

#'df_to_txt was my first function I made by myself it gets a special seat.

#' @export
#' @import tidyverse
df_to_txt <- function(tibble) {
  
  tibble_name <- deparse(substitute(tibble))
  
  file_name <- paste(getwd(),"/",tibble_name,".txt")
  
  file_mod <- gsub("\\s", "", file_name)
  #\\s is a regular expression pattern that matches any white space character
  #gsub replaces all occurance of this pattern with an empty string
  save_fxn <- write_delim(tibble, file = file_mod, delim = "\t", col_names = FALSE)
  
  return(save_fxn)
}




