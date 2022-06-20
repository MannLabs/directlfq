
#ion_file <- "Y:/User/Constantin/iq_directlfq_comparison/iq_formatted_ions_bruderer_15.tsv"
protein_file <- "Y:/User/Constantin/iq_directlfq_comparison/iq_formatted_proteins_bruderer_15.tsv"
quant_file <- "Y:/User/Constantin/iq_directlfq_comparison/iq_formatted_quant_table_bruderer_15.tsv.scaled_to_10000.tsv"
sample_file <- "Y:/User/Constantin/iq_directlfq_comparison/iq_formatted_sample_bruderer_15.tsv.scaled_to_10000.tsv"


quant_df <- read.csv(quant_file, sep = '\t')
sample_df <- read.csv(sample_file, sep = '\t')
protein_df <- read.csv(protein_file, sep = '\t')

protein_list <- quant_df$protein_list
sample_list <- quant_df$sample_list
id <- quant_df$id
quant <- quant_df$quant


convert_to_array <- function(input_list) {
  return(array(input_list, c(length(input_list), 1)))
}
quant_df2 = NULL

quant_df2$protein_list <-convert_to_array(protein_list)
quant_df2$sample_list <- convert_to_array(sample_list)
quant_df2$id <- convert_to_array(id)
quant_df2$quant <- convert_to_array(quant_df$quant)


system.time({
  iq_norm_data <- iq::fast_preprocess(quant_df2, pdf_out = NULL)
  
  result_fastest <- iq::fast_MaxLFQ(iq_norm_data, 
                                    row_names = protein_df$PG.ProteinGroups, 
                                    col_names = sample_df$sample_list)
})



write.table(cbind(protein = rownames(result_fastest$estimate),
                  MaxLFQ_annotation = result_fastest$annotation,
                  result_fastest$estimate), 
            "iq-MaxLFQ-fast.txt", sep = "\t", row.names = FALSE)
















