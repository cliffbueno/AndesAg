# Function for quick but high quality taxonomic barplots
# Procedure:
# Take an mctoolsr object
# 1. Get legend title based on taxonomic level
# 2. Summarize by taxonomic level
# 3. Extract data by factor of the top n taxa
# 4. Update data frame for NA and unclassified taxa
# 5. Make ggplot

cliffplot_taxa_bars <- function(input, level, variable) {
  
  if (level == 1) {
    leg_title <- "Domain"
  }
  
  if (level == 2) {
    leg_title <- "Phylum"
  }
  
  if (level == 3) {
    leg_title <- "Class"
  }
  
  if (level == 4) {
    leg_title <- "Order"
  }
  
  if (level == 5) {
    leg_title <- "Family"
  }
  
  if (level == 6) {
    leg_title <- "Genus"
  }
  
  if (level == 7) {
    leg_title <- "Species"
  }
  
  if (level == 8) {
    leg_title <- "ASV"
  }
  
  if (level == 9) {
    leg_title <- "Guild"
  }
  
  tax_sum <- summarize_taxonomy(input = input, 
                                level = level, 
                                report_higher_tax = F)

  if (level == 1) {  
  bars <- plot_taxa_bars(tax_sum,
                         input$map_loaded,
                         variable,
                         num_taxa = 12,
                         data_only = TRUE)
  cliffplot <- ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
    geom_bar(stat = "identity", colour = NA, size = 0.25) +
    labs(x = "Estuary", y = "Relative abundance", fill = leg_title) +
    scale_fill_manual(values = c("red", "blue")) +
    scale_y_continuous(expand = c(0.01, 0.01)) + 
    theme_classic() +
    theme(axis.title.y = element_text(face = "bold", size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  }
  
  if (level == 2) {
    bars <- plot_taxa_bars(tax_sum,
                           input$map_loaded,
                           variable,
                           num_taxa = 12,
                           data_only = TRUE) %>%
      mutate(taxon = fct_rev(taxon))
    cliffplot <- ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
      geom_bar(stat = "identity", colour = NA, size = 0.25) +
      labs(x = "Estuary", y = "Relative abundance", fill = leg_title) +
      scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
      scale_y_continuous(expand = c(0.01, 0.01)) + 
      theme_classic() +
      theme(axis.title.y = element_text(face = "bold", size = 12),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  }
  
  if (level > 2) {
    bars <- plot_taxa_bars(tax_sum,
                           input$map_loaded,
                           variable,
                           num_taxa = 12, 
                           data_only = TRUE) %>%
      mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
      mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
      mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
      mutate(taxon = fct_rev(taxon))
    cliffplot <- ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
      geom_bar(stat = "identity", colour = NA, size = 0.25) +
      labs(x = "Estuary", y = "Relative abundance", fill = leg_title) +
      scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
      scale_y_continuous(expand = c(0.01, 0.01)) + 
      theme_classic() +
      theme(axis.title.y = element_text(face = "bold", size = 12),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  }
  
  cliffplot
  
}
