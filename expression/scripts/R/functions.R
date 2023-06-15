#######################################
## Quality Control Functions ##
#######################################

## This function plots a PCA and colors it by a specified metadata column
## It takes in a rld object and a character string of the metadata column to color by, and outputs a PCA plot 
## If label = TRUE, all points are labelled.

make_pca <- function(rld, intgroup,  
                     title = "title",         
                     xlimits = c(-15, 15), 
                     ylimits = c(-10, 10),
                     label = FALSE){
  # Calculations
  ntop = 500
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  data <- plotPCA(rld, intgroup = intgroup, 
                  returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"), digits = 2)
  
  # Plot 
  pca_plot <- ggplot(as.data.frame(data), aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = as.data.frame(data)[[intgroup]]), 
               size = 2, alpha = 0.8) +
    labs(title = title,
         subtitle = paste0("By ", intgroup), 
         colour = intgroup) +
    
    # Add scale annotations
    scale_x_continuous(paste0("PC1: ",percentVar[1],"% variance"), 
                       limits = xlimits) +
    scale_y_continuous(paste0("PC2: ",percentVar[2],"% variance"), 
                       limits = ylimits) +
    scale_color_brewer(palette = "Paired") + 
    
    # Make 1 unit on x-axis equal to 1 unit on y-axis
    coord_fixed(ratio = 1) +
    theme_classic()
  
  if(label == TRUE){
    pca_plot <- pca_plot + 
      geom_text_repel(data = data, aes(PC1,PC2, label = name), 
                      hjust = 0.5, box.padding = 0.5, size = 3,
                      max.overlaps = Inf)
    return(pca_plot)
  } else {
    return(pca_plot)
  }
  
}

## This function plots PCA, colors by a metadata column and allows shape by another column
## It takes in a rld object, a character string of the metadata column to color by,
## another character string of the metadata to shape by. Outputs a PCA plot. 
## If label = TRUE, all points are labelled.

make_pca_shape <- function(rld, color_by, shape_by, 
                           title = "title",   
                           subtitle = element_blank(), 
                           xlimits = c(-15, 15), 
                           ylimits = c(-10, 10),
                           label = FALSE){
  # Calculations
  ntop = 500
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  data <- plotPCA(rld, intgroup = c(color_by, shape_by),
                  returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"), digits = 2)
  
  # Plot 
  pca_plot <- ggplot(as.data.frame(data), aes(x = PC1, y = PC2)) + 
    geom_point(aes(
      color = as.data.frame(data)[[color_by]],
      shape = as.data.frame(data)[[shape_by]]
    ), 
    size = 2, alpha = 0.8) +
    labs(title = title,
         subtitle = subtitle, 
         colour = color_by, 
         shape = shape_by) +
    
    # Add scale annotations
    scale_x_continuous(paste0("PC1: ",percentVar[1],"% variance"), 
                       limits = xlimits) +
    scale_y_continuous(paste0("PC2: ",percentVar[2],"% variance"), 
                       limits = ylimits) +
    scale_color_brewer(palette = "Set1") + 
    
    # Make 1 unit on x-axis equal to 1 unit on y-axis
    coord_fixed(ratio = 1) +
    theme_classic()
  
  if(label == TRUE){
    pca_plot <- pca_plot + 
      geom_text_repel(data = data, aes(PC1,PC2, label = name), 
                      hjust = 0.5, box.padding = 0.5, size = 3,
                      max.overlaps = Inf)
    return(pca_plot)
  } else {
    return(pca_plot)
  }
  
}

############################################
## Differential Expression Functions ##
############################################


## This function gets a list of lfcShrink fold changes from a DESeq object after nbinomWaldTest has been called.
## It takes in an nbinomWaldTest object, a character string contrast, and an ensembl.genes object containing gene biotype for annotations. 

get_dds_res <- function(wald_dds, contrast, ensembl.genes, 
                        lfcshrinktype = "apeglm", 
                        shrink = TRUE, parallel = TRUE) {
  
  # Run: get_dds_res(wald_dds, 
  #                  contrast = c("condition", "treated", "ctrl"),
  #                  ensembl.genes = ensembl.genes,
  #                  shrink = TRUE)
  
  # Get the results
  dds_res = results(wald_dds, 
                contrast = contrast,  
                filter = rowMeans(counts(wald_dds, 
                normalized = TRUE)), 
                test = "Wald", alpha = 0.1, 
                independentFiltering = TRUE)
  
  # Make the condition_mdd_vs_ctrl string
  deseq_coef <- paste(contrast[1], contrast[2], "vs", 
                      contrast[3], sep = "_")
  print(deseq_coef)
  
  # Shrink
  if (shrink == TRUE){
    dds_res <- lfcShrink(wald_dds, 
                     coef = deseq_coef, 
                     res = dds_res, 
                     type = lfcshrinktype, parallel = TRUE)
  }
  
  # Add gene annotations
  dds_res$gene_biotype = ensembl.genes$gene_biotype[match(row.names(dds_res), ensembl.genes$gene_id)]
  dds_res$external_gene_name = ensembl.genes$external_gene_name[match(row.names(dds_res), ensembl.genes$gene_id)]
  
  return(dds_res)
}



## This function adds log2fc and pvalues to results(dds), and returns a tibble

add_res_info <- function(dds_res){
  # Takes in a results(dds), converts to a tibble, adds log2fc and pvalue info 
  
  df_export <- as_tibble(dds_res, rownames = "gene_name")
  
  # Add pvalue and log2fc
  df_export$log2fc_info <- dds_res@elementMetadata@listData$description[2]
  df_export$pval_info <- dds_res@elementMetadata@listData$description[5]
  
  return(df_export)
}
