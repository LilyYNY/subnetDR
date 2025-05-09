#' @title Perform GO Enrichment Analysis
#' @description Conducts GO enrichment analysis using MSigDB gene sets for specified ontology.
#' @importFrom msigdbr msigdbr
#' @importFrom magrittr %>%
#' @importFrom clusterProfiler enricher
#' @importFrom tidyr separate
#' @importFrom ggplot2 element_blank
#' @param gene Character vector of gene symbols for analysis.
#' @param Character string specifying GO ontology ("BP", "CC", "MF", or "all"; default: "all").
#' @param Numeric threshold for p-value (default: 0.05).
#' @param Character string specifying p-value adjustment method (default: "BH").
#' @param Numeric threshold for q-value (default: 0.2).
#' @param universe  Character vector of background genes (default: NULL).
#' @param minGSSize Minimum gene set size (default: 10).
#' @param maxGSSize Maximum gene set size (default: 500).
#' @return Data frame of enrichment results or NULL if no significant results.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @rdname enrich.GO
#' @export
#' @importFrom dplyr select
enrich.GO <- function(gene, ont = "all",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                      universe = NULL,
                      minGSSize = 10,
                      maxGSSize = 500) {
  # Validate ontology input
  if (ont == "all") {
    BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
      dplyr::select(gs_name, gene_symbol)
    CC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>%
      dplyr::select(gs_name, gene_symbol)
    MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
      dplyr::select(gs_name, gene_symbol)
    GO <- rbind(BP, CC, MF)  # Combine all GO ontologies
  } else {
    valid_ont <- c("BP", "CC", "MF")
    if (!ont %in% valid_ont) {
      stop("Invalid ont value: ", ont, ". Must be one of: ", paste(valid_ont, collapse = ", "))
    }
    GO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = paste0("GO:", ont)) %>%
      dplyr::select(gs_name, gene_symbol)
  }

  # Perform enrichment analysis
  enrich <- enricher(gene = gene,
                     TERM2GENE = GO,
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff = pvalueCutoff,
                     qvalueCutoff = qvalueCutoff,
                     universe = universe,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize)

  # Process results if enrichment is successful
  if (!is.null(enrich)) {
    enrichRes <- enrich@result
    enrichRes <- enrichRes[enrichRes$p.adjust < pvalueCutoff, ]
    if (nrow(enrichRes) > 0) {
      # Separate Description into Category and Description
      enrichRes <- separate(enrichRes, col = "Description", into = c("Category", "Description"), sep = "_", extra = "merge")
      enrichRes$Category <- gsub("GO", "", enrichRes$Category)  # Clean Category names
      enrichRes <- enrichRes[order(enrichRes$Category, enrichRes$pvalue), ]  # Sort by Category and p-value
      # Convert GeneRatio to numeric
      GeneRatio <- data.frame(enrichRes$GeneRatio)
      enrichRes$geneRatio <- apply(GeneRatio, 1, function(x) eval(parse(text = x)))
      # Compute RichFactor and FoldEnrichment
      enrichRes$RichFactor <- (as.numeric(sapply(strsplit(enrichRes$GeneRatio, "/"), `[`, 1)) /
                                 as.numeric(sapply(strsplit(enrichRes$BgRatio, "/"), `[`, 1)))
      enrichRes$FoldEnrichment <- enrichRes$RichFactor /
        (as.numeric(sapply(strsplit(enrichRes$GeneRatio, "/"), `[`, 2)) /
           as.numeric(sapply(strsplit(enrichRes$BgRatio, "/"), `[`, 2)))
      enrichRes$zScore <- NA  # Placeholder for z-score
      return(enrichRes)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' @title Perform Pathway Enrichment Analysis
#' @description Conducts pathway enrichment analysis using MSigDB gene sets for specified pathway database.
#' @importFrom msigdbr         msigdbr
#' @importFrom magrittr        %>%
#' @importFrom dplyr            select
#' @importFrom clusterProfiler  enricher
#' @importFrom tidyr            separate
#' @importFrom ggplot2 element_blank
#' @param gene  Character vector of gene symbols for analysis.
#' @param CP Character string specifying pathway database ("KEGG", "BIOCARTA", "PID", "REACTOME", "WIKIPATHWAYS", or "all"; default: "KEGG").
#' @param pvalueCutoff Numeric threshold for p-value (default: 0.05).
#' @param pAdjustMethod Character string specifying p-value adjustment method (default: "BH").
#' @param qvalueCutoff Numeric threshold for q-value (default: 0.2).
#' @param universe Character vector of background genes (default: NULL).
#' @param minGSSize Minimum gene set size (default: 10).
#' @param maxGSSize Maximum gene set size (default: 500).
#' @return Data frame of enrichment results or NULL if no significant results.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @rdname enrich.Pathway
#' @export
#' @importFrom dplyr select
enrich.Pathway <- function(gene, CP = "KEGG",
                           pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                           universe = NULL,
                           minGSSize = 10,
                           maxGSSize = 500) {
  # Load pathway gene sets from MSigDB
  if (CP == "all") {
    KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
      dplyr::select(gs_name, gene_symbol)
    Biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") %>%
      dplyr::select(gs_name, gene_symbol)
    PID <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID") %>%
      dplyr::select(gs_name, gene_symbol)
    Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
      dplyr::select(gs_name, gene_symbol)
    WikiPathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
      dplyr::select(gs_name, gene_symbol)
    Pathway <- rbind(KEGG, Biocarta, PID, Reactome, WikiPathways)  # Combine all pathways
  } else {
    # Map pathway input to MSigDB subcategory
    subcollection_map <- list(
      KEGG = "KEGG_LEGACY",
      BIOCARTA = "BIOCARTA_LEGACY",
      PID = "PID",
      REACTOME = "REACTOME",
      WIKIPATHWAYS = "WIKIPATHWAYS"
    )
    subcollection <- subcollection_map[[toupper(CP)]]
    if (is.null(subcollection)) {
      stop("Invalid CP value: ", CP)
    }
    Pathway <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = subcollection) %>%
      dplyr::select(gs_name, gene_symbol)
  }

  # Perform enrichment analysis
  enrich <- enricher(gene = gene,
                     TERM2GENE = Pathway,
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff = pvalueCutoff,
                     qvalueCutoff = qvalueCutoff,
                     universe = universe,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize)

  # Process results if enrichment is successful
  if (!is.null(enrich)) {
    enrichRes <- enrich@result
    enrichRes <- enrichRes[enrichRes$p.adjust < pvalueCutoff, ]
    if (nrow(enrichRes) > 0) {
      # Separate Description into Category and Description
      enrichRes <- separate(enrichRes, col = "Description", into = c("Category", "Description"), sep = "_", extra = "merge")
      enrichRes <- enrichRes[order(enrichRes$Category, enrichRes$pvalue), ]  # Sort by Category and p-value
      # Convert GeneRatio to numeric
      GeneRatio <- data.frame(enrichRes$GeneRatio)
      enrichRes$geneRatio <- apply(GeneRatio, 1, function(x) eval(parse(text = x)))
      # Compute RichFactor and FoldEnrichment
      enrichRes$RichFactor <- (as.numeric(sapply(strsplit(enrichRes$GeneRatio, "/"), `[`, 1)) /
                                 as.numeric(sapply(strsplit(enrichRes$BgRatio, "/"), `[`, 1)))
      enrichRes$FoldEnrichment <- enrichRes$RichFactor /
        (as.numeric(sapply(strsplit(enrichRes$GeneRatio, "/"), `[`, 2)) /
           as.numeric(sapply(strsplit(enrichRes$BgRatio, "/"), `[`, 2)))
      enrichRes$zScore <- NA  # Placeholder for z-score
      return(enrichRes)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' @title Perform Hallmark Enrichment Analysis
#' @description Conducts Hallmark enrichment analysis using MSigDB gene sets.
#' @importFrom msigdbr         msigdbr
#' @importFrom magrittr        %>%
#' @importFrom dplyr           select
#' @importFrom clusterProfiler enricher
#' @importFrom tidyr           separate
#' @importFrom ggplot2 element_blank
#' @param gene Character vector of gene symbols for analysis.
#' @param pvalueCutoff Numeric threshold for p-value (default: 0.05).
#' @param pAdjustMethod Character string specifying p-value adjustment method (default: "BH")
#' @param qvalueCutoff Numeric threshold for q-value (default: 0.2).
#' @param universe Character vector of background genes (default: NULL).
#' @param minGSSize Minimum gene set size (default: 10).
#' @param maxGSSize Maximum gene set size (default: 500).
#' @return Data frame of enrichment results or NULL if no significant results.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @rdname enrich.Hallmark
#' @export
#' @importFrom dplyr select
enrich.Hallmark <- function(gene,
                            pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
                            universe = NULL,
                            minGSSize = 10,
                            maxGSSize = 500) {
  # Load Hallmark gene sets from MSigDB
  Hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol)

  # Perform enrichment analysis
  enrich <- enricher(gene = gene,
                     TERM2GENE = Hallmark,
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff = pvalueCutoff,
                     qvalueCutoff = qvalueCutoff,
                     universe = universe,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize)

  # Process results if enrichment is successful
  if (!is.null(enrich)) {
    enrichRes <- enrich@result
    enrichRes <- enrichRes[enrichRes$p.adjust < pvalueCutoff, ]
    if (nrow(enrichRes) > 0) {
      # Separate Description into Category and Description
      enrichRes <- separate(enrichRes, col = "Description", into = c("Category", "Description"), sep = "_", extra = "merge")
      enrichRes <- enrichRes[order(enrichRes$Category, enrichRes$pvalue), ]  # Sort by Category and p-value
      # Convert GeneRatio to numeric
      GeneRatio <- data.frame(enrichRes$GeneRatio)
      enrichRes$geneRatio <- apply(GeneRatio, 1, function(x) eval(parse(text = x)))
      # Compute RichFactor and FoldEnrichment
      enrichRes$RichFactor <- (as.numeric(sapply(strsplit(enrichRes$GeneRatio, "/"), `[`, 1)) /
                                 as.numeric(sapply(strsplit(enrichRes$BgRatio, "/"), `[`, 1)))
      enrichRes$FoldEnrichment <- enrichRes$RichFactor /
        (as.numeric(sapply(strsplit(enrichRes$GeneRatio, "/"), `[`, 2)) /
           as.numeric(sapply(strsplit(enrichRes$BgRatio, "/"), `[`, 2)))
      enrichRes$zScore <- NA  # Placeholder for z-score
      return(enrichRes)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' @title Plot Enrichment Results
#' @description Generates a dot plot for enrichment results and saves it as a PDF.
#' @importFrom magrittr        %>%
#' @importFrom dplyr           filter if_all everything group_by slice_min ungroup select bind_rows all_of
#' @importFrom ggplot2         ggplot aes geom_point scale_color_gradient xlab ylab theme_bw theme element_text guides guide_legend guide_colorbar ggsave
#' @importFrom RColorBrewer    brewer.pal
#' @importFrom stringr         str_length
#' @importFrom grid            unit
#' @importFrom openxlsx        read.xlsx write.xlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom ggplot2 element_blank
#' @param plotData Data frame containing enrichment results with required columns (Category, Description, p.adjust, Count, geneRatio).
#' @param saveName Character string specifying the output file path for the plot.
#' @return Saves a dot plot to the specified file path.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[egg]{set_panel_size}}
#' @rdname enrich.Plot
#' @export
#' @importFrom egg set_panel_size
enrich.Plot <- function(plotData, saveName) {
  # Remove rows with any NA values
  plotData <- plotData %>% filter(!if_all(.cols = everything(), .fns = is.na))
  # Convert Category to factor with sorted levels
  plotData$Category <- factor(plotData$Category, levels = sort(unique(plotData$Category)))

  # Create dot plot
  Plot <- ggplot(plotData, aes(x = geneRatio, y = reorder(Description, -log10(p.adjust)),
                               color = -log10(p.adjust), size = Count, shape = Category)) +
    geom_point() +
    scale_color_gradient(high = "red", low = brewer.pal(6, "Reds")[2]) +
    xlab("GeneRatio") +
    ylab("") +
    theme_bw(base_size = 13) +
    theme(
      axis.text = element_text(size = 12, face = "plain", color = "black"),
      axis.text.x = element_text(angle = 0, size = 12),
      axis.title.y = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(t = 2, b = 2, r = 2, l = 2, unit = "cm")
    ) +
    guides(shape = guide_legend(order = 1, nrow = 2, byrow = TRUE),
           size = guide_legend(order = 2, nrow = 2, byrow = TRUE),
           colour = guide_colorbar(order = 3))

  # Save plot with dynamic sizing
  maxLength <- max(str_length(plotData$Description))
  ggsave(filename = saveName,
         egg::set_panel_size(Plot, width = unit(6, "cm"), height = unit(nrow(plotData) * 0.6, "cm")),
         width = 20, height = 18, dpi = 300, units = "cm", limitsize = FALSE, scale = 2)
}

#' @title Perform Functional Annotation
#' @description Performs GO, KEGG, and Hallmark enrichment analyses for gene modules and saves results.
#' @importFrom ggplot2  element_blank
#' @param subtype_file Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").
#' @param base_input_path Character string specifying the base directory containing module selection files (default: "./ModuleSelection/").
#' @param base_output_path Character string specifying the base directory for output files (default: "./FunctionalAnnotation/").
#' @param network_methodCharacter vector specifying the PPI network databases (options: "String", "physicalPPIN", "chengF").
#' @param module_method Character vector specifying the module detection methods (options: "Louvain", "WF").
#' @return Saves enrichment results and plots to the specified output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname functional_annotation
#' @export
functional_annotation <- function(subtype_file = "./subtype.xlsx",
                                  base_input_path = "./ModuleSelection/",
                                  base_output_path = "./FunctionalAnnotation/",
                                  network_method = c("String", "physicalPPIN", "chengF"),
                                  module_method = c("Louvain", "WF")) {

  # Validate network and module methods
  network_method <- match.arg(network_method, several.ok = TRUE)
  module_method <- match.arg(module_method, several.ok = TRUE)

  # Read and validate subtype file
  if (!file.exists(subtype_file)) {
    stop("Subtype file does not exist: ", subtype_file)
  }
  subtype_data <- read.xlsx(subtype_file)
  subtypes <- unique(subtype_data$Subtype)
  if (length(subtypes) == 0) {
    stop("No subtypes found in the 'Subtype' column of the subtype file.")
  }
  cat("Detected subtypes:", paste(subtypes, collapse = ", "), "\n")

  # Create base output directory
  if (!dir.exists(base_output_path)) {
    dir.create(base_output_path, recursive = TRUE)
    cat("Created output directory:", base_output_path, "\n")
  }

  # Define required columns for enrichment results
  required_columns <- c("ID", "Category", "Description", "GeneRatio", "BgRatio",
                        "RichFactor", "FoldEnrichment", "zScore", "pvalue",
                        "p.adjust", "qvalue", "geneID", "Count", "geneRatio")

  # Iterate over all combinations of network and module methods
  for (network in network_method) {
    for (method in module_method) {
      cat("\n=== Processing combination: Network method =", network, "| Module method =", method, "===\n")

      # Process each subtype
      for (subtype in subtypes) {
        cat("\nProcessing subtype:", subtype, "\n")

        # Construct input file path
        node_module_file <- file.path(base_input_path, network, method, subtype,
                                      paste0("node_Module_select_", network, "_", subtype, ".txt"))

        # Check if input file exists
        if (!file.exists(node_module_file)) {
          cat("Warning: Node-module file not found:", node_module_file, "\n")
          next
        }

        # Create subtype output directory
        subtype_path <- file.path(base_output_path, network, method, subtype)
        if (!dir.exists(subtype_path)) {
          dir.create(subtype_path, recursive = TRUE)
          cat("Created output directory:", subtype_path, "\n")
        }

        # Read node-module assignments
        Mnode <- tryCatch({
          read.table(node_module_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        }, error = function(e) {
          cat("Error reading file:", node_module_file, "\n")
          cat("Error message:", e$message, "\n")
          NULL
        })

        if (is.null(Mnode)) next

        colnames(Mnode) <- c("node", "module")

        # Check if there are any modules
        if (nrow(Mnode) == 0) {
          cat("No modules found in file:", node_module_file, "\n")
          next
        }

        # Get unique modules
        Module <- unique(Mnode$module)
        cat("Found", length(Module), "modules for subtype", subtype, "\n")

        # Initialize workbook for subtype's enrichment results
        wb <- createWorkbook()
        has_results <- FALSE

        # Process each module
        for (i in seq_along(Module)) {
          perM <- Module[i]
          cat("  Processing module", i, "of", length(Module), ":", perM, "\n")

          # Create module-specific output directory
          topathM <- file.path(subtype_path, perM)
          if (!dir.exists(topathM)) {
            dir.create(topathM, recursive = TRUE)
          }

          perMnode <- Mnode[Mnode$module == perM, "node"]

          # Initialize list to store top 15 results for this module
          module_results <- list()

          # GO Enrichment (BP)
          cat("    Running GO enrichment...\n")
          go <- enrich.GO(gene = perMnode, ont = "BP")
          if (!is.null(go)) {
            write.xlsx(go, file.path(topathM, paste0(perM, "_MsigDB.GO.", network, ".", subtype, ".xlsx")))
            go.top <- go %>%
              group_by(Category) %>%
              slice_min(order_by = p.adjust, n = 15) %>%
              ungroup() %>%
              select(all_of(required_columns))
            write.xlsx(go.top, file.path(topathM, paste0(perM, "_MsigDB.GO_top15.", network, ".", subtype, ".xlsx")))

            go.PlotName <- file.path(topathM, paste0(perM, "_MsigDB.GO_top15.", network, ".", subtype, ".pdf"))
            tryCatch({
              enrich.Plot(plotData = go.top, saveName = go.PlotName)
            }, error = function(e) {
              cat("Error plotting GO results:", e$message, "\n")
            })

            module_results[["GO"]] <- go.top
          }

          # KEGG Pathway Enrichment
          cat("    Running Pathway enrichment...\n")
          pathway <- enrich.Pathway(gene = perMnode, CP = "KEGG")
          if (!is.null(pathway)) {
            write.xlsx(pathway, file.path(topathM, paste0(perM, "_MsigDB.Pathway.", network, ".", subtype, ".xlsx")))
            pw.top <- pathway %>%
              group_by(Category) %>%
              slice_min(order_by = p.adjust, n = 15) %>%
              ungroup() %>%
              select(all_of(required_columns))
            write.xlsx(pw.top, file.path(topathM, paste0(perM, "_MsigDB.Pathway_top15.", network, ".", subtype, ".xlsx")))

            pw.PlotName <- file.path(topathM, paste0(perM, "_MsigDB.Pathway_top15.", network, ".", subtype, ".pdf"))
            tryCatch({
              enrich.Plot(plotData = pw.top, saveName = pw.PlotName)
            }, error = function(e) {
              cat("Error plotting Pathway results:", e$message, "\n")
            })

            module_results[["KEGG"]] <- pw.top
          }

          # Hallmark Enrichment
          cat("    Running Hallmark enrichment...\n")
          hallmark <- enrich.Hallmark(gene = perMnode)
          if (!is.null(hallmark)) {
            write.xlsx(hallmark, file.path(topathM, paste0(perM, "_MsigDB.Hallmark.", network, ".", subtype, ".xlsx")))
            hm.top <- hallmark %>%
              group_by(Category) %>%
              slice_min(order_by = p.adjust, n = 15) %>%
              ungroup() %>%
              select(all_of(required_columns))
            write.xlsx(hm.top, file.path(topathM, paste0(perM, "_MsigDB.Hallmark_top15.", network, ".", subtype, ".xlsx")))

            hm.PlotName <- file.path(topathM, paste0(perM, "_MsigDB.Hallmark_top15.", network, ".", subtype, ".pdf"))
            tryCatch({
              enrich.Plot(plotData = hm.top, saveName = hm.PlotName)
            }, error = function(e) {
              cat("Error plotting Hallmark results:", e$message, "\n")
            })

            module_results[["Hallmark"]] <- hm.top
          }

          # Combine top 15 results for this module
          combined_results <- bind_rows(module_results)
          if (nrow(combined_results) > 0) {
            addWorksheet(wb, sheetName = perM)
            writeData(wb, sheet = perM, x = combined_results)
            has_results <- TRUE
          }
        }

        # Save subtype's enrichment workbook if it contains sheets
        if (has_results) {
          output_file <- file.path(subtype_path, paste0(subtype, "_fun.xlsx"))
          saveWorkbook(wb, output_file, overwrite = TRUE)
          cat("Saved enrichment results to:", output_file, "\n")
        } else {
          cat("No enrichment results for subtype:", subtype, "\n")
        }
      }
    }
  }
  cat("\n=== Functional annotation completed ===\n")
}

#' @title Plot Enrichment Results for All Modules
#' @description  Generates dot plots for enrichment results stored in Excel files for each module.
#' @importFrom magrittr        %>%
#' @importFrom dplyr           mutate filter arrange desc group_by slice ungroup case_when
#' @importFrom stringr         str_trunc
#' @importFrom purrr           walk
#' @importFrom readxl          excel_sheets read_excel
#' @importFrom ggplot2         ggplot aes geom_point
#' @importFrom ggplot2         scale_color_gradientn scale_size_continuous labs
#' @importFrom ggplot2         theme_minimal theme element_text element_rect element_line
#' @importFrom ggplot2         margin rel
#' @importFrom ggplot2  element_blank
#' @importFrom RColorBrewer    brewer.pal
#' @importFrom grDevices       colorRampPalette
#' @param subtype_file Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").
#' @param base_input_path  Character string specifying the base directory containing functional annotation files (default: "./FunctionalAnnotation/").
#' @param network_method Character vector specifying the PPI network databases (options: "String", "physicalPPIN", "chengF").
#' @param module_method Character vector specifying the module detection methods (options: "Louvain", "WF").
#' @return Saves dot plots as PNG files to the specified input directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#' setwd("your_workspace")
#' Run functional annotation for all combinations
#' functional_annotation(subtype_file = "./subtype.xlsx",
#'                      network_method = c("String", "physicalPPIN", "chengF"),
#'                      module_method = c("Louvain", "WF"))
#'  }
#' }
#' @seealso
#'  \code{\link[scales]{pretty_breaks}}
#' @rdname plot_enrichment
#' @export
#' @importFrom scales pretty_breaks
plot_enrichment <- function(subtype_file = "./subtype.xlsx",
                            base_input_path = "./FunctionalAnnotation/",
                            network_method = c("String", "physicalPPIN", "chengF"),
                            module_method = c("Louvain", "WF")) {

  # Validate network and module methods
  network_method <- match.arg(network_method, several.ok = TRUE)
  module_method <- match.arg(module_method, several.ok = TRUE)

  # Read and validate subtype file
  if (!file.exists(subtype_file)) {
    stop("Subtype file does not exist: ", subtype_file)
  }
  subtype_data <- read.xlsx(subtype_file)
  subtypes <- unique(subtype_data$Subtype)
  if (length(subtypes) == 0) {
    stop("No subtypes found in the 'Subtype' column of the subtype file.")
  }
  cat("Detected subtypes:", paste(subtypes, collapse = ", "), "\n")

  # Internal Function: Preprocess Enrichment Data
  # Description: Preprocesses enrichment data for plotting by computing GeneRatio and log p.adjust.
  # Parameters:
  #   df: Data frame containing enrichment results.
  # Returns: Processed data frame.
  preprocess_data <- function(df) {
    df %>%
      mutate(
        GeneRatio = sapply(strsplit(GeneRatio, "/"),
                           function(x) as.numeric(x[1])/as.numeric(x[2])),
        log_padj = -log10(p.adjust),
        log_padj = ifelse(is.infinite(log_padj), NA, log_padj)
      ) %>%
      filter(!is.na(log_padj)) %>%
      arrange(desc(log_padj)) %>%
      group_by(Description) %>%
      slice(1) %>%
      ungroup()
  }

  # Internal Function: Create Single Enrichment Plot
  # Description: Generates a dot plot for a single module's enrichment results.
  # Parameters:
  #   df: Data frame containing preprocessed enrichment data.
  #   base_name: Character string for the base name of the output file.
  #   sheet_name: Character string for the sheet/module name.
  #   output_dir: Character string specifying the output directory.
  # Returns: Saves a PNG plot or NULL if the data frame is empty.
  create_single_plot <- function(df, base_name, sheet_name, output_dir) {
    # Skip empty data frames
    if (nrow(df) == 0) {
      message("Empty data for sheet: ", sheet_name)
      return(NULL)
    }

    # Truncate long descriptions
    max_desc_length <- 100
    df <- df %>%
      mutate(Description = str_trunc(Description, max_desc_length))

    # Adjust point size based on number of rows
    point_size <- case_when(
      nrow(df) > 30 ~ 2,
      nrow(df) > 15 ~ 3,
      TRUE ~ 4
    )

    # Dynamic color range
    color_range <- range(df$log_padj, na.rm = TRUE)

    # Create dot plot
    p <- ggplot(df, aes(x = GeneRatio, y = reorder(Description, log_padj),
                        color = log_padj, size = Count, shape = Category)) +
      geom_point(alpha = 0.8) +
      scale_color_gradientn(
        colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
        limits = color_range,
        name = "-log10(p.adjust)"
      ) +
      scale_size_continuous(
        range = c(point_size, point_size * 2),
        breaks = scales::pretty_breaks(n = 4)
      ) +
      labs(
        title = paste("Enrichment Analysis:", sheet_name),
        x = "GeneRatio",
        y = ""
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(color = "gray40"),
        legend.box.background = element_rect(color = "gray80"),
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 2, 1, 1, "cm")
      )

    # Adjust plot dimensions
    plot_height <- max(15, nrow(df) * 0.4)
    plot_width <- 15 + max(nchar(df$Description)) * 0.15

    # Generate unique filename
    clean_name <- gsub("[^[:alnum:]]", "_", sheet_name)
    file_name <- file.path(output_dir, paste0(basename(base_name), "_", clean_name, ".png"))

    # Save plot
    ggsave(
      filename = file_name,
      plot = p,
      width = plot_width,
      height = plot_height,
      units = "cm",
      dpi = 300,
      limitsize = FALSE
    )
  }

  # Iterate over all combinations of network and module methods
  for (network in network_method) {
    for (method in module_method) {
      cat("Processing plots: Network method:", network, "Module method:", method, "\n")

      # Process each subtype
      for (subtype in subtypes) {
        cat("Processing plots for subtype:", subtype, "\n")

        # Construct input file path
        enrich_file <- file.path(base_input_path, network, method, subtype,
                                 paste0(subtype, "_fun.xlsx"))

        # Check if enrichment file exists
        if (!file.exists(enrich_file)) {
          cat("Enrichment file missing for network:", network, "method:", method,
              "subtype:", subtype, "at", enrich_file, "\n")
          next
        }

        # Get sheet names (modules)
        sheets <- excel_sheets(enrich_file)

        # Process each sheet (module)
        walk(sheets, ~{
          data <- read_excel(enrich_file, sheet = .x) %>%
            preprocess_data()

          create_single_plot(data, paste0(subtype, "_fun"), .x,
                             file.path(base_input_path, network, method, subtype))
        })

        cat("Plotting completed for subtype:", subtype, "\n")
      }
    }
  }
}
