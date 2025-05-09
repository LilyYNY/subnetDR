#' @title Module Selection Function for PPI Network Modules
#' @description This function selects modules from PPI networks based on node count thresholds and generates output files and visualizations for selected modules.
#' @importFrom openxlsx read.xlsx
#' @importFrom utils read.table write.table write.csv
#' @importFrom stats aggregate
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip geom_text labs theme_classic theme element_text
#' @importFrom grDevices tiff dev.off
#' @param subtype_file Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").
#' @param base_input_path Character string specifying the base directory containing module division files (default: "./ModuleDivision/").
#' @param base_output_path Character string specifying the base directory for output files (default: "./ModuleSelection/").
#' @param network_method PCharacter vector specifying the PPI network databases (options: "String", "physicalPPIN", "chengF"; default: all).
#' @param module_method Character vector specifying the module detection methods (options: "Louvain", "WF"; default: both).
#' @param numberCutoff Numeric threshold for minimum number of nodes in a module (default: 9).
#' @return Saves selected module data (nodes, edges, counts) and bar plots to the specified output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#' setwd("your_workspace")
#' Run module selection for all combinations of network and module methods
#' module_selection(subtype_file = "./subtype.xlsx",
#'                 network_method = c("String", "physicalPPIN", "chengF"),
#'                 module_method = c("Louvain", "WF"),
#'                 numberCutoff = 9)
#'  }
#' }
#' @rdname module_selection
#' @export
module_selection <- function(subtype_file = "./subtype.xlsx",
                             base_input_path = "./ModuleDivision/",
                             base_output_path = "./ModuleSelection/",
                             network_method = c("String", "physicalPPIN", "chengF"),
                             module_method = c("Louvain", "WF"),
                             numberCutoff = 9) {

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
  dir.create(base_output_path, showWarnings = FALSE, recursive = TRUE)

  # Iterate over all combinations of network and module methods
  for (network in network_method) {
    for (method in module_method) {
      cat("Processing: Network method:", network, "Module method:", method, "\n")

      # Process each subtype
      for (subtype in subtypes) {
        cat("Processing module selection for subtype:", subtype, "\n")

        # Construct input file paths
        node_module_file <- file.path(base_input_path, network, method, subtype,
                                      paste0("node_Module_", network, "_", method, ".txt"))
        edge_file <- file.path(base_input_path, network, method, subtype,
                               paste0("edges_", network, "_", method, ".txt"))

        # Check if input files exist
        if (!file.exists(node_module_file) || !file.exists(edge_file)) {
          cat("Input file missing for network:", network, "method:", method, "subtype:", subtype, "\n")
          cat("Node file:", node_module_file, "\n")
          cat("Edge file:", edge_file, "\n")
          next
        }

        # Create subtype-specific output directory
        output_path <- file.path(base_output_path, network, method, subtype)
        dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

        # Step 1: Read node-module assignments
        node_module <- read.table(node_module_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        colnames(node_module) <- c("node", "module")
        node_module$module <- paste("M", node_module$module, sep = "")  # Add "M" prefix to module IDs

        # Step 2: Count nodes per module
        perModule_num <- aggregate(node_module$node, by = list(module = node_module$module), length)
        colnames(perModule_num) <- c("module", "count")
        perModule_num <- perModule_num[order(-perModule_num$count), ]  # Sort by node count (descending)

        # Save module node counts to CSV
        write.csv(perModule_num,
                  file = file.path(output_path, paste0("Module_Number_", network, "_", subtype, ".csv")),
                  row.names = FALSE, quote = FALSE)

        # Step 3: Select modules with node count greater than numberCutoff
        module_select <- subset(perModule_num, count > numberCutoff)
        module_select <- module_select[order(module_select$count, decreasing = TRUE), ]  # Sort by count

        # Save selected modules to CSV
        write.csv(module_select,
                  file = file.path(output_path, paste0("Module_select_", network, "_", subtype, ".csv")),
                  row.names = FALSE, quote = FALSE)

        # Filter nodes in selected modules
        node_module_sel <- subset(node_module, node_module$module %in% module_select$module)
        # Save selected nodes to text file
        write.table(node_module_sel,
                    file = file.path(output_path, paste0("node_Module_select_", network, "_", subtype, ".txt")),
                    row.names = FALSE, sep = "\t", quote = FALSE)

        # Step 4: Read and filter edge data
        edgeALL <- read.table(edge_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        edgeALL$module <- paste("M", edgeALL$module, sep = "")  # Add "M" prefix to module IDs

        # Filter edges where both nodes are in selected modules
        edge_module_sel <- subset(edgeALL, edgeALL$node1 %in% node_module_sel$node & edgeALL$node2 %in% node_module_sel$node)
        # Save selected edges to text file
        write.table(edge_module_sel,
                    file = file.path(output_path, paste0("edges_select_", network, "_", subtype, ".txt")),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Step 5: Create bar plot of module node counts
        perModule_num$color <- ifelse(perModule_num$count > numberCutoff, "lightslateblue", "gray")  # Color based on cutoff
        perModule_num$module <- factor(perModule_num$module, levels = perModule_num$module[order(perModule_num$count)])  # Order modules
        yColor <- rev(perModule_num$color)  # Reverse color for y-axis text

        # Generate bar plot
        p <- ggplot(perModule_num, aes(x = module, y = count)) +
          geom_bar(stat = 'identity', fill = perModule_num$color, colour = perModule_num$color, width = 0.8) +
          coord_flip() +  # Flip coordinates for horizontal bars
          geom_text(aes(label = count, hjust = -0.2), size = 4) +  # Add count labels
          labs(x = "", y = "", fill = "") +  # Remove axis labels
          theme_classic() +  # Use classic theme
          theme(axis.text.y = element_text(size = 12, colour = yColor),  # Customize y-axis text
                axis.text = element_text(size = 11, face = "plain"))  # Customize general text

        # Save plot as TIFF
        file_name <- file.path(output_path, paste0("Module_Number_", network, "_", subtype, ".tiff"))
        tiff(file_name, width = 2500, height = 2500, res = 300)
        print(p)
        dev.off()

        cat("Module selection for subtype", subtype, "completed.\n\n")
      }
    }
  }
}
