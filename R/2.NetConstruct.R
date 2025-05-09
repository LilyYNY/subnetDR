#' @title Protein-Protein Interaction (PPI) Network Construction Function
#' @description This function constructs PPI networks for differentially expressed proteins across subtypes using specified PPI databases (String, PhysicalPPIN, ChengF).
#' @importFrom openxlsx read.xlsx
#' @param diff_file  Character string specifying the path to the differential expression results Excel file.
#' @param ppiScore Numeric threshold for PPI interaction score (default: 400, used for String and ChengF methods).
#' @param subtype_file Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").
#' @param output_dir Character string specifying the output directory for PPI network files (default: "./Netconstruct_result/").
#' @param ppi_method Character vector specifying the PPI databases to use (options: "String", "PhysicalPPIN", "ChengF").
#' @return Saves PPI network files for each subtype and method to the specified output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' setwd("your_workspace")
#' run_network_construction(
#'   diff_file = "./DEPS_result/diff_expression_results_all.xlsx",
#'   ppi_method = c("String", "PhysicalPPIN", "ChengF"),
#'   output_dir = "./Netconstruct_result/")
#'}
#' @rdname run_network_construction
#' @export
run_network_construction <- function(diff_file,  # Input differential expression results file
                                     ppiScore = 400,
                                     subtype_file = "./subtype.xlsx",  # Phenotype information file
                                     output_dir = "./Netconstruct_result/",
                                     ppi_method = c("String", "PhysicalPPIN", "ChengF")) {

  # Read phenotype data from Excel file
  subtype_data <- read.xlsx(subtype_file)

  # Extract unique subtypes from phenotype data
  subtypes <- unique(subtype_data$Subtype)  # List of all unique subtype names

  # Internal Function: Process PPI network construction for a single subtype
  # Parameters:
  #   subtype: Character string specifying the subtype to process.
  #   ppiScore: Numeric threshold for PPI interaction score.
  #   diff_file: Path to the differential expression results Excel file.
  #   output_dir: Output directory for saving PPI network files.
  #   ppi_method: Vector of PPI databases to use.
  # Returns: Saves PPI network files for the specified subtype and methods.
  process_subtype <- function(subtype, ppiScore, diff_file, output_dir, ppi_method) {

    # Read differential expression results for the specified subtype
    dep_data <- read.xlsx(diff_file, sheet = paste0(subtype, "_DiffResults"))

    # Extract upregulated proteins (where label is "up")
    dep_data_up <- subset(dep_data, label == "up")  # Filter for upregulated proteins
    proteins <- dep_data_up$Gene  # Extract gene/protein names

    # Create subtype-specific output directory
    subtype_path <- paste0(output_dir, subtype, "/")
    dir.create(subtype_path, showWarnings = FALSE, recursive = TRUE)  # Create directory recursively

    # Process PPI network construction for String method
    if ("String" %in% ppi_method) {
      # Create String method-specific subdirectory
      subtype_path_String <- paste0(subtype_path, "String/")
      dir.create(subtype_path_String, showWarnings = FALSE, recursive = TRUE)

      # Read String PPI data files
      ppiID_file <- "./PPI/String/9606.protein.info.v11.5_new.txt"  # Protein ID to symbol mapping
      ppiALL_file <- "./PPI/String/9606.protein.links.v11.5.txt"    # PPI interaction data
      ppiID <- read.table(ppiID_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      ppiALL <- read.table(ppiALL_file, header = TRUE, stringsAsFactors = FALSE)

      # Prepare ID mapping table and filter for differentially expressed proteins
      ID.ann <- ppiID
      colnames(ID.ann) <- c("ID", "Symbol")
      ID.ann <- ID.ann[ID.ann$Symbol %in% proteins, ]
      rownames(ID.ann) <- ID.ann$ID

      # Filter PPI interactions based on score threshold
      ppi <- ppiALL
      colnames(ppi) <- c("node1", "node2", "score")
      ppi <- subset(ppi, ppi$score > ppiScore)

      # Filter PPI interactions to include only those with both proteins in the input list
      ppi$match <- apply(ppi, 1, function(x) ifelse(x[1] %in% ID.ann$ID & x[2] %in% ID.ann$ID, "yes", "no"))
      ppiMatch <- ppi[ppi$match == "yes", c("node1", "node2")]

      # Convert node IDs to gene symbols
      node1_symbol <- ID.ann[as.character(ppiMatch$node1), "Symbol"]
      node2_symbol <- ID.ann[as.character(ppiMatch$node2), "Symbol"]
      ppiMatch$node1 <- node1_symbol
      ppiMatch$node2 <- node2_symbol

      # Remove duplicate PPIs (e.g., A-B and B-A)
      ppiMatch <- ppiMatch[row.names(unique(t(apply(ppiMatch[, 1:2], 1, sort)))), ]

      # Save PPI network to file
      output_file <- paste0(subtype_path_String, "ppi_", subtype, ".txt")
      write.table(ppiMatch, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("PPI network for", subtype, "saved to:", output_file, "\n")
    }

    # Process PPI network construction for PhysicalPPIN method
    if ("PhysicalPPIN" %in% ppi_method) {
      # Create PhysicalPPIN method-specific subdirectory
      subtype_path_Physical <- paste0(subtype_path, "PhysicalPPIN/")
      dir.create(subtype_path_Physical, showWarnings = FALSE, recursive = TRUE)

      # Read PhysicalPPIN data file
      ppiALL_file <- "./PPI/physicalPPIN/physicalppi.txt"
      ppiALL <- read.table(ppiALL_file, fill = TRUE, stringsAsFactors = FALSE)

      # Handle missing values in PPI data
      ppiALL[is.na(ppiALL)] <- NA

      # Assign column names to PPI data
      colnames(ppiALL) <- c("node1", "node2")

      # Filter PPI interactions to include only those with both proteins in the input list
      ppiALL$match <- apply(ppiALL, 1, function(x) ifelse(x[1] %in% proteins & x[2] %in% proteins, "yes", "no"))
      ppiMatch <- ppiALL[ppiALL$match == "yes", c("node1", "node2")]

      # Remove duplicate PPIs (e.g., A-B and B-A)
      ppiMatch <- ppiMatch[row.names(unique(t(apply(ppiMatch[, 1:2], 1, sort)))), ]

      # Save PPI network to file
      output_file <- paste0(subtype_path_Physical, "ppi_", subtype, ".txt")
      write.table(ppiMatch, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("PPI network for", subtype, "saved to:", output_file, "\n")
    }

    # Process PPI network construction for ChengF method
    if ("ChengF" %in% ppi_method) {
      # Create ChengF method-specific subdirectory
      subtype_path_ChengF <- paste0(subtype_path, "ChengF/")
      dir.create(subtype_path_ChengF, showWarnings = FALSE, recursive = TRUE)

      # Read ChengF PPI data files
      ppiID_file <- "./PPI/ChengF/Homo_sapiens_gene_info.txt"  # Gene ID to symbol mapping
      ppiALL_file <- "./PPI/ChengF/Human Interactome.txt"      # PPI interaction data
      gene_info <- read.table(ppiID_file, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
      ppiALL <- read.table(ppiALL_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

      # Extract Gene ID to Symbol mapping and remove missing values
      ppiID <- gene_info[, c("GeneID", "Symbol")]
      ppiID <- na.omit(ppiID)

      # Prepare ID mapping table and filter for differentially expressed proteins
      ID.ann <- ppiID
      colnames(ID.ann) <- c("ID", "Symbol")
      ID.ann <- ID.ann[ID.ann$Symbol %in% proteins, ]
      rownames(ID.ann) <- ID.ann$ID

      # Select relevant PPI columns
      ppi <- ppiALL[, c("Protein.A", "Protein.B")]
      colnames(ppi) <- c("node1", "node2")

      # Filter PPI interactions to include only those with both proteins in the input list
      ppi$match <- apply(ppi, 1, function(x) ifelse(x[1] %in% ID.ann$ID & x[2] %in% ID.ann$ID, "yes", "no"))
      ppiMatch <- ppi[ppi$match == "yes", c("node1", "node2")]

      # Convert node IDs to gene symbols
      node1_symbol <- ID.ann[as.character(ppiMatch$node1), "Symbol"]
      node2_symbol <- ID.ann[as.character(ppiMatch$node2), "Symbol"]
      ppiMatch$node1 <- node1_symbol
      ppiMatch$node2 <- node2_symbol

      # Remove duplicate PPIs (e.g., A-B and B-A)
      ppiMatch <- ppiMatch[row.names(unique(t(apply(ppiMatch[, 1:2], 1, sort)))), ]

      # Save PPI network to file
      output_file <- paste0(subtype_path_ChengF, "ppi_", subtype, ".txt")
      write.table(ppiMatch, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("PPI network for", subtype, "saved to:", output_file, "\n")
    }
  }

  # Process each subtype to construct PPI networks
  for (subtype in subtypes) {
    process_subtype(subtype, ppiScore, diff_file, output_dir, ppi_method)
  }
}
