#' @title Module Division Function for PPI Networks
#' @description This function performs module detection on PPI networks for specified subtypes using selected network databases and clustering methods (Louvain or WF).
#' @importFrom openxlsx read.xlsx
#' @importFrom igraph graph_from_data_frame V vcount ecount
#' @importFrom igraph cluster_louvain cluster_edge_betweenness
#' @importFrom igraph cluster_label_prop membership
#' @importFrom stats phyper
#' @importFrom utils read.table write.table
#' @param subtype_file Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").
#' @param base_input_path Character string specifying the base directory containing PPI network files (default: "./Netconstruct_result/").
#' @param output_base_path Character string specifying the base directory for output files (default: "./ModuleDivision/").
#' @param network_method Character vector specifying the PPI network databases (options: "String", "physicalPPIN", "chengF").
#' @param module_method Character vector specifying the module detection methods (options: "Louvain", "WF").
#' @return Saves node and edge module assignments to text files in the specified output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#'Set working directory and run module division for all combinations
#'setwd("your_workspace")
#'module_division(network_method = c("String", "physicalPPIN", "chengF"),
#'                 module_method = c("Louvain", "WF"))
#' }
#' @rdname module_division
#' @export
module_division <- function(subtype_file = "./subtype.xlsx",
                            base_input_path = "./Netconstruct_result/",
                            output_base_path = "./ModuleDivision/",
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
  module_division_path <- output_base_path
  dir.create(module_division_path, showWarnings = FALSE, recursive = TRUE)

  # Iterate over all combinations of network and module methods
  for (network in network_method) {
    for (method in module_method) {
      cat("Processing: Network method:", network, " Module method:", method, "\n")

      # Process each subtype
      for (subtype in subtypes) {
        # Construct path to PPI file
        ppi_file_path <- file.path(base_input_path, subtype, network, paste0("ppi_", subtype, ".txt"))

        # Check if PPI file exists
        if (!file.exists(ppi_file_path)) {
          cat("PPI file missing for network:", network, "subtype:", subtype, "at", ppi_file_path, "\n")
          next
        }

        # Create subtype-specific output directory
        subtype_output_path <- file.path(module_division_path, network, method, subtype)
        dir.create(subtype_output_path, showWarnings = FALSE, recursive = TRUE)

        # Analyze subtype for module detection
        subtype_module(ppi_file_path, subtype_output_path, network, method)
      }
    }
  }
}

#' @title Analyze Subtype for Module Detection
#' @description Performs module detection on a PPI network for a given subtype using the specified method.
#' @importFrom igraph graph_from_data_frame V vcount ecount
#' @importFrom igraph cluster_louvain cluster_edge_betweenness cluster_label_prop membership
#' @importFrom stats phyper
#' @importFrom utils read.table write.table
#' @param ppi_file_path Character string specifying the path to the PPI network file.
#' @param output_base_path  Character string specifying the output directory for results.
#' @param network_method Character string specifying the network database used.
#' @param module_method Character string specifying the module detection method ("Louvain" or "WF").
#' @return Saves node and edge module assignments to text files.
#' @details DETAILS
#' @examples
#' \dontrun{
#' # setwd("your_workspace")
#' module_division(network_method = c("String", "physicalPPIN", "chengF"),
#' module_method = c("Louvain", "WF"))
#'  }
#' @rdname subtype_module
#' @export
subtype_module <- function(ppi_file_path, output_base_path, network_method, module_method) {
  # Read PPI network data
  ppi_data <- read.table(ppi_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # Create igraph object from PPI data
  edges <- ppi_data[, 1:2]
  colnames(edges) <- c("node1", "node2")
  nodes <- unique(c(edges$node1, edges$node2))  # Extract unique nodes
  net <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

  # Log network details
  cat("igraph object created successfully.\n")
  cat("Number of nodes:", vcount(net), "\n")
  cat("Number of edges:", ecount(net), "\n")

  if (module_method == "Louvain") {
    # Perform Louvain clustering for module detection
    set.seed(123)  # Ensure reproducibility
    cluster_louvain_result <- cluster_louvain(net, weights = NULL, resolution = 1)

    # Log number of detected modules
    num_modules <- length(unique(membership(cluster_louvain_result)))
    cat("Louvain clustering detected", num_modules, "modules.\n")

    # Create node-module assignment data frame
    node_module <- data.frame(
      node = cluster_louvain_result$names,
      module = cluster_louvain_result$membership,
      stringsAsFactors = FALSE
    )

    # Save node-module assignments
    write.table(node_module,
                file = file.path(output_base_path, paste0("node_Module_", network_method, "_", module_method, ".txt")),
                row.names = FALSE, sep = "\t", quote = FALSE)

    # Annotate edges with module IDs
    edges$module <- apply(edges, 1, function(x) {
      node1_module <- node_module$module[node_module$node == x[1]]
      node2_module <- node_module$module[node_module$node == x[2]]
      if (length(node1_module) == 0 || length(node2_module) == 0) {
        return(0)  # Node not assigned to a module
      } else if (node1_module == node2_module) {
        return(node1_module)  # Same module
      } else {
        return(0)  # Different modules
      }
    })

    # Save edge-module assignments
    write.table(edges,
                file = file.path(output_base_path, paste0("edges_", network_method, "_", module_method, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = c("node1", "node2", "module"))

    # Filter and save intramodule edges
    edges_module <- subset(edges, module != 0)
    cat("Number of intramodule edges:", nrow(edges_module), "\n")
    write.table(edges_module,
                file = file.path(output_base_path, paste0("edge_Module_", network_method, "_", module_method, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = c("node1", "node2", "module"))

  } else if (module_method == "WF") {
    # Perform WF (Walktrap-Fan) clustering: Combines Edge Betweenness and Label Propagation
    set.seed(123)  # Ensure reproducibility

    # Edge Betweenness (Girvan-Newman) clustering
    fan_GN <- cluster_edge_betweenness(net, weights = NULL)
    cat("Edge Betweenness clustering detected", length(unique(membership(fan_GN))), "modules.\n")

    # Extract node-module assignments for Edge Betweenness
    fan_GN_label <- data.frame(
      node = V(net)$name,
      module = membership(fan_GN),
      stringsAsFactors = FALSE
    )

    # Label Propagation clustering
    fan_LP <- cluster_label_prop(net, weights = NULL)
    cat("Label Propagation clustering detected", length(unique(membership(fan_LP))), "modules.\n")

    # Extract node-module assignments for Label Propagation
    fan_LP_label <- data.frame(
      node = V(net)$name,
      module = membership(fan_LP),
      stringsAsFactors = FALSE
    )

    # Split nodes by module for both methods
    fan_GN_label_list <- split(fan_GN_label$node, fan_GN_label$module)
    fan_LP_label_list <- split(fan_LP_label$node, fan_LP_label$module)

    # Filter out single-node modules
    fan_GN_label_list2 <- fan_GN_label_list[lengths(fan_GN_label_list) > 1]
    fan_LP_label_list2 <- fan_LP_label_list[lengths(fan_LP_label_list) > 1]

    # Log effective module counts
    cat("Effective Edge Betweenness modules (>1 node):", length(fan_GN_label_list2), "\n")
    cat("Effective Label Propagation modules (>1 node):", length(fan_LP_label_list2), "\n")

    # Initialize matrices for intersection, union, and hypergeometric p-values
    fan_GN_LP2 <- matrix(0, nrow = length(fan_LP_label_list2), ncol = length(fan_GN_label_list2))
    fan_GN_LP_union2 <- matrix(0, nrow = length(fan_LP_label_list2), ncol = length(fan_GN_label_list2))
    fan_GN_LP_phyper2 <- matrix(0, nrow = length(fan_LP_label_list2), ncol = length(fan_GN_label_list2))

    # Set matrix row and column names
    colnames(fan_GN_LP2) <- paste("GN", seq_len(length(fan_GN_label_list2)), sep = "_")
    rownames(fan_GN_LP2) <- paste("LP", seq_len(length(fan_LP_label_list2)), sep = "_")

    # Compute intersection and union sizes for module pairs
    for (i in seq_len(length(fan_GN_label_list2))) {
      for (j in seq_len(length(fan_LP_label_list2))) {
        intersection_size <- length(intersect(fan_GN_label_list2[[i]], fan_LP_label_list2[[j]]))
        union_size <- length(union(fan_GN_label_list2[[i]], fan_LP_label_list2[[j]]))
        fan_GN_LP2[j, i] <- intersection_size
        fan_GN_LP_union2[j, i] <- union_size
      }
    }

    # Perform hypergeometric test for module overlap significance
    for (i in seq_len(ncol(fan_GN_LP_phyper2))) {
      for (j in seq_len(nrow(fan_GN_LP_phyper2))) {
        fan_GN_LP_phyper2[j, i] <- 1 - phyper(
          fan_GN_LP2[j, i],
          length(fan_GN_label_list2[[i]]),
          fan_GN_LP_union2[j, i],
          length(fan_LP_label_list2[[j]])
        )
      }
    }

    # Identify significant module pairs (p < 0.05)
    sig_fan <- as.data.frame(which(fan_GN_LP_phyper2 < 0.05, arr.ind = TRUE))
    if (nrow(sig_fan) == 0) {
      cat("No significant module pairs found (p < 0.05).\n")
      return()  # Exit if no significant pairs
    } else {
      cat("Number of significant module pairs:", nrow(sig_fan), "\n")
    }

    # Extract nodes for significant module pairs
    sig_fan_list <- list()
    for (i in seq_len(nrow(sig_fan))) {
      gn_index <- sig_fan[i, "col"]
      lp_index <- sig_fan[i, "row"]
      intersect_nodes <- intersect(fan_GN_label_list2[[gn_index]], fan_LP_label_list2[[lp_index]])
      sig_fan_list[[i]] <- intersect_nodes
    }

    # Assign new module IDs to nodes
    total_nodes <- sum(sapply(sig_fan_list, length))
    modules <- as.data.frame(matrix(NA, nrow = total_nodes, ncol = 2))
    w <- f <- 0
    for (i in seq_len(length(sig_fan_list))) {
      w <- f + 1
      f <- f + length(sig_fan_list[[i]])
      modules[w:f, 1] <- unlist(sig_fan_list[[i]])
      modules[w:f, 2] <- i
    }
    colnames(modules) <- c("node", "module")

    # Save node-module assignments
    write.table(modules,
                file = file.path(output_base_path, paste0("node_Module_", network_method, "_", module_method, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = c("node", "module"))

    # Annotate edges with module IDs
    edges$module <- 0  # Initialize
    for (i in seq_len(length(sig_fan_list))) {
      pattern <- as.vector(modules[which(modules$module == i), 1])
      for (j in seq_len(nrow(edges))) {
        if (edges[j, 1] %in% pattern && edges[j, 2] %in% pattern) {
          edges$module[j] <- i
        }
      }
    }

    # Save edge-module assignments
    write.table(edges,
                file = file.path(output_base_path, paste0("edges_", network_method, "_", module_method, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = c("node1", "node2", "module"))

    # Filter and save intramodule edges
    edges_module <- subset(edges, module != 0)
    cat("Number of intramodule edges:", nrow(edges_module), "\n")
    write.table(edges_module,
                file = file.path(output_base_path, paste0("edge_Module_", network_method, "_", module_method, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = c("node1", "node2", "module"))

    # Log final statistics
    cat("Total edges:", nrow(edges), "\n")
    cat("Total modules:", length(sig_fan_list), "\n")
  }

  cat("Module detection completed!\n")
}
