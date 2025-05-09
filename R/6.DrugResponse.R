#' @title Drug Response Analysis for Subtype-Specific Gene Modules
#' @description This function performs drug response analysis using oncoPredict and ConsensusClusterPlus  to predict drug responses and cluster samples based on module-specific gene expression.
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of rename mutate_all left_join distinct
#' @importFrom readxl read_excel
#' @importFrom openxlsx read.xlsx write.xlsx getSheetNames
#' @importFrom utils read.delim write.table read.csv write.csv packageVersion
#' @importFrom oncoPredict calcPhenotype
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom stats wilcox.test
#' @param base_path Character string specifying the base directory for all input/output paths.
#' @param subtype_file Character string specifying the path to the subtype phenotype data file.
#' @param expr_base_path Character string specifying the directory for expression data output.
#' @param expr_data_file Character string specifying the path to the expression data file.
#' @param diff_expr_file Character string specifying the path to the differential expression results file.
#' @param module_base_path Character string specifying the directory containing module selection files.
#' @param drug_response_path Character string specifying the directory for drug response output.
#' @param gdsc_expr_file Character string specifying the path to the GDSC expression data file.
#' @param gdsc_res_file Character string specifying the path to the GDSC response data file.
#' @param network_methods Character vector specifying the PPI network databases (e.g., "String", "physicalPPIN", "chengF").
#' @param module_methods Character vector specifying the module detection methods (e.g., "Louvain", "WF").
#' @param gene_col_name Character string specifying the column name for genes in expression data.
#' @return Saves drug response predictions, clustering results, and drug response networks (DRNs) to the specified output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  drug_response_analysis(
#'  base_path = "your_workspace",
#'  subtype_file = "subtype.xlsx",
#'  expr_base_path = "DrugResponse/subtype_expression",
#'  expr_data_file = "expression.xlsx",
#'  diff_expr_file = "DEPS_result/diff_expression_results_all.xlsx",
#'  module_base_path = "ModuleSelection",
#'  drug_response_path = "DrugResponse",
#'  gdsc_expr_file = "GDSC2_Expr.RData",
#'  gdsc_res_file = "GDSC2_Res.RData",
#'  network_methods = c("String", "physicalPPIN", "chengF"),
#'  module_methods = c("Louvain", "WF"),
#'  gene_col_name = "Gene")
#'  }
#' }
#' @rdname drug_response_analysis
#' @export
drug_response_analysis <- function(
    base_path,
    subtype_file,
    expr_base_path,
    expr_data_file,
    diff_expr_file,
    module_base_path,
    drug_response_path,
    gdsc_expr_file,
    gdsc_res_file,
    network_methods,
    module_methods,
    gene_col_name
) {

  # Step 1: Check package versions
  oncoPredict_version <- packageVersion("oncoPredict")
  cat("Using oncoPredict version:", as.character(oncoPredict_version), "\n")
  cc_version <- packageVersion("ConsensusClusterPlus")
  cat("Using ConsensusClusterPlus version:", as.character(cc_version), "\n")

  # Step 2: Resolve and normalize file paths
  if (!is.null(base_path)) {
    base_path <- normalizePath(base_path, mustWork = FALSE)
    subtype_file <- file.path(base_path, subtype_file)
    expr_base_path <- file.path(base_path, expr_base_path)
    expr_data_file <- file.path(base_path, expr_data_file)
    diff_expr_file <- file.path(base_path, diff_expr_file)
    module_base_path <- file.path(base_path, module_base_path)
    drug_response_path <- file.path(base_path, drug_response_path)
    gdsc_expr_file <- file.path(base_path, gdsc_expr_file)
    gdsc_res_file <- file.path(base_path, gdsc_res_file)
  }

  # Normalize all paths
  subtype_file <- normalizePath(subtype_file, mustWork = FALSE)
  expr_base_path <- normalizePath(expr_base_path, mustWork = FALSE)
  expr_data_file <- normalizePath(expr_data_file, mustWork = FALSE)
  diff_expr_file <- normalizePath(diff_expr_file, mustWork = FALSE)
  module_base_path <- normalizePath(module_base_path, mustWork = FALSE)
  drug_response_path <- normalizePath(drug_response_path, mustWork = FALSE)
  gdsc_expr_file <- normalizePath(gdsc_expr_file, mustWork = FALSE)
  gdsc_res_file <- normalizePath(gdsc_res_file, mustWork = FALSE)

  # Create output directories
  dir.create(expr_base_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(drug_response_path, showWarnings = FALSE, recursive = TRUE)

  # Step 3: Read subtype file and generate differential expression sheet names
  cat("Reading subtype file and generating differential expression sheet names...\n")
  if (!file.exists(subtype_file)) {
    stop("Subtype file does not exist: ", subtype_file)
  }
  subtype_data <- read.xlsx(subtype_file)
  if (!all(c("Sample", "Subtype") %in% colnames(subtype_data))) {
    stop("subtype.xlsx must contain 'Sample' and 'Subtype' columns")
  }
  subtypes <- unique(subtype_data$Subtype)
  if (length(subtypes) == 0) {
    stop("No subtypes found in the 'Subtype' column of the subtype file.")
  }
  cat("Detected subtypes:", paste(subtypes, collapse = ", "), "\n")

  # Generate sheet names for differential expression results
  diff_sheets <- setNames(paste0(subtypes, "_DiffResults"), subtypes)
  cat("Differential expression sheets:", paste(names(diff_sheets), "=", diff_sheets, collapse = ", "), "\n")

  # Step 4: Generate subtype-specific expression files
  cat("Generating subtype-specific expression files...\n")
  if (!file.exists(expr_data_file)) {
    stop("Expression data file does not exist: ", expr_data_file)
  }
  expr_data <- read_excel(expr_data_file)

  # Validate gene column name
  if (!gene_col_name %in% colnames(expr_data)) {
    alternatives <- c("Gene", "gene", "Gene.name", "gene_name", "Gene_name", "Gene name")
    found <- FALSE
    for (alt in alternatives) {
      if (alt %in% colnames(expr_data)) {
        gene_col_name <- alt
        found <- TRUE
        break
      }
    }
    if (!found) {
      stop("Gene column not found in expression data. Tried: ",
           paste(c(gene_col_name, alternatives), collapse = ", "))
    }
  }
  cat("Using gene column name:", gene_col_name, "\n")

  # Create expression files for each subtype
  for (subtype in subtypes) {
    subtype_samples <- subtype_data$Sample[subtype_data$Subtype == subtype]
    if (length(subtype_samples) == 0) {
      cat("No samples found for subtype:", subtype, "Skipping...\n")
      next
    }

    valid_samples <- subtype_samples[subtype_samples %in% colnames(expr_data)]
    if (length(valid_samples) == 0) {
      cat("No valid samples in expression data for subtype:", subtype, "Skipping...\n")
      next
    }

    subtype_expr <- expr_data %>%
      select(all_of(gene_col_name), all_of(valid_samples))

    expr_file <- file.path(expr_base_path, paste0(subtype, "_expression.xlsx"))
    tryCatch({
      write.xlsx(subtype_expr, expr_file, rowNames = FALSE)
      cat(subtype, "expression data saved to:", expr_file, "\n")
    }, error = function(e) {
      cat("Failed to save expression data to:", expr_file, "Error:", e$message, "\n")
      next
    })
  }

  # Step 5: Generate filtered expression files for upregulated genes
  cat("Generating filtered expression files for upregulated genes...\n")
  if (!file.exists(diff_expr_file)) {
    stop("Differential expression file does not exist: ", diff_expr_file)
  }

  for (subtype in subtypes) {
    sheet_name <- diff_sheets[subtype]
    if (is.na(sheet_name) || !(sheet_name %in% getSheetNames(diff_expr_file))) {
      cat("Differential expression sheet missing for subtype:", subtype, "Skipping...\n")
      next
    }
    diff_df <- read.xlsx(diff_expr_file, sheet = sheet_name)

    if (!all(c("Gene", "label") %in% colnames(diff_df))) {
      cat("Required columns (Gene, label) missing in differential expression for:", subtype, "Skipping...\n")
      next
    }

    up_genes <- diff_df$Gene[diff_df$label == "up"]
    if (length(up_genes) == 0) {
      cat("No upregulated genes found for subtype:", subtype, "Skipping...\n")
      next
    }

    expr_file <- file.path(expr_base_path, paste0(subtype, "_expression.xlsx"))
    if (!file.exists(expr_file)) {
      cat("Expression file missing for subtype:", subtype, "at", expr_file, "\n")
      next
    }
    expr_df <- read_excel(expr_file)

    filtered_expr <- expr_df[expr_df[[gene_col_name]] %in% up_genes, ]
    if (nrow(filtered_expr) == 0) {
      cat("No matching genes in expression data for subtype:", subtype, "\n")
      next
    }

    output_file <- file.path(expr_base_path, paste0(subtype, "_up_expression.txt"))
    tryCatch({
      write.table(filtered_expr, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat(subtype, "upregulated gene expression data saved to:", output_file, "\n")
    }, error = function(e) {
      cat("Failed to save upregulated expression data to:", output_file, "Error:", e$message, "\n")
      next
    })
  }

  # Step 6: Load GDSC data
  cat("Loading GDSC data...\n")
  if (!file.exists(gdsc_expr_file) || !file.exists(gdsc_res_file)) {
    stop("GDSC data files missing: ", gdsc_expr_file, " or ", gdsc_res_file)
  }
  load(gdsc_expr_file)  # Loads GDSC2_Expr
  load(gdsc_res_file)   # Loads GDSC2_Res

  # Step 7: Drug response analysis
  cat("Starting drug response analysis...\n")

  processed_count <- 0
  total_combinations <- length(subtypes) * length(network_methods) * length(module_methods)

  # Iterate over network and module methods
  for (network in network_methods) {
    for (module_method in module_methods) {
      cat("=====================================\n")
      cat("Processing network:", network, "with module method:", module_method, "\n")
      cat("=====================================\n")

      # Process each subtype
      for (subtype in subtypes) {
        cat("Processing subtype:", subtype, "\n")

        # Create directory structure
        subtype_dir <- file.path(drug_response_path, subtype)
        dir.create(subtype_dir, showWarnings = FALSE, recursive = TRUE)

        module_dir <- file.path(subtype_dir, module_method)
        dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)

        network_dir <- file.path(module_dir, tolower(network))
        dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)
        Sys.sleep(0.2)  # Delay to ensure directory creation
        if (!dir.exists(network_dir)) {
          cat("Failed to create network directory:", network_dir, "Skipping...\n")
          next
        }
        if (file.access(network_dir, mode = 2) != 0) {
          cat("No write access to network directory:", network_dir, "Skipping...\n")
          next
        }
        cat("Network directory:", network_dir, "\n")

        # Read filtered expression data
        expr_file <- file.path(expr_base_path, paste0(subtype, "_up_expression.txt"))
        if (!file.exists(expr_file)) {
          cat("Filtered expression file missing for subtype:", subtype, "at", expr_file, "\n")
          next
        }
        expression_data <- read.delim(expr_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>%
          rename(node = all_of(gene_col_name)) %>%
          mutate_all(~ifelse(is.na(.), 0, .))
        expression_data$node <- as.character(expression_data$node)

        # Read module file
        module_file <- file.path(module_base_path, network, module_method, subtype,
                                 paste0("node_Module_select_", network, "_", subtype, ".txt"))
        if (!file.exists(module_file)) {
          cat("Module file missing for network:", network, "module method:", module_method,
              "subtype:", subtype, "at", module_file, "\n")
          next
        }
        module_data <- read.delim(module_file, header = TRUE, stringsAsFactors = FALSE)

        # Validate module data
        if (!all(c("node", "module") %in% colnames(module_data))) {
          cat("Module file missing required columns (node, module):", module_file, "\n")
          next
        }
        module_data$node <- as.character(module_data$node)
        module_data$module <- as.character(module_data$module)
        if (any(is.na(module_data$node) | module_data$node == "")) {
          cat("Module file contains NA or empty nodes:", module_file, "\n")
          next
        }
        if (any(is.na(module_data$module) | module_data$module == "")) {
          cat("Module file contains NA or empty modules:", module_file, "\n")
          next
        }
        non_key_cols <- setdiff(colnames(module_data), c("node", "module"))
        if (length(non_key_cols) > 0) {
          module_data[non_key_cols][is.na(module_data[non_key_cols])] <- 0
        }

        # Merge expression and module data
        tryCatch({
          merged_data <- module_data %>%
            left_join(expression_data, by = "node") %>%
            distinct(node, .keep_all = TRUE)
        }, error = function(e) {
          cat("Error in merging data for subtype:", subtype, "network:", network,
              "module method:", module_method, "Error:", e$message, "\n")
          next
        })

        if (!exists("merged_data") || nrow(merged_data) == 0) {
          cat("Merged data is empty or not created for subtype:", subtype, "\n")
          next
        }

        if (any(is.na(merged_data$module))) {
          cat("Merged data contains NA modules for subtype:", subtype, "\n")
          next
        }

        cat("Original module proteins:", nrow(module_data), "\n")
        cat("Matched proteins:", sum(module_data$node %in% expression_data$node), "\n")
        cat("Unique modules:", length(unique(merged_data$module)), "\n")

        # Save merged data
        merged_output <- file.path(network_dir,
                                   paste0("expr_file_", tolower(module_method), "_", tolower(network), "_", subtype, ".txt"))
        tryCatch({
          write.table(merged_data, merged_output, sep = "\t", row.names = FALSE, quote = FALSE)
          cat("Merged data saved:", merged_output, "\n")
        }, error = function(e) {
          cat("Error saving merged data to:", merged_output, "Error:", e$message, "\n")
          next
        })

        # Process each module
        unique_modules <- unique(merged_data$module)
        if (length(unique_modules) == 0) {
          cat("No valid modules found for subtype:", subtype, "network:", network,
              "module method:", module_method, "\n")
          next
        }

        for (current_module in unique_modules) {
          cat("Analyzing module:", current_module, "\n")

          # Create module analysis directory
          module_analysis_dir <- file.path(network_dir, current_module)
          dir.create(module_analysis_dir, showWarnings = FALSE, recursive = TRUE)
          Sys.sleep(0.2)
          if (!dir.exists(module_analysis_dir)) {
            cat("Failed to create module analysis directory:", module_analysis_dir, "Skipping...\n")
            next
          }
          if (file.access(module_analysis_dir, mode = 2) != 0) {
            cat("No write access to module analysis directory:", module_analysis_dir, "Skipping...\n")
            next
          }
          cat("Module analysis directory:", module_analysis_dir, "\n")

          # Save module-specific data
          subset_data <- merged_data %>% filter(module == current_module)
          if (nrow(subset_data) == 0 || ncol(subset_data) <= 1) {
            cat("No valid data for module:", current_module, "Skipping...\n")
            next
          }

          module_expr_file <- file.path(module_analysis_dir,
                                        paste0("expr_file_", tolower(module_method), "_", tolower(network), "_", subtype, "_", current_module, ".csv"))
          tryCatch({
            write.csv(subset_data, module_expr_file, row.names = FALSE)
            cat("Module data saved:", module_expr_file, "\n")
          }, error = function(e) {
            cat("Error saving module data to:", module_expr_file, "Error:", e$message, "\n")
            next
          })

          # Read module expression data
          if (!file.exists(module_expr_file)) {
            cat("Module expression file missing:", module_expr_file, "Skipping module\n")
            next
          }
          testData <- read.csv(module_expr_file, header = TRUE, row.names = 1)
          if ("module" %in% colnames(testData)) {
            testData <- testData[, !colnames(testData) %in% c("module"), drop = FALSE]
          }
          colnames(testData) <- make.names(colnames(testData), unique = TRUE)
          testExpr <- as.matrix(testData)

          # Validate testExpr
          min_samples <- 5
          if (nrow(testExpr) == 0 || ncol(testExpr) < min_samples || any(is.na(testExpr))) {
            cat("Invalid testExpr for module:", current_module,
                "Rows:", nrow(testExpr), "Cols:", ncol(testExpr),
                "NA present:", any(is.na(testExpr)), "Skipping...\n")
            next
          }
          if (!all(sapply(testExpr, is.numeric))) {
            cat("Non-numeric values in testExpr for module:", current_module,
                "Columns:", paste(colnames(testExpr)[!sapply(testExpr, is.numeric)], collapse = ", "),
                "Skipping...\n")
            next
          }

          # Run oncoPredict
          calcPhenotype_output_dir <- file.path(module_analysis_dir, "calcPhenotype_Output")
          dir.create(calcPhenotype_output_dir, showWarnings = FALSE, recursive = TRUE)
          if (file.access(calcPhenotype_output_dir, mode = 2) != 0) {
            cat("No write access to calcPhenotype output directory:", calcPhenotype_output_dir, "Skipping...\n")
            next
          }

          original_wd <- getwd()
          success <- FALSE

          tryCatch({
            setwd(module_analysis_dir)
            calcPhenotype(
              trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = min_samples,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData',
              cc = TRUE
            )
            temp_drug_pred_file <- file.path(module_analysis_dir, "DrugPredictions.csv")
            if (file.exists(temp_drug_pred_file)) {
              file.rename(temp_drug_pred_file, file.path(calcPhenotype_output_dir, "DrugPredictions.csv"))
            }
            success <- TRUE
          }, error = function(e) {
            cat("Error in oncoPredict for module:", current_module, "Error:", e$message, "\n")
          }, finally = {
            setwd(original_wd)
          })

          if (!success) {
            cat("Skipping module", current_module, "due to calcPhenotype failure\n")
            next
          }

          # Check for DrugPredictions.csv
          drug_pred_file <- file.path(calcPhenotype_output_dir, "DrugPredictions.csv")
          if (!file.exists(drug_pred_file)) {
            cat("Drug predictions missing for module:", current_module, "at:", drug_pred_file, "\n")
            next
          }

          # Prepare data for clustering
          numeric_cols <- sapply(testData, is.numeric)
          testData_numeric <- testData[, numeric_cols, drop = FALSE]
          moduleData <- testData_numeric[, colSums(testData_numeric, na.rm = TRUE) != 0, drop = FALSE]
          moduleData <- as.matrix(moduleData)

          # Validate moduleData
          if (nrow(moduleData) == 0 || ncol(moduleData) == 0) {
            cat("Invalid moduleData for module:", current_module,
                "Rows:", nrow(moduleData), "Cols:", ncol(moduleData), "Skipping...\n")
            next
          }

          # Run ConsensusClusterPlus
          cluster_output_dir <- file.path(module_analysis_dir, "sample_cluster")
          dir.create(cluster_output_dir, showWarnings = FALSE, recursive = TRUE)
          if (file.access(cluster_output_dir, mode = 2) != 0) {
            cat("No write access to cluster output directory:", cluster_output_dir, "Skipping...\n")
            next
          }

          success <- FALSE
          tryCatch({
            setwd(cluster_output_dir)
            ConsensusClusterPlus(
              moduleData,
              maxK = 3,
              reps = 1000,
              pItem = 0.8,
              pFeature = 1,
              seed = 10000,
              clusterAlg = "pam",
              distance = "pearson",
              title = "cluster_results",
              plot = "png",
              writeTable = TRUE
            )
            temp_subdir <- file.path(cluster_output_dir, "cluster_results")
            if (dir.exists(temp_subdir)) {
              temp_files <- list.files(temp_subdir, full.names = TRUE)
              for (temp_file in temp_files) {
                target_file <- file.path(cluster_output_dir, basename(temp_file))
                file.rename(temp_file, target_file)
              }
              temp_cluster_file <- file.path(cluster_output_dir, "cluster_results.2.consensusClass.csv")
              target_cluster_file <- file.path(cluster_output_dir, "cluster_results.k=2.consensusClass.csv")
              if (file.exists(temp_cluster_file)) {
                file.rename(temp_cluster_file, target_cluster_file)
              }
              unlink(temp_subdir, recursive = TRUE)
            }
            success <- TRUE
          }, error = function(e) {
            cat("Error in ConsensusClusterPlus for module:", current_module, "Error:", e$message, "\n")
          }, finally = {
            setwd(original_wd)
          })

          if (!success) {
            cat("Skipping module", current_module, "due to ConsensusClusterPlus failure\n")
            next
          }

          # Check for clustering results
          cluster_file <- file.path(cluster_output_dir, "cluster_results.k=2.consensusClass.csv")
          if (!file.exists(cluster_file)) {
            cat("Clustering result missing for module:", current_module, "at:", cluster_file, "\n")
            next
          }

          # Read drug predictions and clustering results
          drugData <- read.csv(drug_pred_file, check.names = FALSE)
          colnames(drugData)[1] <- "Sample"

          cluster <- read.csv(cluster_file, header = FALSE)
          colnames(cluster) <- c("Sample", "Cluster")

          # Merge data
          dat <- merge(cluster, drugData, by = "Sample")
          drug_names <- colnames(dat)[3:ncol(dat)]

          # Perform Wilcoxon test for drug response differences
          pvalue <- c()
          label <- c()
          for (i in 3:ncol(dat)) {
            group1 <- dat[dat$Cluster == 1, i]
            group2 <- dat[dat$Cluster == 2, i]
            res <- wilcox.test(group1, group2)
            pvalue <- c(pvalue, res$p.value)
            label <- c(label, ifelse(res$p.value < 0.05, "sig", "non-sig"))
          }

          drugSig <- data.frame(drug = drug_names, pvalue = pvalue, label = label)
          countSig <- c(sum(drugSig$label == "non-sig"), sum(drugSig$label == "sig"))
          sigCount <- data.frame(label = c("non-sig", "sig"), count = countSig)

          # Save drug response level
          drug_level_file <- file.path(module_analysis_dir,
                                       paste0("drug_response_level_", tolower(module_method), "_", tolower(network), "_", subtype, "_", current_module, ".csv"))
          tryCatch({
            write.csv(sigCount, drug_level_file, row.names = FALSE, quote = FALSE)
            cat("Drug response level saved:", drug_level_file, "\n")
          }, error = function(e) {
            cat("Error saving drug response level to:", drug_level_file, "Error:", e$message, "\n")
          })

          # Generate Drug Response Network (DRN)
          moduleData_for_DRN <- data.frame(t(testData))
          moduleLabel <- apply(moduleData_for_DRN, 2, function(x) {
            ifelse(x > median(x, na.rm = TRUE), "high", "low")
          })
          moduleLabel <- as.data.frame(moduleLabel)

          drugData2 <- read.csv(drug_pred_file, row.names = 1, check.names = FALSE)

          dat_DRN <- merge(moduleLabel, drugData2, by = "row.names")
          rownames(dat_DRN) <- dat_DRN$Row.names
          dat_DRN <- dat_DRN[, -1]

          protein_names <- colnames(moduleLabel)
          drug_names2 <- colnames(drugData2)

          proteins_vec <- c()
          drugs_vec <- c()
          pvalue_vec <- c()

          for (pro in protein_names) {
            for (dg in drug_names2) {
              proteins_vec <- c(proteins_vec, pro)
              drugs_vec <- c(drugs_vec, dg)
              temp <- dat_DRN[, c(pro, dg)]
              dat_high <- temp[temp[,1] == "high", 2]
              dat_low <- temp[temp[,1] == "low", 2]
              if (all(is.na(dat_high)) || all(is.na(dat_low)) || length(dat_high) < 2 || length(dat_low) < 2) {
                pvalue_vec <- c(pvalue_vec, NA)
              } else {
                res_temp <- wilcox.test(dat_high, dat_low)
                pvalue_vec <- c(pvalue_vec, res_temp$p.value)
              }
            }
          }

          DPI_all <- data.frame(protein = proteins_vec, drug = drugs_vec, pvalue = pvalue_vec)
          DPI_all$label <- ifelse(DPI_all$pvalue < 0.05, "sig", "non-sig")
          DRN <- subset(DPI_all, pvalue < 0.05)

          # Save DRN
          drn_file <- file.path(module_analysis_dir,
                                paste0("DRN_", tolower(module_method), "_", tolower(network), "_", subtype, "_", current_module, ".txt"))
          tryCatch({
            write.table(DRN, drn_file, sep = "\t", quote = FALSE, row.names = FALSE)
            cat("DRN saved:", drn_file, "\n")
          }, error = function(e) {
            cat("Error saving DRN to:", drn_file, "Error:", e$message, "\n")
          })

          # Save DRN info
          drn_info_file <- file.path(module_analysis_dir,
                                     paste0("DRN_info_", tolower(module_method), "_", tolower(network), "_", subtype, "_", current_module, ".txt"))
          tryCatch({
            DRN.info <- data.frame(node = c(protein_names, drug_names2),
                                   type = c(rep("protein", length(protein_names)),
                                            rep("drug", length(drug_names2))))
            write.table(DRN.info, drn_info_file, sep = "\t", quote = FALSE, row.names = FALSE)
            cat("DRN_info saved:", drn_info_file, "\n")
          }, error = function(e) {
            cat("Error saving DRN_info to:", drn_info_file, "Error:", e$message, "\n")
          })

          cat("Drug response analysis completed for module:", current_module, "\n")
          processed_count <- processed_count + 1
        }

        cat("All modules for subtype", subtype, "completed!\n")
      }
    }
  }

  # Step 8: Summarize processed combinations
  cat("Processed", processed_count, "of", total_combinations * sum(sapply(subtypes, function(s) {
    sum(sapply(module_methods, function(mm) {
      sum(sapply(network_methods, function(n) {
        module_file <- file.path(module_base_path, n, mm, s, paste0("node_Module_select_", n, "_", s, ".txt"))
        if (file.exists(module_file)) {
          modules <- unique(read.delim(module_file)$module)
          length(modules)
        } else {
          0
        }
      }))
    }))
  })), "module combinations\n")
}
