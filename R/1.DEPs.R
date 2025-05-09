#' @title Differential Expression Analysis Function
#' @description This function performs differential expression analysis on gene expression data across different subtypes, using Wilcoxon rank-sum tests and BH correction.
#' @importFrom readxl read_excel
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom stats wilcox.test p.adjust
#' @param input_dir Character string specifying the input directory containing expression and phenotype data.
#' @param expr_name Character string for the expression data file name (default: "expression.xlsx").
#' @param pheno_name Character string for the phenotype data file name (default: "subtype.xlsx").
#' @param detection_threshold Numeric threshold for gene detection rate (default: 0.75).
#' @param fc_threshold Numeric threshold for fold change (default: 1.5).
#' @param p_threshold Numeric threshold for adjusted p-value (default: 0.05).
#' @param out_excel Character string for the output Excel file name (default: "diff_expression_results.xlsx").
#' @return Saves differential expression results and statistics to an Excel file in the output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' run_diff_expr_analysis(
#'   input_dir          = "your_workspace",
#'   expr_name          = "expression.xlsx",
#'   pheno_name         = "subtype.xlsx",
#'   detection_threshold = 0.8,
#'   fc_threshold       = 2,
#'   p_threshold        = 0.01
#' )
#' or
#' run_diff_expr_analysis("your file path")}
#' @rdname run_diff_expr_analysis
#' @export
run_diff_expr_analysis <- function(input_dir,
                                   expr_name = "expression.xlsx",
                                   pheno_name = "subtype.xlsx",
                                   detection_threshold = 0.75,
                                   fc_threshold = 1.5,
                                   p_threshold = 0.05,
                                   out_excel = "diff_expression_results.xlsx") {

  # Create output directory for results
  out_dir <- file.path(input_dir, "DEPS_result")
  dir.create(out_dir, showWarnings = FALSE)  # Suppress warnings if directory already exists

  # Read expression data from Excel file
  expr_path <- file.path(input_dir, expr_name)
  expr_data <- read_excel(expr_path)
  expr_data <- as.data.frame(expr_data)  # Convert to data frame
  rownames(expr_data) <- expr_data[[1]]  # Set first column as row names (gene IDs)
  expr_data <- expr_data[-1]             # Remove the gene ID column from data

  # Read phenotype data from Excel file
  pheno_path <- file.path(input_dir, pheno_name)
  subtype_data <- as.data.frame(read_excel(pheno_path))  # Convert to data frame

  # Extract unique subtypes from phenotype data
  subtypes <- unique(subtype_data$Subtype)
  results_list <- list()  # Initialize list to store differential expression results
  stats_list <- list()    # Initialize list to store label statistics

  # Internal Function: Perform differential expression analysis for a single subtype
  # Parameters:
  #   group1: Data frame of expression data for the experimental group (specific subtype).
  #   group2: Data frame of expression data for the control group (non-specific subtype).
  #   detection_threshold: Numeric threshold for gene detection rate.
  #   fc_threshold: Numeric threshold for fold change.
  #   p_threshold: Numeric threshold for adjusted p-value.
  # Returns: Data frame containing differential expression results for the subtype.
  diff_expr_test_single <- function(group1, group2, detection_threshold = 0.75,
                                    fc_threshold = 1.5, p_threshold = 0.05) {
    n_genes <- nrow(group1)  # Number of genes in the expression data

    # Initialize vectors to store results
    mean1  <- rep(NA, n_genes)  # Mean expression for experimental group
    mean2  <- rep(NA, n_genes)  # Mean expression for control group
    fc     <- rep(NA, n_genes)  # Fold change
    log2FC <- rep(NA, n_genes)  # Log2 fold change
    pvals  <- rep(NA, n_genes)  # P-values from Wilcoxon test
    detect1 <- rep(NA, n_genes) # Detection rate for experimental group
    detect2 <- rep(NA, n_genes) # Detection rate for control group

    # Loop through each gene to compute statistics
    for (i in 1:n_genes) {
      gene_expr1 <- as.numeric(group1[i, ])  # Expression values for experimental group
      gene_expr2 <- as.numeric(group2[i, ])  # Expression values for control group

      # Calculate detection rate: proportion of non-NA and non-zero values
      detect1[i] <- sum(!is.na(gene_expr1) & gene_expr1 != 0) / length(gene_expr1)
      detect2[i] <- sum(!is.na(gene_expr2) & gene_expr2 != 0) / length(gene_expr2)

      # Skip gene if detection rate in either group is below threshold
      if (detect1[i] < detection_threshold || detect2[i] < detection_threshold) {
        next
      }

      # Calculate mean expression for each group, ignoring NA values
      mean1[i] <- mean(gene_expr1, na.rm = TRUE)
      mean2[i] <- mean(gene_expr2, na.rm = TRUE)

      # Calculate fold change and log2 fold change, handling division by zero
      if (mean2[i] == 0) {
        fc[i] <- NA
        log2FC[i] <- NA
      } else {
        fc[i] <- mean1[i] / mean2[i]
        log2FC[i] <- log2(mean1[i]) - log2(mean2[i])
      }

      # Perform Wilcoxon rank-sum test
      test_res <- wilcox.test(gene_expr1, gene_expr2)
      pvals[i] <- test_res$p.value
    }

    # Apply Benjamini-Hochberg (BH) correction for multiple testing
    adj_p <- p.adjust(pvals, method = "BH")

    # Create results data frame
    res <- data.frame(
      Gene = rownames(group1),          # Gene IDs
      mean_group1 = mean1,             # Mean expression in experimental group
      mean_group2 = mean2,             # Mean expression in control group
      fc = fc,                         # Fold change
      log2FC = log2FC,                 # Log2 fold change
      p_value = pvals,                 # Raw p-value
      adj_p = adj_p,                   # Adjusted p-value
      detection_group1 = detect1,       # Detection rate in experimental group
      detection_group2 = detect2,       # Detection rate in control group
      stringsAsFactors = FALSE
    )

    # Assign labels based on significance and fold change
    res$label <- "non-significant"  # Default label
    # Label as "up" for upregulated genes
    res$label[
      res$detection_group1 >= detection_threshold &
        res$detection_group2 >= detection_threshold &
        res$adj_p < p_threshold & !is.na(res$fc) & res$fc > fc_threshold
    ] <- "up"
    # Label as "down" for downregulated genes
    res$label[
      res$detection_group1 >= detection_threshold &
        res$detection_group2 >= detection_threshold &
        res$adj_p < p_threshold & !is.na(res$fc) & res$fc < (1 / fc_threshold)
    ] <- "down"

    # Sort results by adjusted p-value
    res <- res[order(res$adj_p), ]
    return(res)
  }

  # Loop through each subtype to perform differential expression analysis
  for (g in subtypes) {
    # Extract samples for the current subtype (experimental) and others (control)
    exp_samples <- subtype_data$Sample[subtype_data$Subtype == g]
    ctrl_samples <- subtype_data$Sample[subtype_data$Subtype != g]

    # Subset expression data for experimental and control groups
    group_exp <- expr_data[, colnames(expr_data) %in% exp_samples]
    group_ctrl <- expr_data[, colnames(expr_data) %in% ctrl_samples]

    # Run differential expression analysis for the current subtype
    res <- diff_expr_test_single(group_exp, group_ctrl)

    # Store results in the results list
    results_list[[g]] <- res

    # Calculate and store label statistics (e.g., counts of up/down/non-significant)
    stats <- as.data.frame(table(res$label))
    colnames(stats) <- c("Label", "Count")
    stats_list[[g]] <- stats
  }

  # Save results to an Excel file with separate sheets for each subtype
  output_file <- file.path(out_dir, "diff_expression_results_all.xlsx")
  wb <- createWorkbook()  # Create a new workbook

  # Loop through subtypes to write results and statistics to Excel sheets
  for (g in subtypes) {
    sheet_diff <- paste0(g, "_DiffResults")   # Sheet for differential results
    sheet_stats <- paste0(g, "_Statistics")   # Sheet for label statistics

    # Add and write differential expression results
    addWorksheet(wb, sheet_diff)
    writeData(wb, sheet_diff, results_list[[g]])

    # Add and write label statistics
    addWorksheet(wb, sheet_stats)
    writeData(wb, sheet_stats, stats_list[[g]])
  }

  # Save the workbook to the output file
  saveWorkbook(wb, file = output_file, overwrite = TRUE)

  # Print completion message with output file path
  cat("All subtype differential expression analyses completed. Results saved to:", output_file, "\n")
}

