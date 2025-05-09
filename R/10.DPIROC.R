#' @title DPI PS Score ROC Analysis
#' @description Perform ROC analysis on all dpi_ps_score.txt files under a specified directory, using known drug–target data to evaluate prediction performance across score cutoffs.
#' @importFrom utils read.table write.table
#' @importFrom pROC roc auc
#' @param dpi_base Character. Base directory in which to search for dpi_ps_score.txt files (recursively).
#' @param ttd_path Character. File path to the DTI_TTD.RData file (loads object TTD). If NULL, defaults to file.path(dpi_base, "DTI_TTD.RData").
#' @param drugbank_path Character. File path to the DTI_DrugBank.RData file
#'   (loads object ``DrugBank``). If NULL, defaults to file.path(dpi_base, "DTI_DrugBank.RData").
#' @return ROC_<dpi_ps_score.txt> file in the same directory containing cutoff, AUC, drug count, and DPI count.
#' @details DETAILS
#' @examples
#' \dontrun{
#' dpi_ps_roc(
#'     dpi_base      = "./PRS/PS",
#'     ttd_path      = "./PRS/PS/DTI_TTD.RData",
#'     drugbank_path = "./PRS/PS/DTI_DrugBank.RData"
#' }
#' @rdname dpi_ps_roc
#' @export
dpi_ps_roc <- function(
    dpi_base,
    ttd_path      = NULL,
    drugbank_path = NULL
) {

  # Set DTI file paths based on dpi_base (if not provided)
  if (is.null(ttd_path)) {
    ttd_path <- file.path(dpi_base, "DTI_TTD.RData")
  }
  if (is.null(drugbank_path)) {
    drugbank_path <- file.path(dpi_base, "DTI_DrugBank.RData")
  }

  # Load known drug-target data (DTI_TTD and DTI_DrugBank)
  load(ttd_path)
  load(drugbank_path)
  drug2target <- rbind(TTD, DrugBank)
  drug2target$DrugName <- toupper(drug2target$DrugName)
  drug2target <- unique(drug2target)

  # Find dpi_ps_score.txt files
  dpi_files <- list.files(
    path       = dpi_base,
    pattern    = "^dpi_ps_score\\.txt$",
    recursive  = TRUE,
    full.names = TRUE
  )
  if (length(dpi_files) == 0) {
    stop("No dpi_ps_score.txt files found. Please check the dpi_base path and file naming convention.")
  }

  # Perform ROC analysis for each dpi_ps_score.txt file
  for (file in dpi_files) {
    cat("Processing:", file, "\n")

    pred_dpi <- read.table(
      file,
      header          = TRUE,
      sep             = "\t",
      stringsAsFactors = FALSE
    )

    # Check if the file has at least 3 columns
    if (ncol(pred_dpi) < 3) {
      cat("File ", file, " has insufficient columns, skipping\n")
      next
    }

    # Rename the first 3 columns to "Drug.Name", "Target.Name", "score"
    colnames(pred_dpi)[1:3] <- c("Drug.Name", "Target.Name", "score")
    pred_dpi$Drug.Name <- toupper(pred_dpi$Drug.Name)
    pred_dpi$score <- as.numeric(pred_dpi$score)

    # Initialize vectors for ROC analysis results
    cutoff_ls <- c()
    auc_ls    <- c()
    drug_ls   <- c()
    dpi_ls    <- c()

    # Iterate over each score value as a candidate cutoff
    for (cutoff in pred_dpi$score) {
      temp <- pred_dpi[pred_dpi$score >= cutoff, ]

      # Assign type based on known drug-target relationships
      temp$type <- ifelse(
        temp$Drug.Name   %in% drug2target$DrugName &
          temp$Target.Name %in% drug2target$TargetName,
        1, 0
      )

      # Calculate ROC only if both high and low groups have data
      if (length(unique(temp$type)) == 2) {
        roc_obj <- roc(temp$type, temp$score)
        auc_val <- auc(roc_obj)

        cutoff_ls <- c(cutoff_ls, cutoff)
        auc_ls    <- c(auc_ls, auc_val)
        drug_ls   <- c(drug_ls, length(unique(temp$Drug.Name)))
        dpi_ls    <- c(dpi_ls, nrow(temp))
      }
    }

    # Create ROC analysis results data frame
    roc_df <- data.frame(
      cutoff    = cutoff_ls,
      auc       = auc_ls,
      drugCount = drug_ls,
      dpiCount  = dpi_ls,
      stringsAsFactors = FALSE
    )

    # Explicitly convert to numeric (in case of factors, etc.)
    roc_df$cutoff <- as.numeric(as.character(roc_df$cutoff))

    # Sort roc_df by cutoff in descending order
    roc_df <- roc_df[order(roc_df$cutoff, decreasing = TRUE), ]

    # Save results
    output_file <- file.path(dirname(file), paste0("ROC_", basename(file)))
    write.table(
      roc_df,
      file      = output_file,
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
    cat("Saved ROC file:", output_file, "\n")
  }
}
