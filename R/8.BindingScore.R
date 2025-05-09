#' @title Drug Response Prediction Using DeepPurpose
#' @description This function predicts drug-target interactions using the DeepPurpose Python library for specified subtypes, network methods, and module methods, integrating DRN and SEQ data.
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
#' @importFrom dplyr pull mutate select arrange rename row_number inner_join
#' @importFrom stringr str_extract str_replace
#' @importFrom readr read_delim write_delim cols
#' @importFrom reticulate use_condaenv py_run_string import_main py_to_r
#' @param subtype_file Character string specifying the path to the subtype phenotype data file.
#' @param drn_base Character string specifying the base directory containing Drug Response Network (DRN) files.
#' @param prs_base Character string specifying the base directory containing Protein and Drug sequence files.
#' @param output_base Character string specifying the base directory for output files.
#' @param py_env Character string specifying the path to the Python Conda environment with DeepPurpose installed.
#' @param network_methods  Character vector specifying the PPI network databases (e.g., "string", "physicalppin", "chengf").
#' @param module_methods Character vector specifying the module detection methods (e.g., "Louvain", "WF").
#' @return Saves DPI information and drug-target interaction predictions to the specified output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  predict_drug_response(
#'  subtype_file = "./subtype.xlsx",
#'  drn_base = "./DrugResponse",
#'  prs_base = "./PRS/SMILES_SEQ",
#'  output_base = "./PRS/PS",
#'  py_env = "./envs/py3.9",
#'  network_methods = c("string", "physicalppin", "chengf"),
#'  module_methods = c("Louvain", "WF")
#'  )
#'  }
#' }
#' @rdname predict_drug_response
#' @export
predict_drug_response <- function(
    subtype_file,
    drn_base,
    prs_base,
    output_base,
    py_env,
    network_methods,
    module_methods
) {

  # 1. Read subtype file
  if (!file.exists(subtype_file)) {
    stop("subtype.xlsx does not exist: ", subtype_file)
  }
  subtypes <- read_excel(subtype_file, skip = 1) %>%
    pull(2) %>%
    unique() %>%
    na.omit()

  if (length(subtypes) == 0) {
    stop("No valid subtype information found in subtype.xlsx file")
  }

  # 2. Set up Python environment
  use_condaenv(py_env)
  py_run_string("from DeepPurpose import DTI as models")
  py <- import_main()
  pretrained_model <- "MPNN_CNN_BindingDB_IC50"
  model <- py$models$model_pretrained(model = pretrained_model)

  # 3. Iterate over all subtypes, modules, and network methods
  for (phenotype in subtypes) {
    for (module_dir in module_methods) {
      for (network_method in network_methods) {
        module_dir_lower <- tolower(module_dir)
        network_method_lower <- tolower(network_method)

        # Get DRN file list
        drn_pattern <- paste0(
          "DRN_", module_dir_lower, "_", network_method_lower,
          "_", phenotype, "_M\\d+\\.txt$"
        )
        drn_files <- list.files(
          path = file.path(drn_base, phenotype, module_dir, network_method),
          pattern = drn_pattern,
          recursive = TRUE,
          full.names = TRUE
        )

        if (length(drn_files) == 0) {
          cat("No DRN files found at path:", file.path(drn_base, phenotype, module_dir, network_method), "\n")
          next
        }

        # Process each DRN file
        for (drn_file in drn_files) {
          cat("Processing:", drn_file, "\n")

          # Extract module name
          module <- str_extract(basename(drn_file), "M\\d+")

          # Read DRN data
          drn_df <- read_delim(drn_file, delim = "\t", col_types = cols())
          colnames(drn_df) <- tolower(colnames(drn_df))
          if (!all(c("protein", "drug") %in% colnames(drn_df))) {
            cat("DRN file missing 'protein' or 'drug' column\n")
            next
          }
          drn_df <- drn_df %>%
            mutate(drug = str_replace(drug, "_\\d+$", "")) %>%
            select(protein, drug)

          # Construct SEQ file paths
          seq_file <- file.path(
            prs_base, phenotype, module_dir, network_method, module,
            paste0(
              "seq_file_", module_dir_lower, "_", network_method_lower,
              "_", phenotype, "_", module, ".csv"
            )
          )
          smiles_file <- file.path(
            prs_base, phenotype, module_dir, network_method, module,
            paste0(
              "smiles_file_", module_dir_lower, "_", network_method_lower,
              "_", phenotype, "_", module, ".csv"
            )
          )

          if (!file.exists(seq_file) || !file.exists(smiles_file)) {
            cat("Missing SEQ data files\n")
            next
          }

          # Read sequence data
          seq_df <- read_csv(seq_file, show_col_types = FALSE) %>%
            rename(protein = node, sequence = sequence) %>%
            distinct(protein, .keep_all = TRUE) %>%
            rename(Sequence = sequence)

          smiles_df <- read_csv(smiles_file, show_col_types = FALSE) %>%
            rename(drug = node, SMILES = SMILES) %>%
            distinct(drug, .keep_all = TRUE)

          # Integrate data
          dpi_info <- drn_df %>%
            left_join(seq_df, by = "protein") %>%
            left_join(smiles_df, by = "drug") %>%
            rename(
              target_name = protein,
              drug_name = drug,
              target_seq = Sequence,
              drug_smiles = SMILES
            ) %>%
            filter(!is.na(target_seq) & !is.na(drug_smiles))

          # Create output directory
          current_out_dir <- file.path(output_base, phenotype, module_dir, network_method, module)
          dir.create(current_out_dir, recursive = TRUE, showWarnings = FALSE)

          # Save DPI information
          dpi_info_file <- file.path(
            current_out_dir,
            paste0(
              "DPI_info_", module_dir, "_", network_method, "_",
              phenotype, "_", module, ".txt"
            )
          )
          write_delim(dpi_info, dpi_info_file, delim = "\t")

          # Call DeepPurpose for prediction
          df_pred <- read_delim(dpi_info_file, delim = "\t", col_types = cols())
          target_seq_list <- as.list(df_pred$target_seq)
          drug_smiles_list <- as.list(df_pred$drug_smiles)
          target_name_list <- as.list(df_pred$target_name)
          drug_name_list <- as.list(df_pred$drug_name)

          if (length(target_seq_list) == 0 || length(drug_smiles_list) == 0) {
            cat("Empty data, skipping prediction\n")
            next
          }

          # Perform virtual screening
          my_predict <- py$models$virtual_screening(
            drug_smiles_list,
            target_seq_list,
            model,
            drug_name_list,
            target_name_list,
            verbose = FALSE
          )

          # Convert and save prediction results
          pred_scores <- py_to_r(my_predict)
          pred_df <- data.frame(
            drug_name = unlist(drug_name_list),
            target_name = unlist(target_name_list),
            binding_score = pred_scores
          )

          pred_file <- file.path(
            current_out_dir,
            paste0(
              "virtual_screening_", module_dir, "_", network_method, "_",
              phenotype, "_", module, ".txt"
            )
          )
          write_delim(pred_df, pred_file, delim = "\t")

          Sys.sleep(1)
        } # end drn_file loop
      } # end network_method loop
    } # end module_dir loop
  } # end phenotype loop

  # Process result files
  prs_ps_base <- output_base
  dpi_info_files <- list.files(prs_ps_base, pattern = "^DPI_info_.*\\.txt$", recursive = TRUE, full.names = TRUE)

  for (dpi_file in dpi_info_files) {
    cat("Processing DPI_info file:", dpi_file, "\n")

    current_dir <- dirname(dpi_file)

    # Construct exact vs_file name
    vs_file_name <- gsub("DPI_info", "virtual_screening", basename(dpi_file))
    vs_file <- file.path(current_dir, vs_file_name)

    if (!file.exists(vs_file)) {
      cat("vs_file does not exist:", vs_file, "\n")
      next
    }

    # Read files
    vs_df <- read_delim(vs_file, delim = "\t", col_types = cols())
    cat("vs_df columns:", names(vs_df), "\n")

    dpi_df <- read_delim(dpi_file, delim = "\t", col_types = cols())
    cat("dpi_df columns:", names(dpi_df), "\n")

    # Check for required columns
    if (!all(c("drug_name", "target_name") %in% names(vs_df))) {
      cat("Error: vs_df is missing required columns for", vs_file, "\n")
      next
    }
    if (!all(c("drug_name", "target_name") %in% names(dpi_df))) {
      cat("Error: dpi_df is missing required columns for", dpi_file, "\n")
      next
    }

    # Merge data
    merged_df <- dpi_df %>%
      inner_join(vs_df, by = c("drug_name", "target_name")) %>%
      arrange(binding_score) %>%  # Sort by Binding Score in ascending order (lower Binding Score ranks higher)
      mutate(Rank = row_number()) %>%
      select(Rank, target_name, drug_name, binding_score) %>%
      rename(
        'Target Name' = target_name,
        'Drug Name' = drug_name,
        'Binding Score' = binding_score
      )

    # Generate output file name
    output_file <- file.path(
      current_dir,
      gsub("^DPI_info_", "virtual_screening__", basename(dpi_file))
    )

    # Write output file
    write_delim(merged_df, output_file, delim = "\t")
  }
}
