#' @title Split Edge Module Data
#' @description Splits edge module data into separate files based on module numbers and saves them to specified directories.
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv write_tsv read_delim write_delim read_csv cols
#' @importFrom readxl read_excel
#' @importFrom dplyr rename pull filter select mutate arrange inner_join row_number distinct
#' @importFrom stringr str_extract str_replace
#' @importFrom reticulate use_condaenv py_run_string import_main py_to_r
#' @importFrom purrr walk
#' @param input_file  Character string specifying the path to the input edge module file.
#' @param out_base Character string specifying the base directory for output files.
#' @return Saves edge module files for each valid module to the specified output directory; returns NULL if no valid modules.
#' @details DETAILS
#' @examples
#' \dontrun{
#' files <- list.files(
#'     path       = "./ModuleDivision",
#'     pattern    = "^edge_Module_.*\\.txt$",
#'     recursive  = TRUE,
#'     full.names = TRUE
#' )
#' purrr::walk(files, split_edge_module,out_base)
#' }
#' @rdname split_edge_module
#' @export
split_edge_module <- function(
    input_file,
    out_base
) {

  # Step 1: Parse path information to extract metadata
  file_dir       <- dirname(input_file)
  subtype        <- basename(file_dir)
  module_method  <- basename(dirname(file_dir))
  network_method <- basename(dirname(dirname(file_dir)))

  # Step 2: Read input data
  df <- read_tsv(input_file, col_types = cols())

  # Step 3: Rename node1/node2 to gene1/gene2 if present
  if (all(c("node1","node2") %in% names(df))) {
    df <- df %>% rename(gene1 = node1, gene2 = node2)
  }
  if (!"module" %in% names(df)) {
    stop("Input file missing 'module' column: ", input_file)
  }

  # Step 4: Filter modules greater than 0
  modules <- df %>% pull(module) %>% unique() %>% sort()
  modules <- modules[modules > 0]
  if (length(modules) == 0) {
    message("No modules greater than 0, skipping: ", input_file)
    return(invisible(NULL))
  }

  # Step 5: Split data by module and write to output files
  purrr::walk(modules, function(mod) {
    # Construct output directory path
    out_dir <- file.path(
      out_base,
      subtype,
      module_method,
      tolower(network_method),
      paste0("M", mod)
    )

    # Check if output directory exists
    if (!dir.exists(out_dir)) {
      message("Output directory does not exist, skipping module M", mod, ": ", out_dir)
      return(invisible(NULL))
    }

    # Filter data for the current module
    sub_df <- df %>%
      filter(module == mod) %>%
      select(gene1, gene2)

    # Construct output file path
    out_file <- file.path(
      out_dir,
      paste0("edge_Module_",
             network_method, "_",
             subtype, "_M", mod, ".txt")
    )

    # Write filtered data to file
    write_tsv(sub_df, out_file)
    message("Saved: ", out_file)
  })
}

#' @title Process PRS-DTI Data
#' @description Perturbation Response Scanning (PRS) and Drug-Target Interaction (DTI) data using DeepPurpose and compute protein-drug interaction scores for each subtype, network, and module.
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv write_tsv read_csv write_csv cols
#' @importFrom readxl read_excel
#' @importFrom dplyr rename pull filter select mutate arrange inner_join first
#' @importFrom reticulate use_condaenv py
#' @importFrom utils write.csv
#' @importFrom stats na.omit
#' @param subtype_file Character string specifying the path to the subtype phenotype data file.
#' @param prs_ps_base Character string specifying the base directory containing PRS and DTI files.
#' @param network_methods Character vector specifying the PPI network databases (default: c("string", "physicalppin", "chengf")).
#' @param module_methods Character vector specifying the module detection methods (default: c("Louvain", "WF")).
#' @param py_env  Character string specifying the path to the Python Conda environment with DeepPurpose installed.
#' @return Saves PRS-DTI scores and related files to the specified output directories.
#' @details DETAILS
#' @examples
#' \dontrun{
#'  process_prs_dti(
#'     subtype_file    = "./subtype.xlsx",
#'     prs_ps_base     = "./PRS/PS",
#'     network_methods = c("string","physicalppin","chengf"),
#'     module_methods  = c("Louvain","WF"),
#'     py_env          = "./envs/py3.9"
#'  )
#'  }
#' @rdname process_prs_dti
#' @export
process_prs_dti <- function(
    subtype_file,
    prs_ps_base,
    network_methods = c("string","physicalppin","chengf"),
    module_methods  = c("Louvain","WF"),
    py_env
) {
  # 1. Load and check parameters
  if (!file.exists(subtype_file)) stop("Cannot find subtype file: ", subtype_file)
  if (!dir.exists(prs_ps_base))    stop("The PRS PS base directory does not exist: ", prs_ps_base)

  # 2. Switch and initialize the Python environment
  # Switch and initialize the Python environment
  use_condaenv(py_env, required = TRUE)
  # Load DeepPurpose DTI
  py_run_string("from DeepPurpose import DTI as models")
  pretrained_model <- "MPNN_CNN_BindingDB_IC50"
  model <- py$models$model_pretrained(model = pretrained_model)
  # Set the local source code path
  enm_path <- "F:\\sample_test\\python"
  # Add the local path to Python's sys.path
  py_run_string(sprintf("import sys; sys.path.append(r'%s')", enm_path))
  # load enm
  py_run_string("from enm.Enm import *")

  # 3. Read subtype list
  subtypes <- read_excel(subtype_file, skip=1) %>% pull(2) %>% unique() %>% na.omit()
  if (length(subtypes)==0) stop("No subtypes detected, please check the subtype file.")

  # 4. Traverse various subtypes/networks/modular approaches
  for (phenotype in subtypes) {
    cat("\n==== Handling subtypes: ", phenotype, "====\n")
    for (net in network_methods) for (modl in module_methods) {
      base_dir <- file.path(prs_ps_base, phenotype, modl, net)
      if (!dir.exists(base_dir)) {
        cat("Skip: Directory does not exist", base_dir, "\n"); next
      }
      # Each module subdirectory
      module_dirs <- list.dirs(base_dir, full.names=TRUE, recursive=FALSE)
      if (length(module_dirs)==0) {
        cat("Not in", base_dir, "Find the module directory\n"); next
      }
      for (current_dir in module_dirs) {
        mod_name <- basename(current_dir)
        cat("-> Module：", mod_name, "\n")
        BA_pattern   <- paste0("^virtual_screening__", modl, "_", net, "_", phenotype, "_", mod_name, "\\.txt$")
        PPIN_pattern <- "^edge_Module.*\\.txt$"
        BA_file   <- list.files(current_dir, BA_pattern, full.names=TRUE) %>% first()
        PPIN_file <- list.files(current_dir, PPIN_pattern, full.names=TRUE) %>% first()
        if (is.na(BA_file) || is.na(PPIN_file)) {
          cat("  The document is incomplete, skipping the module.\n"); next
        }

        # 5. Network Analysis
        ok <- try({
          py_run_string("enm = Enm('PPIN')")
          py_run_string(sprintf("enm.read_network(r'%s', sep='\\t')", PPIN_file))
          py_run_string("enm.gnm_analysis(normalized=False)")
          py_run_string("enm.cluster_matrix(enm.prs_mat)")
        }, silent=TRUE)
        if (inherits(ok, "try-error")) {
          cat("  Python analysis failed, skipping\n"); next
        }
        #  sensitivity
        pcc_df_path     <- file.path(current_dir, "pcc_df.csv")
        prs_mat_df_path <- file.path(current_dir, "prs_mat_df.txt")
        py_run_string(sprintf("enm.df.to_csv(r'%s', index=True)", pcc_df_path))
        py_run_string(sprintf(
          "enm.prs_mat_df.to_csv(r'%s', sep='\\t', quoting=3)",
          prs_mat_df_path
        ))

        # 6. Read BA
        BA <- read_tsv(BA_file, col_types=cols()) %>%
          rename(
            Drug.Name    = `Drug Name`,
            Target.Name  = `Target Name`,
            Binding.Score= `Binding Score`
          ) %>%
          select(Rank, Target.Name, Drug.Name, Binding.Score) %>%
          mutate(
            Binding.Score = as.numeric(Binding.Score)
          ) %>%
          filter(!is.na(Binding.Score))

        # 7. Read sensitivity
        Sens <- read_csv(pcc_df_path, col_types=cols()) %>%
          rename(Target.Name = orf_name) %>%
          select(Target.Name, sens) %>%
          mutate(sens = as.numeric(sens)) %>%
          filter(!is.na(sens))

        # 8. Merge calculation ps
        ps_df <- BA %>%
          inner_join(Sens, by="Target.Name") %>%
          mutate(ps = Binding.Score * sens) %>%
          arrange(desc(ps))

        # 9. save results
        write_csv(ps_df, file.path(current_dir, "prs_dti_score.csv"))
        ps_df %>%
          select(Drug.Name, Target.Name, score = ps) %>%
          write_tsv(file.path(current_dir, "dpi_ps_score.txt"))
        cat("  Module processing completed：", current_dir, "\n")
      }
    }
    cat("Subtype completed：", phenotype, "\n")
  }
  cat("\nAll processing completed!\n")
}
