#' @title Generate Protein and Drug Sequence Files for Drug Response Prediction
#' @description Processes drug response network data to generate protein sequence (FASTA) and
#' drug SMILES files required for DeepPurpose-based drug-target interaction prediction.
#' This function serves as an R wrapper for Python sequence processing functionality.
#' @importFrom reticulate use_condaenv py_run_string import_main py_to_r
#' @param input_base Character. Path to directory containing Drug Response Network (DRN) input files. Each file should be named in the format "*_DRN_info.txt".
#' @param output_base Character. Path to directory where output sequence files will be saved.
#' @param subtype_file Character. Path to Excel file containing subtype phenotype data.
#' @param py_env Character. Name or path of Conda environment containing required Python  dependencies (default: "py3.9").
#' @return Saves sequence output directory.
#' @details DETAILS
#' @examples
#' \dontrun{
#' run_SEQCre(
#'   input_base = "data/DRN_files",
#'   output_base = "output/sequence_files",
#'   subtype_file = "metadata/subtype_info.xlsx",
#'   py_env = "drug_discovery"
#' )
#' }
#' @rdname get_seq
#' @export
run_SEQCre <- function(input_base,
                         output_base,
                         subtype_file,
                         py_env = "py3.9") {
  library(reticulate)

  # 2) Switch to the specified conda environment
  use_condaenv(py_env, required = TRUE)

  # 3) Find and load SEQCre.py from outside the package or in the project directory. 
  py_path <- "F:/sample_test/python/SEQCre/7.SEQCre.py"  # Adjust according to your actual path.
  if (!file.exists(py_path)) {
    stop("cannot found SEQCre.py，please confirm py_path is correct：", py_path)
  }
  source_python(py_path)

  # 4) Call the Python interface
  #  run_seqcre will traverse all DRN_info files under input_base, generating the corresponding seq_file and smiles_file to output_base.
  py$run_seqcre(input_base, output_base, subtype_file)

  invisible(NULL)
}
