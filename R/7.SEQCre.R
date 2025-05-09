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

  # 2) 切换到指定 conda 环境
  use_condaenv(py_env, required = TRUE)

  # 3) 找到并加载包外或项目目录下的 SEQCre.py
  #    这里假设你在当前工作目录或包的 inst/python 下放了脚本
  py_path <- "F:/sample_test/python/SEQCre/7.SEQCre.py"  # <- 根据你的实际路径调整
  if (!file.exists(py_path)) {
    stop("找不到 SEQCre.py，请确认 py_path 正确：", py_path)
  }
  source_python(py_path)

  # 4) 调用 Python 接口
  #    run_seqcre 会遍历 input_base 下的所有 DRN_info 文件，
  #    生成对应的 seq_file 和 smiles_file 到 output_base
  py$run_seqcre(input_base, output_base, subtype_file)

  invisible(NULL)
}
