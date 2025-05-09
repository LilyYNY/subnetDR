#' @title Run SubnetDR Pipeline
#' @description Executes the complete SubnetDR pipeline for differential expression analysis, network construction, module division, selection, functional annotation, enrichment plotting, drug response analysis, and drug-target interaction processing.
#' @param base_path Character string specifying the base directory path containing input files and where output files will be saved.
#' @param py_env Character string specifying the path to the Python environment used for certain pipeline steps (e.g., drug response prediction).
#' @return The function generates various output files in the specified directories and prints a completion message.
#' @details details
#' @examples
#' \dontrun{
#' run_subnetDR_pipeline(
#'   base_path = "F:/subnetDR/sample_test",
#'   py_env = "C:/Users/YangMiao/.conda/envs/py3.9")
#' }
#' @rdname run_subnetDR_pipeline
#' @export
run_subnetDR_pipeline <- function(base_path , py_env) {
  setwd(base_path)

  # Differential expression analysis
  run_diff_expr_analysis(base_path)

  # Network construction
  run_network_construction(
    diff_file = "./DEPS_result/diff_expression_results_all.xlsx",
    ppi_method = c("String", "PhysicalPPIN", "ChengF"),
    output_dir = "./Netconstruct_result/"
  )

  # Module division
  module_division(
    network_method = c("String", "physicalPPIN", "chengF"),
    module_method = c("Louvain", "WF")
  )

  # Module selection
  module_selection(
    subtype_file = "./subtype.xlsx",
    network_method = c("String", "physicalPPIN", "chengF"),
    module_method = c("Louvain", "WF"),
    numberCutoff = 9
  )

  # Functional annotation analysis
  functional_annotation(
    subtype_file = "./subtype.xlsx",
    network_method = c("String", "physicalPPIN", "chengF"),
    module_method = c("Louvain", "WF")
  )

  # enrichment analysis
  plot_enrichment(
    subtype_file = "./subtype.xlsx",
    network_method = c("String", "physicalPPIN", "chengF"),
    module_method = c("Louvain", "WF")
  )

  # Drug response analysis
  drug_response_analysis(
    base_path = base_path,
    subtype_file = "subtype.xlsx",
    expr_base_path = "DrugResponse/subtype_expression",
    expr_data_file = "expression.xlsx",
    diff_expr_file = "DEPS_result/diff_expression_results_all.xlsx",
    module_base_path = "ModuleSelection",
    drug_response_path = "DrugResponse",
    gdsc_expr_file = "GDSC2_Expr.RData",
    gdsc_res_file = "GDSC2_Res.RData",
    network_methods = c("String", "physicalPPIN", "chengF"),
    module_methods = c("Louvain", "WF"),
    gene_col_name = "Gene"
  )

  run_SEQCre(
    input_base = ile.path(base_path, "DrugResponse"),
    output_base = file.path(base_path, "PRS/SMILES_SEQ"),
    subtype_file = file.path(base_path, "subtype.xlsx"),
    py_env = py_env
  )

  # Binding score analysis
  predict_BA(
    subtype_file = file.path(base_path, "subtype.xlsx"),
    drn_base = file.path(base_path, "DrugResponse"),
    prs_base = file.path(base_path, "PRS/SMILES_SEQ"),
    output_base = file.path(base_path, "PRS/PS"),
    py_env = py_env,
    network_methods = c("string", "physicalppin", "chengf"),
    module_methods = c("Louvain", "WF")
  )

  # Getting edge module files
  files <- list.files(
    path = file.path(base_path, "ModuleDivision"),
    pattern = "^edge_Module_.*\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  purrr::walk(files, split_edge_module, out_base = file.path(base_path, "PRS/PS"))

  # PRS and DTI analysis
  process_prs_dti(
    subtype_file = file.path(base_path, "subtype.xlsx"),
    prs_ps_base = file.path(base_path, "PRS/PS"),
    network_methods = c("string", "physicalppin", "chengf"),
    module_methods = c("Louvain", "WF"),
    py_env = py_env
  )

  #  ROC analysis
  dpi_ps_roc(
    dpi_base = file.path(base_path, "PRS/PS"),
    ttd_path = file.path(base_path, "PRS/PS", "DTI_TTD.RData"),
    drugbank_path = file.path(base_path, "PRS/PS", "DTI_DrugBank.RData")
  )

  message("SubnetDR pipeline completed successfully!")
}

#or Step by step
# detach(package: subnetDR)
# install.packages("subnetDR")
# library("subnetDR")
# setwd("F:/subnetDR/sample_test")
# #auto-pipeline
# run_subnetDR_pipeline(base_path = "F:/subnetDR/sample_test", py_env = "C:/Users/YangMiao/.conda/envs/py3.9")
#
# #step by step
# run_diff_expr_analysis("F:/subnetDR/sample_test")
#
# run_network_construction(
#   diff_file = "./DEPS_result/diff_expression_results_all.xlsx",
#   ppi_method = c("String", "PhysicalPPIN", "ChengF"),
#   output_dir = "./Netconstruct_result/"
# )
#
# module_division(network_method = c("String", "physicalPPIN", "chengF"),
#                 module_method = c("Louvain", "WF"))
# module_selection(subtype_file = "./subtype.xlsx",
#                  network_method = c("String", "physicalPPIN", "chengF"),
#                  module_method = c("Louvain", "WF"),
#                  numberCutoff = 9)
#
# functional_annotation(subtype_file = "./subtype.xlsx",
#                       network_method = c("String", "physicalPPIN", "chengF"),
#                       module_method = c("Louvain", "WF"))
# plot_enrichment(subtype_file = "./subtype.xlsx",
#                 network_method = c("String", "physicalPPIN", "chengF"),
#                 module_method = c("Louvain", "WF"))
#
# drug_response_analysis(
#   base_path = "F:/subnetDR/sample_test",
#   subtype_file = "subtype.xlsx",
#   expr_base_path = "DrugResponse/subtype_expression",
#   expr_data_file = "expression.xlsx",
#   diff_expr_file = "DEPS_result/diff_expression_results_all.xlsx",
#   module_base_path = "ModuleSelection",
#   drug_response_path = "DrugResponse",
#   gdsc_expr_file = "GDSC2_Expr.RData",
#   gdsc_res_file = "GDSC2_Res.RData",
#   network_methods = c("String", "physicalPPIN", "chengF"),
#   module_methods = c("Louvain", "WF"),
#   gene_col_name = "Gene"
# )
#
#
# run_SEQCre_R(
#   input_base   = "F:/subnetDR/sample_test/DrugResponse",
#   output_base  = "F:/subnetDR/sample_test/PRS/SMILES_SEQ",
#   subtype_file = "F:/subnetDR/sample_test/subtype.xlsx",
#   py_env       = "C:/Users/YangMiao/.conda/envs/py3.9"
# )
#
# predict_drug_response(
#   subtype_file = "F:/subnetDR/sample_test/subtype.xlsx",
#   drn_base = "F:/subnetDR/sample_test//DrugResponse",
#   prs_base = "F:/subnetDR/sample_test/PRS/SMILES_SEQ",
#   output_base = "F:/subnetDR/sample_teste/PRS/PS",
#   py_env = "C:/Users/YangMiao/.conda/envs/py3.9",
#   network_methods = c("string", "physicalppin", "chengf"),
#   module_methods = c("Louvain", "WF")
# )
#
# files <- list.files(
#   path       = "F:/subnetDR/sample_test/ModuleDivision",
#   pattern    = "^edge_Module_.*\\.txt$",
#   recursive  = TRUE,
#   full.names = TRUE
# )
# purrr::walk(files, split_edge_module)
#
# process_prs_dti(
#   subtype_file    = "F:/subnetDR/sample_test/subtype.xlsx",
#   prs_ps_base     = "F:/subnetDR/sample_test/PRS/PS",
#   network_methods = c("string","physicalppin","chengf"),
#   module_methods  = c("Louvain","WF"),
#   py_env          = "C:/Users/YangMiao/.conda/envs/py3.9"
# )
#
# dpi_ps_roc(
#   dpi_base      = "F:/subnetDR/sample_test/PRS/PS",
#   ttd_path      = file.path("F:/subnetDR/sample_test/PRS/PS", "DTI_TTD.RData"),
#   drugbank_path = file.path("F:/subnetDR/sample_test/PRS/PS", "DTI_DrugBank.RData")
# )
#
