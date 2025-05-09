% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.NetConstruct.R
\name{run_network_construction}
\alias{run_network_construction}
\title{Protein-Protein Interaction (PPI) Network Construction Function}
\usage{
run_network_construction(
  diff_file,
  ppiScore = 400,
  subtype_file = "./subtype.xlsx",
  output_dir = "./Netconstruct_result/",
  ppi_method = c("String", "PhysicalPPIN", "ChengF")
)
}
\arguments{
\item{diff_file}{Character string specifying the path to the differential expression results Excel file.}

\item{ppiScore}{Numeric threshold for PPI interaction score (default: 400, used for String and ChengF methods).}

\item{subtype_file}{Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").}

\item{output_dir}{Character string specifying the output directory for PPI network files (default: "./Netconstruct_result/").}

\item{ppi_method}{Character vector specifying the PPI databases to use (options: "String", "PhysicalPPIN", "ChengF").}
}
\value{
Saves PPI network files for each subtype and method to the specified output directory.
}
\description{
This function constructs PPI networks for differentially expressed proteins across subtypes using specified PPI databases (String, PhysicalPPIN, ChengF).
}
\details{
DETAILS
}
\examples{
\dontrun{
setwd("your_workspace")
run_network_construction(
  diff_file = "./DEPS_result/diff_expression_results_all.xlsx",
  ppi_method = c("String", "PhysicalPPIN", "ChengF"),
  output_dir = "./Netconstruct_result/")
}
}
