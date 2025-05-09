% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.ModuleAnnotation.R
\name{plot_enrichment}
\alias{plot_enrichment}
\title{Plot Enrichment Results for All Modules}
\usage{
plot_enrichment(
  subtype_file = "./subtype.xlsx",
  base_input_path = "./FunctionalAnnotation/",
  network_method = c("String", "physicalPPIN", "chengF"),
  module_method = c("Louvain", "WF")
)
}
\arguments{
\item{subtype_file}{Character string specifying the path to the subtype phenotype data file (default: "./subtype.xlsx").}

\item{base_input_path}{Character string specifying the base directory containing functional annotation files (default: "./FunctionalAnnotation/").}

\item{network_method}{Character vector specifying the PPI network databases (options: "String", "physicalPPIN", "chengF").}

\item{module_method}{Character vector specifying the module detection methods (options: "Louvain", "WF").}
}
\value{
Saves dot plots as PNG files to the specified input directory.
}
\description{
Generates dot plots for enrichment results stored in Excel files for each module.
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
setwd("your_workspace")
Run functional annotation for all combinations
functional_annotation(subtype_file = "./subtype.xlsx",
                     network_method = c("String", "physicalPPIN", "chengF"),
                     module_method = c("Louvain", "WF"))
 }
}
}
\seealso{
\code{\link[scales]{pretty_breaks}}
}
