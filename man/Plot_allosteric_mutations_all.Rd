% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_allosteric_mutations_all.R
\name{Plot_allosteric_mutations_all}
\alias{Plot_allosteric_mutations_all}
\title{Plot allosteric mutations}
\usage{
Plot_allosteric_mutations_all(
  ddG = ddG,
  anno = anno,
  assay_sele = assay_sele,
  threshold = threshold,
  rect_input = rect_input,
  rect_alpha = rect_alpha,
  y_min = y_min,
  y_max = y_max
)
}
\arguments{
\item{ddG}{free energy data}

\item{anno}{structure information annotation data}

\item{assay_sele}{assay:"RAF","PI3"...}

\item{threshold}{threshold combined}

\item{rect_input}{rect_input}

\item{rect_alpha}{rect_alpha}

\item{y_min}{y_min}

\item{y_max}{y_max}
}
\value{
Nothing
}
\description{
This function allows you to plot allosteric mutations.
}
