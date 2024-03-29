% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/krasddpcams__plot_cor_gr_tecan_binding_MoCHi_prediction_double.R
\name{krasddpcams__plot_cor_gr_tecan_binding_MoCHi_prediction_double}
\alias{krasddpcams__plot_cor_gr_tecan_binding_MoCHi_prediction_double}
\title{validation of double mutation effects of tecan test vs MoCHi prediction}
\usage{
krasddpcams__plot_cor_gr_tecan_binding_MoCHi_prediction_double(
  tecandata = tecandata,
  plasmidid = plasmidid,
  ddG1 = ddG1,
  ddG2 = ddG2,
  colour_scheme
)
}
\arguments{
\item{tecandata}{tecandata output}

\item{plasmidid}{plasmid id}

\item{ddG1}{single mutation free energy change from DiMSum output of AbundancePCA}

\item{ddG2}{single mutation free energy change from DiMSum output of BindingPCA}

\item{colour_scheme}{colour scheme list}
}
\value{
Nothing
}
\description{
This function allows you to plot correlation between tecan culture growthrate and MoCHi prediction for double mutations.
}
