
#' calculate weighted mean absolute free energy
#' 
#' This function allows you to calculate weighted mean absolute free energy change.
#' @param ddG free energy data
#' @return Nothing
#' @export 
#' @import data.table
krasddpcams__get_weighted_mean_abs_ddG<-function(
  ddG=ddG,
  assay_sele=assay_sele
  ){
  ddG<-krasddpcams__read_ddG(ddG = ddG,
                assay_sele = assay_sele)
  output<-ddG[Pos_real>1,sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(output,"V1","mean")
  output_sigma<-ddG[Pos_real>1,sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(output_sigma,"V1","sigma")
  weighte_mean<-merge(output,output_sigma,by="Pos_real")
  weighte_mean
}