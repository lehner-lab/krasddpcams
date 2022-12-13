#' calculate weighted mean free energy
#' 
#' This function allows you to calculate weighted mean free energy change.
#' @param ddG free energy data
#' @return weighted mean free energy changes
#' @export 
Get_weighted_mean_ddG<-function(ddG=ddG,
                                assay_sele=assay_sele){
  ddG<-Read_ddG(ddG = ddG,
                assay_sele = assay_sele)
  output<-ddG[Pos_real>1,sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(output,"V1","mean")
  output_sigma<-ddG[Pos_real>1,sqrt(1/sum(1/.SD[[2]]^2, na.rm = T)),.SDcols = c("mean_kcal/mol","std_kcal/mol"),by="Pos_real"]
  setnames(output_sigma,"V1","sigma")
  weighte_mean<-merge(output,output_sigma,by="Pos_real")
  weighte_mean
}