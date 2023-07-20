
#' prediction abundance fitness from folding energy change
#' 
#' This function allows you to plot fitness distribution.
#' @param input dimsum data
#' @param assay_name folding or binding
#' @param anno_input anno_input
#' @param RT R*T
#'  
#' @return predicted abundance fitness
#' @export
#' @import data.table
krasddpcams__abundance_predict_fitness<-function(
  mochi_parameters=mochi_parameters,
  fold_n=fold_n,
  folding_energy=folding_energy,
  RT=0.001987*(273+30)
  ){
  mochi_model<-fread(mochi_parameters)
  folding_fraction<-krasddpcams__fraction_folded(folding_energy = folding_energy,RT=RT)
  fitness_folding<- folding_fraction * mochi_model[fold==fold_n,kernel] + mochi_model[fold==fold_n,bias]
  fitness_folding
}