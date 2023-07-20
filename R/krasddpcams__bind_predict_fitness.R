
#' prediction binding fitness from folding and binding energy change
#' 
#' This function allows you to plot fitness distribution.
#' @param mochi_parameters mochi_parameters 
#' @param folding_energy folding_energy 
#' @param binding_energy binding_energy
#' @param RT 0.001987*(273+30)
#' 
#' @return predicted binding fitness
#' @export
#' @import data.table
krasddpcams__bind_predict_fitness<-function(
  mochi_parameters=mochi_parameters,
  folding_energy=folding_energy,
  binding_energy=binding_energy,
  RT=0.001987*(273+30)
  ){
  mochi_model<-fread(mochi_parameters)
  fraction_bound<-krasddpcams__fraction_bound(folding_energy = folding_energy,binding_energy = binding_energy)
  fitness_binding<- fraction_bound * mochi_model[fold==1,kernel] + mochi_model[fold==1,bias]
  fitness_binding
}