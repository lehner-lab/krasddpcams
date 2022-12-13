#' normalise overlapping variants between blocks
#' 
#' This function allows you to normalise overlapping variants between blocks.
#' @param input dimsum data
#' 
#' @return normalized single mutation data
#' @export
Nor_overlap_single_mt_fitness<-function(input=input){
  input_single<-input[Nham_aa==1,]
  input_single[,nor_fitness_nooverlap:=sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("nor_fitness","nor_fitness_sigma"),by="aa_seq"]
  input_single[,nor_fitness_nooverlap_sigma:=sqrt(1/sum(1/.SD[[1]]^2, na.rm = T)),.SDcols = c("nor_fitness_sigma"),by="aa_seq"]
  input_single[,nor_gr_nooverlap:=sum(.SD[[1]]/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("nor_gr","nor_gr_sigma"),by="aa_seq"]
  input_single[,nor_gr_nooverlap_sigma:=sqrt(1/sum(1/.SD[[1]]^2, na.rm = T)),.SDcols = c("nor_gr_sigma"),by="aa_seq"]
  input_single<-input_single[!duplicated(aa_seq),]
  input_single
}