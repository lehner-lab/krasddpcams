#' A Function to get background mutations for each blocks
#' 
#' This function allows you to get background mutations for each blocks
#' @param input dimsum data
#' @param threshold count number to identify background
#' 
#' @return list by order of fitness of single mutations
#' @export
Get_background_order_fitness<-function(input=input,
                                       threshold=threshold){
  count_params <- function(dt) {
    counts <- lapply(dt[, -1, with = FALSE], function(col) table(col))
    return(counts)
  }
  param_counts <- count_params(input)
  mt1<-data.table(param_counts$mt1)
  list<-mt1[N>threshold,col]
  input_background<-input[Nham_aa==1&mt1%in%list,]
  output<-input_background[order(input_background[Nham_aa==1&mt1%in%list,-nor_fitness]),mt1]
  output
}
 
 
 