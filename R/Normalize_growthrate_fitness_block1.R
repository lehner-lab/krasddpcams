#' read fitness input for block1 only
#' 
#' This function allows you to input fitness for block1.
#' @param block1_dimsum_df dimsum data
#' 
#' @return Nothing
#' @export
Normalize_growthrate_fitness_block1<-function(block1_dimsum_df = block1_dimsum_df){
  ### load data
  load(block1_dimsum_df)
  block1<-as.data.table(all_variants)
  ### combine all data set and create final df
  data_before_nor<-block1
  data_before_nor[,block:="block1"]
  
  data_after_nor<-data_before_nor
  ### linear transform
  data_after_nor[block=="block1",nor_gr:=growthrate]
  data_after_nor[block=="block1",nor_gr_sigma:=growthrate_sigma]
  data_after_nor[block=="block1",nor_fitness:=fitness]
  data_after_nor[block=="block1",nor_fitness_sigma:=sigma]
  return(data_after_nor)
}