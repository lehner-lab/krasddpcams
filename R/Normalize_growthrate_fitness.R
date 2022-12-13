#' A Function to normalize growth rates from dimsum output
#' 
#' This function allows you to normalize growth rates of different blocks from dimsum output according to the weighted mean stop and wt growth rates.
#' @param block1_dimsum_df block1's dimsum output
#' @param block2_dimsum_df block2's dimsum output
#' @param block3_dimsum_df block3's dimsum output
#' 
#' @return normalized fitness
#' @export
Normalize_growthrate_fitness<-function(block1_dimsum_df = block1_dimsum_df,block2_dimsum_df = block2_dimsum_df,block3_dimsum_df = block3_dimsum_df){
  ### load data
  load(block1_dimsum_df)
  block1<-as.data.table(all_variants)
  load(block2_dimsum_df)
  block2<-as.data.table(all_variants)
  load(block3_dimsum_df)
  block3<-as.data.table(all_variants)
  ### combine all data set and create final df
  data_before_nor<-rbind("block1"=block1,"block2"=block2,"block3"=block3,idcol="block",fill=T)
  
  ### normalized by error estimate (error estimated WT + error estimated dead variants)
  data_before_nor$gr_over_sigmasquared<-data_before_nor$growthrate/(data_before_nor$growthrate_sigma)**2
  data_before_nor$one_over_sigmasquared<-1/(data_before_nor$growthrate_sigma)**2
  dead_gr_block1<-data_before_nor[STOP==T&block=="block1",]
  dead_gr_block2<-data_before_nor[STOP==T&block=="block2",]
  dead_gr_block3<-data_before_nor[STOP==T&block=="block3",]
  stop1<-sum(dead_gr_block1$gr_over_sigmasquared, na.rm = TRUE)/sum(dead_gr_block1$one_over_sigmasquared, na.rm = TRUE)
  stop2<-sum(dead_gr_block2$gr_over_sigmasquared, na.rm = TRUE)/sum(dead_gr_block2$one_over_sigmasquared, na.rm = TRUE)
  stop3<-sum(dead_gr_block3$gr_over_sigmasquared, na.rm = TRUE)/sum(dead_gr_block3$one_over_sigmasquared, na.rm = TRUE)
  wt_gr_block1<-data_before_nor[WT==T&block=="block1",]
  wt_gr_block2<-data_before_nor[WT==T&block=="block2",]
  wt_gr_block3<-data_before_nor[WT==T&block=="block3",]
  wt1<-sum(wt_gr_block1$gr_over_sigmasquared, na.rm = TRUE)/sum(wt_gr_block1$one_over_sigmasquared, na.rm = TRUE)
  wt2<-sum(wt_gr_block2$gr_over_sigmasquared, na.rm = TRUE)/sum(wt_gr_block2$one_over_sigmasquared, na.rm = TRUE)
  wt3<-sum(wt_gr_block3$gr_over_sigmasquared, na.rm = TRUE)/sum(wt_gr_block3$one_over_sigmasquared, na.rm = TRUE)
  ### get coeffienct for linear transformation
  scaling_data<-data.frame(cbind(c(stop1,wt1),c(stop2,wt2),c(stop3,wt3)))
  colnames(scaling_data)<-c("block1","block2","block3")
  b1b2<-lm(formula = block1~block2,data = scaling_data)
  b1b2<-summary(b1b2)
  a2<-b1b2$coefficients[[2]]
  b2<-b1b2$coefficients[[1]]
  b1b3<-lm(formula = block1~block3,data = scaling_data)
  a3<-b1b3$coefficients[[2]]
  b3<-b1b3$coefficients[[1]]
  ### normalized by error estimate (error estimated WT + error estimated dead variants)
  data_before_nor$fitness_over_sigmasquared<-data_before_nor$fitness/(data_before_nor$sigma)**2
  data_before_nor$one_over_fitness_sigmasquared<-1/(data_before_nor$sigma)**2
  dead_fitness_block1<-data_before_nor[STOP==T&block=="block1",]
  dead_fitness_block2<-data_before_nor[STOP==T&block=="block2",]
  dead_fitness_block3<-data_before_nor[STOP==T&block=="block3",]
  stop1_fitness<-sum(dead_fitness_block1$fitness_over_sigmasquared, na.rm = TRUE)/sum(dead_fitness_block1$one_over_fitness_sigmasquared, na.rm = TRUE)
  stop2_fitness<-sum(dead_fitness_block2$fitness_over_sigmasquared, na.rm = TRUE)/sum(dead_fitness_block2$one_over_fitness_sigmasquared, na.rm = TRUE)
  stop3_fitness<-sum(dead_fitness_block3$fitness_over_sigmasquared, na.rm = TRUE)/sum(dead_fitness_block3$one_over_fitness_sigmasquared, na.rm = TRUE)
  wt_fitness_block1<-data_before_nor[WT==T&block=="block1",]
  wt_fitness_block2<-data_before_nor[WT==T&block=="block2",]
  wt_fitness_block3<-data_before_nor[WT==T&block=="block3",]
  wt1_fitness<-sum(wt_fitness_block1$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness_block1$one_over_fitness_sigmasquared, na.rm = TRUE)
  wt2_fitness<-sum(wt_fitness_block2$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness_block2$one_over_fitness_sigmasquared, na.rm = TRUE)
  wt3_fitness<-sum(wt_fitness_block3$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness_block3$one_over_fitness_sigmasquared, na.rm = TRUE)
  ### get coeffienct for linear transfor
  scaling_data_fitness<-data.frame(cbind(c(stop1_fitness,wt1_fitness),c(stop2_fitness,wt2_fitness),c(stop3_fitness,wt3_fitness)))
  colnames(scaling_data_fitness)<-c("block1","block2","block3")
  c1c2<-lm(formula = block1~block2,data = scaling_data_fitness)
  d2<-c1c2$coefficients[[2]]
  e2<-c1c2$coefficients[[1]]
  c1c3<-lm(formula = block1~block3,data = scaling_data_fitness)
  c1c3<-summary(c1c3)
  d3<-c1c3$coefficients[[2]]
  e3<-c1c3$coefficients[[1]]
  data_after_nor<-data_before_nor
  ### linear transfor
  data_after_nor[block=="block1",nor_gr:=growthrate]
  data_after_nor[block=="block2",nor_gr:=growthrate*a2+b2]
  data_after_nor[block=="block3",nor_gr:=growthrate*a3+b3]
  data_after_nor[block=="block1",nor_gr_sigma:=growthrate_sigma]
  data_after_nor[block=="block2",nor_gr_sigma:=growthrate_sigma*a2]
  data_after_nor[block=="block3",nor_gr_sigma:=growthrate_sigma*a3]
  data_after_nor[block=="block1",nor_fitness:=fitness]
  data_after_nor[block=="block2",nor_fitness:=fitness*d2+e2]
  data_after_nor[block=="block3",nor_fitness:=fitness*d3+e3]
  data_after_nor[block=="block1",nor_fitness_sigma:=sigma]
  data_after_nor[block=="block2",nor_fitness_sigma:=sigma*d2]
  data_after_nor[block=="block3",nor_fitness_sigma:=sigma*d3]
  return(data_after_nor)
}