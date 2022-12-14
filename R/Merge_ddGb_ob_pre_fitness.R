#' A Function to merge folding fitness against folding ddG
#' 
#' This function allows you to plot folding fitness against folding ddG.
#' @param prediction fitness data prediction
#' @param block1_dimsum_df block1_dimsum_df
#' @param block2_dimsum_df block2_dimsum_df
#' @param block3_dimsum_df block3_dimsum_df
#' @param assay_sele assay_sele
#' @param wt_aa_input wt_aa_input
#' 
#' @return Nothing
#' @export
Merge_ddGb_ob_pre_fitness<-function(prediction=prediction,
                                     block1_dimsum_df=block1_dimsum_df,
                                     block2_dimsum_df=block2_dimsum_df,
                                     block3_dimsum_df=block3_dimsum_df,
                                     assay_sele=assay_sele,
                                     wt_aa_input=wt_aa_input){
  pre<-fread(prediction)
  pre_pos<-Pos_id(input = pre,wt_aa = wt_aa_input)
  pre_pos
  load(block1_dimsum_df)
  block1<-as.data.table(all_variants)
  load(block2_dimsum_df)
  block2<-as.data.table(all_variants)
  load(block3_dimsum_df)
  block3<-as.data.table(all_variants)
  ### combine all data set and create final df
  data_before_nor<-rbind("block1"=block1,"block2"=block2,"block3"=block3,idcol="block",fill=T)
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
  ### linear transfor
  pre_nor<-pre_pos
  extract_prediction<-function(row){return(row[54+as.numeric(row[68])])}
  pre_nor$predicted_fitness<-apply(pre_nor,MARGIN = 1,FUN = extract_prediction)#create a new column for predicted fitness. In a matrix, 1 in MARGIN means rows
  pre_nor$predicted_fitness<-as.numeric(pre_nor$predicted_fitness)
  assay_sele_df<-data.table(assay=c("folding","DARPin K27","DARPin K55","PI3KCG","RAF1","RALGDS","SOS1"),
                            number=c(0,3,6,9,12,15,18))
  nb<-assay_sele_df[assay==assay_sele,number]
  pre_nor[phenotype==1+nb,pre_nor_mean_fitness:=mean]
  pre_nor[phenotype==2+nb,pre_nor_mean_fitness:=mean*d2+e2]
  pre_nor[phenotype==3+nb,pre_nor_mean_fitness:=mean*d3+e3]
  pre_nor[phenotype==1+nb,pre_nor_fitness_sigma:=std]
  pre_nor[phenotype==2+nb,pre_nor_fitness_sigma:=std*d2]
  pre_nor[phenotype==3+nb,pre_nor_fitness_sigma:=std*d3]
  pre_nor[phenotype==1+nb,ob_nor_fitness:=fitness]
  pre_nor[phenotype==2+nb,ob_nor_fitness:=fitness*d2+e2]
  pre_nor[phenotype==3+nb,ob_nor_fitness:=fitness*d3+e3]
  pre_nor[phenotype==1+nb,ob_nor_fitness_sigma:=sigma]
  pre_nor[phenotype==2+nb,ob_nor_fitness_sigma:=sigma*d2]
  pre_nor[phenotype==3+nb,ob_nor_fitness_sigma:=sigma*d3]
  pre_nor[phenotype==1+nb,pre_nor_fitness:=predicted_fitness]
  pre_nor[phenotype==2+nb,pre_nor_fitness:=predicted_fitness*d2+e2]
  pre_nor[phenotype==3+nb,pre_nor_fitness:=predicted_fitness*d3+e3]
  pre_nor
  # ggplot2::ggplot(data=pre_nor[phenotype==1+nb|phenotype==2+nb|phenotype==3+nb,],aes(x=ob_nor_fitness,y=pre_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black") +
  #   geom_smooth(method="lm") +
  #   stat_cor(method = "pearson",aes(label=..rr.label..), label.x=-1)+
  #   theme_bw()
  
  
}