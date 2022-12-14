#' A Function to plot2D folding fitness against folding ddG
#' 
#' This function allows you to plot folding fitness against folding ddG.
#' @param prediction fitness data prediction
#' @param folding free energy of folding
#' @param binding free energy of binding
#' @param binding_assay binding assay's name
#' @param block_sele block's name
#' @param wt_aa_input wt_aa_input
#' 
#' @return Nothing
#' @export
Merge_ddGf_fitness<-function(prediction=prediction,
                                     folding_ddG=folding_ddG,
                                     block1_dimsum_df=block1_dimsum_df,
                                     block2_dimsum_df=block2_dimsum_df,
                                     block3_dimsum_df=block3_dimsum_df,
                                     wt_aa_input=wt_aa_input){
  pre<-fread(prediction)
  folding_ddG<-fread(folding_ddG)
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
  
  extract_additive_trait0<-function(row){return(row[68+as.numeric(row[68])*2-1])}
  extract_additive_trait1<-function(row){return(row[68+as.numeric(row[68])*2])}
  
  pre_nor$additive_trait0<-apply(pre_nor,MARGIN = 1,FUN = extract_additive_trait0)#create a new column for predicted fitness. In a matrix, 1 in MARGIN means rows
  pre_nor$additive_trait0<-as.numeric(pre_nor$additive_trait0)
  pre_nor$additive_trait1<-apply(pre_nor,MARGIN = 1,FUN = extract_additive_trait1)#create a new column for predicted fitness. In a matrix, 1 in MARGIN means rows
  pre_nor$additive_trait1<-as.numeric(pre_nor$additive_trait1)
  pre_nor[,additive_trait:=additive_trait0+additive_trait1]
  pre_nor[phenotype==1,pre_nor_mean_fitness:=mean]
  pre_nor[phenotype==2,pre_nor_mean_fitness:=mean*d2+e2]
  pre_nor[phenotype==3,pre_nor_mean_fitness:=mean*d3+e3]
  pre_nor[phenotype==1,pre_nor_fitness_sigma:=std]
  pre_nor[phenotype==2,pre_nor_fitness_sigma:=std*d2]
  pre_nor[phenotype==3,pre_nor_fitness_sigma:=std*d3]
  pre_nor[phenotype==1,ob_nor_fitness:=fitness]
  pre_nor[phenotype==2,ob_nor_fitness:=fitness*d2+e2]
  pre_nor[phenotype==3,ob_nor_fitness:=fitness*d3+e3]
  pre_nor[phenotype==1,ob_nor_fitness_sigma:=sigma]
  pre_nor[phenotype==2,ob_nor_fitness_sigma:=sigma*d2]
  pre_nor[phenotype==3,ob_nor_fitness_sigma:=sigma*d3]
  pre_nor[phenotype==1,pre_nor_fitness:=predicted_fitness]
  pre_nor[phenotype==2,pre_nor_fitness:=predicted_fitness*d2+e2]
  pre_nor[phenotype==3,pre_nor_fitness:=predicted_fitness*d3+e3]
  pre_nor
  # output<-list()
  # output[["p1"]]<-ggplot2::ggplot(data=pre_nor[Abundance1==1|Abundance2==1|Abundance3==1,],aes(x=additive_trait,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black") +
  #   geom_line(aes(x=additive_trait,y=pre_nor_fitness),color="red",trans = "log")+
  #   theme_bw()
  # output[["p2"]]<-ggplot2::ggplot(data=pre_nor[Abundance1==1,],aes(x=additive_trait,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black") +
  #   geom_line(aes(x=additive_trait,y=pre_nor_fitness),color="red",trans = "log")+
  #   theme_bw()
  # output[["p3"]]<-ggplot2::ggplot(data=pre_nor[Abundance2==1,],aes(x=additive_trait,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black") +
  #   geom_line(aes(x=additive_trait,y=pre_nor_fitness),color="red",trans = "log")+
  #   theme_bw()
  # output[["p4"]]<-ggplot2::ggplot(data=pre_nor[Abundance3==1,],aes(x=additive_trait,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black") +
  #   geom_line(aes(x=additive_trait,y=pre_nor_fitness),color="red",trans = "log")+
  #   theme_bw()
  # output[["p5"]]<-ggplot2::ggplot(data=pre_nor[Abundance1==1&is.na(AA_Pos2),],aes(x=additive_trait0,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black",trans = "log") +
  #   geom_line(aes(x=additive_trait0,y=mean),color="red")+
  #   theme_bw()
  # output[["p6"]]<-ggplot2::ggplot(data=pre_nor[Abundance2==1&is.na(AA_Pos2),],aes(x=additive_trait0,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black") +
  #   geom_line(aes(x=additive_trait0,y=mean),color="red")+
  #   theme_bw()
  # output[["p7"]]<-ggplot2::ggplot(data=pre_nor[Abundance3==1&is.na(AA_Pos2),],aes(x=additive_trait0,y=ob_nor_fitness))+
  #   geom_bin2d(bins = 100) +
  #   scale_fill_gradient(low="white",high="black",trans = "log") +
  #   geom_line(aes(x=additive_trait0,y=mean),color="red")+
  #   theme_bw()
  # output
  
  
}