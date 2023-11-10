#' A Function to plot2D predicted fitness against observed fitness
#' 
#' This function allows you to plot predicted fitness against observed fitness.
#' @param prediction fitness data prediction
#' @param block1_dimsum_df block1_dimsum_df
#' @param block2_dimsum_df block2_dimsum_df
#' @param block3_dimsum_df block3_dimsum_df
#' @param assay_sele assay_sele
#' @param wt_aa_input wt_aa_input
#' 
#' @return model_performance
#' @export
#' @import data.table
krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test<-function(prediction=prediction,
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
  for (i in 1:10){
    x<-as.factor(i)
    eval(parse(text=paste0("pre_nor",x,"=data.table(pre_nor)")
    ))
    eval(parse(text=paste0("pre_nor",x,"[,predicted_fitness:=fold_",x,"]")
    ))
    eval(parse(text=paste0("pre_nor",x,"[,dataset:=",x,"]")
    ))
  }
  
  pre_nor_all<-rbind(pre_nor1,pre_nor2,pre_nor3,pre_nor4,pre_nor5,pre_nor6,pre_nor7,pre_nor8,pre_nor9,pre_nor10)
  pre_nor<-pre_nor_all
  
  
  assay_sele_df<-data.table(assay=c("folding","DARPin K27","DARPin K55","PIK3CG","RAF1","RALGDS","SOS1"),
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
  pre_nor_training<-pre_nor[Fold!=dataset&!is.na(Fold),]
  lm_mochi_training<-lm(predicted_fitness~fitness,pre_nor_training[phenotype==1+nb|phenotype==2+nb|phenotype==3+nb,])
  
  model_performance<-list()
  model_performance[["training"]]<-ggplot(data=pre_nor_training[phenotype==1+nb|phenotype==2+nb|phenotype==3+nb,],aes(x=fitness,y=predicted_fitness))+
    stat_binhex(bins = 50,size=0,color="black") +
    scale_fill_gradient(low="white",high="black",trans="log10",guide = guide_colorbar(barwidth = 0.5,barheight = 1.5)) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(intercept = 0,slope=1,linetype="dashed")+
    annotate("text",x=-1,y=0.5,
             label = paste0("R\u00B2 = ",round(summary(lm_mochi_training)$r.squared,2)),
             size=7*0.35 )+
    theme_classic()+
    xlab("Observed fitness")+
    ylab("Predicted fitness")+
    ggtitle(assay_sele)+
    theme(text = element_text(size=7),
          axis.text = element_text(size=7),
          legend.text = element_text(size=7),
          plot.title = element_text(size=7))+
    coord_fixed()
  pre_nor_test<-pre_nor[Fold==dataset&!is.na(Fold),]
  lm_mochi_test<-lm(predicted_fitness~fitness,pre_nor_test[phenotype==1+nb|phenotype==2+nb|phenotype==3+nb,])
  model_performance[["test"]]<-ggplot(data=pre_nor_test[phenotype==1+nb|phenotype==2+nb|phenotype==3+nb,],aes(x=fitness,y=predicted_fitness))+
    stat_binhex(bins = 50,size=0,color="black") +
    scale_fill_gradient(low="white",high="black",trans="log10",guide = guide_colorbar(barwidth = 0.5,barheight = 1.5)) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(intercept = 0,slope=1,linetype="dashed")+
    annotate("text",x=-1,y=0.5,
             label = paste0("R\u00B2 = ",round(summary(lm_mochi_test)$r.squared,2)),
             size=7*0.35 )+
    theme_classic()+
    xlab("Observed fitness")+
    ylab("Predicted fitness")+
    ggtitle(assay_sele)+
    theme(text = element_text(size=7),
          axis.text = element_text(size=7),
          legend.text = element_text(size=7),
          plot.title = element_text(size=7))+
    coord_fixed()
  model_performance[["R_squared"]]<-data.table(dataset=c("training","test"),R_squared=c(summary(lm_mochi_training)$r.squared,summary(lm_mochi_test)$r.squared))
  return(model_performance)
}
