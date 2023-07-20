
#' Plot ROC-AUC for predicting contact residues by ddG or fitness
#' 
#' This function allows you to plot ROC-AUC for predicting contact residues by ddG or fitness
#' @param fitness_all fitness_all
#' @param weighted_meab_abs_ddGf weighted_meab_abs_ddGf
#' @param weighted_meab_abs_ddGb weighted_meab_abs_ddGb
#' @param anno_input anno_input
#' @param bind bind
#' @param bind_ligand bind_ligand
#' 
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface<-function(
  fitness_all=fitness_all,
  weighted_meab_abs_ddGf=weighted_meab_abs_ddGf,
  weighted_meab_abs_ddGb=weighted_meab_abs_ddGb,
  anno_input=anno_input,
  bind=bind,
  bind_ligand=bind_ligand
  ){
  abundance_fitness<-fitness_all[Nham_aa==1&assay=="stab",.(AA_Pos1,nor_fitness)]
  abundance_fitness_mean<-abundance_fitness[,mean(.SD[[1]],na.rm = T),.SDcols=c("nor_fitness"),by="AA_Pos1"]
  colnames(abundance_fitness_mean)<-c("Pos_real","mean_abundance_fitness")
  binding_fitness<-fitness_all[Nham_aa==1&assay==bind,.(AA_Pos1,nor_fitness)]
  binding_fitness_mean<-binding_fitness[,mean(.SD[[1]],na.rm = T),.SDcols=c("nor_fitness"),by="AA_Pos1"]
  colnames(binding_fitness_mean)<-c("Pos_real","mean_binding_fitness")
  ddGf<-weighted_meab_abs_ddGf[,c(1:2)]
  colnames(ddGf)<-c("Pos_real","weighted_mean_ddGf")
  ddGb<-weighted_meab_abs_ddGb[,c(1:2)]
  colnames(ddGb)<-c("Pos_real","weighted_mean_ddGb")
  anno<-anno_input[,c(1:28)]
  anno[,boostDM:=0]
  anno[get(bind_ligand)<5,boostDM:=1]
  anno<-anno[,.(Pos,boostDM)]
  
  test_file_final_dt<-merge(abundance_fitness_mean,binding_fitness_mean,by="Pos_real")
  test_file_final_dt<-merge(test_file_final_dt,ddGf,by="Pos_real")
  test_file_final_dt<-merge(test_file_final_dt,ddGb,by="Pos_real")
  test_file_final_dt<-merge(test_file_final_dt,anno,by.x="Pos_real",by.y="Pos")
  
  metric_names_plot<-c("mean_binding_fitness","weighted_mean_ddGb")
  subset_dt <- test_file_final_dt[order(Pos_real)][!duplicated(Pos_real)]
  perf_list<-list()
  for(name in metric_names_plot){
    subset_dt[, plot_metric := .SD[[1]],,.SDcols = name]
    roc_df <- data.frame(
      predictions = subset_dt[,plot_metric], 
      labels = subset_dt[,boostDM])
    roc_df[,"predictions_lm"] <- unlist(lm(labels~predictions, data = roc_df)["fitted.values"])
    pred <- ROCR::prediction(roc_df$predictions_lm, roc_df$labels)
    perf <- ROCR::performance(pred,"tpr","fpr")
    auc <- round(ROCR::performance(pred, measure = "auc")@'y.values'[[1]], 2)
    #Save
    perf_list[[name]] <- data.table(
      FPR = perf@'x.values'[[1]],
      TPR = perf@'y.values'[[1]],
      measure = name,
      auc = auc)
  }
  perf_list
  # test_df<-as.data.frame(test_file_final_dt)
  # Y = test_df[,1]
  # split <- sample.split(Y, SplitRatio = 0.8)
  # train <- subset(test_df, split == "TRUE")
  # test <- subset(test_df, split == "FALSE")
  # #ddGf
  # model_ddGf = glm(boostDM~weighted_mean_ddGf,
  #                 train , family="binomial")
  # summary(model_ddGf)
  # pred_test_ddGf <- predict(model_ddGf,test,type="response")
  # test_prob_ddGf = predict(model_ddGf, test, type = "response")
  # plot(roc(test$boostDM ~ test_prob_ddGf), print.auc = TRUE, 
  #      col = "blue", print.auc.y = .4, add = F)
  # #ddGb
  # model_ddGb = glm(boostDM~weighted_mean_ddGb,
  #                  train , family="binomial")
  # summary(model_ddGb)
  # pred_test_ddGb <- predict(model_ddGb,test,type="response")
  # test_prob_ddGb = predict(model_ddGb, test, type = "response")
  # plot(roc(test$boostDM ~ test_prob_ddGb), print.auc = TRUE,
  #      col = "red", print.auc.y = .5, add = T)
  # #abundance fitness
  # model_fitness_abundance = glm(boostDM~mean_abundance_fitness,
  #                  train , family="binomial")
  # summary(model_fitness_abundance)
  # pred_test_fitness_abundance <- predict(model_fitness_abundance,test,type="response")
  # test_prob_fitness_abundance = predict(model_fitness_abundance, test, type = "response")
  # plot(roc(test$boostDM ~ test_prob_fitness_abundance), print.auc = TRUE,
  #      col = "green", print.auc.y = .3, add = T)
  # #binding fitness
  # model_fitness_binding = glm(boostDM~mean_binding_fitness,
  #                  train , family="binomial")
  # summary(model_fitness_binding)
  # pred_test_fitness_binding <- predict(model_fitness_binding,test,type="response")
  # test_prob_fitness_binding = predict(model_fitness_binding, test, type = "response")
  # plot(roc(test$boostDM ~ test_prob_fitness_binding), print.auc = TRUE,
  #      col = "yellow", print.auc.y = .2, add = T)
  # 
  
  plot_dt <- rbindlist(perf_list)
  # plot_dt
  plot_cols <- c("black","red")
  names(plot_cols) <- metric_names_plot
  plot_dt<-within(plot_dt,
                  measure<-factor(measure,
                                  levels=c("weighted_mean_ddGb","mean_binding_fitness")))
  
  #Plot
  auc_dt <- plot_dt[!duplicated(measure)][order(measure, decreasing = T)]
  auc_dt[, FPR := 0.5]
  auc_dt[, TPR := seq(0, 1, 1/(.N+1))[2:(.N+1)]]
  auc_dt<-within(auc_dt,
                 measure<-factor(measure,
                                 levels=c("weighted_mean_ddGb","mean_binding_fitness")))
  
  ggplot2::ggplot(plot_dt,ggplot2::aes(FPR, TPR, color = measure)) +
    ggplot2::geom_line(size=0.35) +
    ggplot2::geom_abline(linetype = 2,size=0.35) +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") +
    ggplot2::geom_text(data = auc_dt, ggplot2::aes(label=paste("AUC = ", auc, sep=""),color=measure),size=7*0.35) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values=plot_cols,name=NULL,labels=c("mean binding fitness","weighted mean of\nabsolute ddGb"))+
    ggplot2::theme(text = ggplot2::element_text(size = 7),
                   axis.text = ggplot2::element_text(size = 7),legend.text = ggplot2::element_text(size = 7),
                   legend.position="top",legend.direction = "vertical")+
    ggplot2::coord_fixed()
}