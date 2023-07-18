#' A Function to plot predicted fitness from replicates 1&2 against observed fitness replicate 3 (held out data)
#' 
#' This function allows you to plot predicted fitness from replicates 1&2 against observed fitness replicate 3.
#' @param prediction prediction
#' @param observed observed 
#' @param phenotypen phenotypen
#' 
#' @return Nothing
#' @export
Plot_pre_12_ob_3_fitness_heldout_cor<-function(prediction=prediction,
                                       observed=observed,
                                       phenotypen=phenotypen){
  #Load prediction data
  pred_dt_orig <- fread(observed)
  pred_dt_norep3 <- fread(prediction)
  #Merge
  pred_dt <- merge(
    pred_dt_norep3[,.SD,,.SDcols = names(pred_dt_norep3)[!names(pred_dt_norep3) %in% c("fitness", "sigma")]],
    pred_dt_orig[, .(aa_seq, phenotype, fitness = fitness3_uncorr, sigma = sigma3_uncorr)],
    by = c("aa_seq", "phenotype"))[!is.na(fitness)]
  #Predicted phenotype for held out test data
  for(i in 1:10){pred_dt[Fold==i, fitness_pred := .SD[[1]],,.SDcols = paste0('fold_', i)]}
  # #Plot performance on held out replicate (blocks separately)
  # for(i in names(pred_dt)[grepl("Abundance|Binding", names(pred_dt))]){
  #   plot_model_performance_heldoutrep(
  #     input_dt = pred_dt[get(i)==1],
  #     report_outpath = report_outpath,
  #     name = i)
  # }
  # #Plot performance on held out replicate (blocks combined)
  # for(i in c("Abundance", sapply(lapply(lapply(strsplit(names(pred_dt)[grepl("Binding", names(pred_dt))], "_"), unlist), rev), '[', 1))){
  #   pred_dt[, all_blocks := apply(.SD, 1, sum),,.SDcols = grepl(i, names(pred_dt))]
  #   name <- i
  #   if(i!="Abundance"){
  #     name <- paste0("Binding_", name)
  #   }
  #   plot_model_performance_heldoutrep(
  #     input_dt = pred_dt[all_blocks!=0],
  #     report_outpath = report_outpath,
  #     name = paste0(name, "_all"))
  # }
  # ####
  plot_dt <- pred_dt[!is.na(fitness_pred)&phenotype==phenotypen,.(observed_fitness = fitness, predicted_fitness = fitness_pred)]
  lm_mochi<-lm(predicted_fitness~observed_fitness,plot_dt)
  ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "black") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10",guide = guide_colorbar(barwidth = 0.5,barheight = 1.5)) +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness (Rep. 3)") +
    ggplot2::ylab("Predicted fitness (Reps. 1&2)") +
    annotate("text",x=-1,y=0.2,
             label = paste0("R\u00B2 = ",round(summary(lm_mochi)$r.squared,2)),
             size=7*0.35 )+
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()+
    ggplot2::geom_abline(linetype="dashed")+
    ggplot2::theme(text = element_text(size=7),
                   axis.text = element_text(size=7),
                   legend.text = element_text(size=7),
                   legend.key.size = unit(1, "cm"),
                   plot.title = element_text(size=7))+
    ggplot2::coord_fixed()
}