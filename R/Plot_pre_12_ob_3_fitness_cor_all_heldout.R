#' A Function to plot predicted fitness from replicates 1&2 against observed fitness replicate 3 for all blocks
#' 
#' This function allows you to plot predicted fitness from replicates 1&2 against observed fitness replicate 3.
#' @param prediction prediction
#' @param observed observed 
#' @param assay_sele assay_sele
#' @param r_y r_y
#' @return Nothing
#' @export
Plot_pre_12_ob_3_fitness_cor_all_heldout<-function(prediction=prediction,
                                           observed=observed,
                                           assay_sele=assay_sele,
                                           r_y=r_y){
  pred_RAF_orig <- fread(observed)
  pred_dt_norep3 <- fread(prediction)
  pred_dt <- merge(
    pred_dt_norep3[,.SD,,.SDcols = names(pred_dt_norep3)[!names(pred_dt_norep3) %in% c("fitness", "sigma")]],
    pred_dt_orig[, .(aa_seq, phenotype, fitness = fitness3_uncorr, sigma = sigma3_uncorr)],
    by = c("aa_seq", "phenotype"))[!is.na(fitness)]
  for(i in 1:10){pred_dt[Fold==i, fitness_pred := .SD[[1]],,.SDcols = paste0('fold_', i)]}
  assay_sele_df<-data.table(assay=c("folding","DARPin K27","DARPin K55","PIK3CG","RAF1","RALGDS","SOS1"),
                            number=c(0,3,6,9,12,15,18))
  nb<-assay_sele_df[assay==assay_sele,number]
  plot_dt <- pred_dt[!is.na(fitness_pred)&(phenotype==1+nb|phenotype==2+nb|phenotype==3+nb),.(observed_fitness = fitness, predicted_fitness = fitness_pred)]
  lm_mochi<-lm(predicted_fitness~observed_fitness,plot_dt)
  ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "black") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10",guide = guide_colorbar(barwidth = 0.5,barheight = 1.5)) +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness (Rep. 3)") +
    ggplot2::ylab("Predicted fitness (Reps. 1&2)") +
    annotate("text",x=-1,y=r_y,
             label = paste0("R\u00B2 = ",round(summary(lm_mochi)$r.squared,2)),
             size=7*0.35 )+
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()+
    ggplot2::geom_abline(linetype="dashed")+
    ggplot2::ggtitle(assay_sele)+
    ggplot2::theme(text = element_text(size=7),
                   axis.text = element_text(size=7),
                   legend.text = element_text(size=7),
                   legend.key.size = unit(1, "cm"),
                   plot.title = element_text(size=7))+
    ggplot2::coord_fixed()
}