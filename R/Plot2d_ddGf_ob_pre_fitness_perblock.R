#' A Function to plot2D folding predicted fitness against observed fitness
#' 
#' This function allows you to plot folding predicted fitness against observed fitness.
#' @param pre_nor pre_nor
#' @param phenotypen phenotypen
#' 
#' @return Nothing
#' @export
Plot2d_ddGf_ob_pre_fitness_perblock<-function(pre_nor=pre_nor,
                                          phenotypen=phenotypen){
  pre_nor
  lm_mochi<-lm(pre_nor_fitness~ob_nor_fitness,pre_nor[phenotype==phenotypen,])
  ggplot()+
    stat_binhex(data=pre_nor[phenotype==phenotypen,],aes(x=ob_nor_fitness,y=pre_nor_fitness),
                bins = 50,size=0,color="black") +
    scale_fill_gradient(low="white",high="black",trans="log10",guide = guide_colorbar(barwidth = 0.5,barheight = 1.5)) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    geom_abline(intercept = 0,slope=1,linetype="dashed")+
    annotate("text",x=-0.8,y=0.3,
             label = paste0("R\u00B2 = ",round(summary(lm_mochi)$r.squared,2)),
             size=7*0.35 )+
    theme_classic()+
    xlab("Observed fitness")+
    ylab("Predicted fitness")+
    theme(text = element_text(size=7),
          axis.text = element_text(size=7),
          legend.text = element_text(size=7),
          legend.key.size = unit(1, "cm"),
          plot.title = element_text(size=7))+
    coord_fixed()
}