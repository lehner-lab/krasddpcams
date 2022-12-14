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
  ggplot2::ggplot()+
    stat_binhex(data=pre_nor[phenotype==phenotypen,],aes(x=ob_nor_fitness,y=pre_nor_fitness),
                bins = 50,size=0,color="black") +
    scale_fill_gradient(low="white",high="black",trans="log10") +
    annotate("text",x=-0.8,y=0.3,
             label = paste0("r = ",round(sqrt(summary(lm_mochi)$r.squared),3)),
             size=7*0.35 )+
    geom_smooth(data=pre_nor[phenotype==phenotypen,],aes(x=ob_nor_fitness,y=pre_nor_fitness),method=lm, se=T,size=0.3,color=colour_scheme[["red"]]) +
    theme_classic()+
    xlab("Observed fitness")+
    ylab("Predicted fitness")+
    theme(text = element_text(size=7),
          axis.text = element_text(size=7),
          legend.text = element_text(size=7),
          plot.title = element_text(size=7))+
    coord_fixed()
}