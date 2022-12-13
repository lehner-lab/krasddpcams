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
Plot2d_ddGf_fitness<-function(pre_nor=pre_nor,
                                   mochi_parameters=mochi_parameters,
                                   fold_n=fold_n,
                                   phenotypen=phenotypen,RT=RT,
                                   bin_input=bin_input){
  pre_nor_input<-pre_nor[phenotype==phenotypen,]
  folding_range<-pre_nor_input[Nham_aa>0,range(fold_1_additive_trait0*RT,na.rm=T)]
  folding_energy_grid<-seq(folding_range[1],
                           folding_range[2],
                           (folding_range[2]-folding_range[1])/500)
  fitness_folding <- doubledeepms__abundance_predict_fitness(mochi_parameters=mochi_parameters,
                                                             fold_n=fold_n,
                                                             folding_energy = folding_energy_grid,
                                                             RT = RT)
  
  pred_fitness_dt <- data.table(
    f_dg_pred = rep(folding_energy_grid, 2),
    observed_fitness = fitness_folding,
    mut_order = rep(c(1, 2), each = length(folding_energy_grid)))
  
  folding_fraction<-doubledeepms__fraction_folded(folding_energy = folding_energy_grid)
  
  
  ggplot()+
    stat_binhex(data=pre_nor_input[Nham_aa>0,],aes(x=fold_1_additive_trait0*RT,y=fitness),bins=bin_input,size=0,color="black")+
    scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    geom_line(data = pred_fitness_dt,aes(x=f_dg_pred,y=observed_fitness), color = "red")+
    geom_vline(xintercept=0,size=0.1,color="black")+
    geom_hline(yintercept=0,size=0.1,color="black")+
    theme_classic2()+
    ylab("Abundance(observed)")+
    xlab("Folding ddG (inferred)")+
    theme(text = element_text(size=7),
          legend.position="right",
          legend.text = element_text(size=7),
          axis.text.x = element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,-0.5,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))
}