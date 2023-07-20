
#' Plot fitness correlation between full RAF1 and RAF1RBD bindingPCA
#' 
#' This function allows you to plot fitness correlation between full RAF1 and RAF1RBD bindingPCA.
#' @param GAP_fitness dimsum output
#' @param RAF_fitness dimsum output
#' @param driven_mt1 driven mutation 1
#' @param driven_mt2 driven mutation 2
#' @param driven_mt3 driven mutation 3
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_density_GAP_RAF_single_gr_6mut<-function(
  GAP_fitness = GAP_fitness,
  RAF_fitness = RAF_fitness,
  driven_mt1 = driven_mt1,
  driven_mt2 = driven_mt2,
  driven_mt3 = driven_mt3,
  driven_mt4 = driven_mt4,
  driven_mt5 = driven_mt5,
  driven_mt6 = driven_mt6,
  colour_scheme
  ){
  ggplot2::ggplot()+
    ggplot2::geom_density(data=GAP_fitness[Nham_aa==1,],ggplot2::aes(x=nor_gr),color=colour_scheme[["red"]])+
    ggplot2::geom_vline(xintercept = GAP_fitness[mt1==driven_mt1&Nham_aa==1,nor_gr],color=colour_scheme[["red"]],linetype="dashed")+
    ggplot2::annotate("text",x=GAP_fitness[mt1==driven_mt1&Nham_aa==1,nor_gr],y=12,
             label = driven_mt1,
             size=7*0.35,color=colour_scheme[["red"]], angle=45)+
    ggplot2::geom_vline(xintercept = GAP_fitness[mt1==driven_mt2&Nham_aa==1,nor_gr],color=colour_scheme[["red"]],linetype="dashed")+
    ggplot2::annotate("text",x=GAP_fitness[mt1==driven_mt2&Nham_aa==1,nor_gr],y=12,
             label = driven_mt2,
             size=7*0.35 ,color=colour_scheme[["red"]], angle=45)+
    ggplot2::geom_vline(xintercept = GAP_fitness[mt1==driven_mt3&Nham_aa==1,nor_gr],color=colour_scheme[["red"]],linetype="dashed")+
    ggplot2::annotate("text",x=GAP_fitness[mt1==driven_mt3&Nham_aa==1,nor_gr],y=12,
             label = driven_mt3,
             size=7*0.35,color=colour_scheme[["red"]], angle=45)+
    ggplot2::geom_vline(xintercept = GAP_fitness[WT==TRUE,nor_gr],color=colour_scheme[["red"]])+
    ggplot2::geom_vline(xintercept = GAP_fitness[mt1==driven_mt4&Nham_aa==1,nor_gr],color=colour_scheme[["red"]],linetype="dashed")+
    ggplot2::annotate("text",x=GAP_fitness[mt1==driven_mt4&Nham_aa==1,nor_gr],y=12,
             label = driven_mt4,
             size=7*0.35,color=colour_scheme[["red"]], angle=45)+
    ggplot2::geom_vline(xintercept = GAP_fitness[mt1==driven_mt5&Nham_aa==1,nor_gr],color=colour_scheme[["red"]],linetype="dashed")+
    ggplot2::annotate("text",x=GAP_fitness[mt1==driven_mt5&Nham_aa==1,nor_gr],y=12,
             label = driven_mt5,
             size=7*0.35 ,color=colour_scheme[["red"]], angle=45)+
    ggplot2::geom_vline(xintercept = GAP_fitness[mt1==driven_mt6&Nham_aa==1,nor_gr],color=colour_scheme[["red"]],linetype="dashed")+
    ggplot2::annotate("text",x=GAP_fitness[mt1==driven_mt6&Nham_aa==1,nor_gr],y=12,
             label = driven_mt6,
             size=7*0.35,color=colour_scheme[["red"]], angle=45)+
    ggplot2::geom_vline(xintercept = GAP_fitness[WT==TRUE,nor_gr],color=colour_scheme[["red"]])+
    ggplot2::geom_density(data=RAF_fitness[block=="block1"&Nham_aa==1,],ggplot2::aes(x=nor_gr),color=colour_scheme[["blue"]])+
    ggplot2::geom_vline(xintercept = RAF_fitness[mt1==driven_mt1&Nham_aa==1,nor_gr],color=colour_scheme[["blue"]],linetype="dashed")+
    ggplot2::annotate("text",x=RAF_fitness[mt1==driven_mt1&Nham_aa==1,nor_gr],y=12,
             label = driven_mt1,
             size=7*0.35,color=colour_scheme[["blue"]], angle=45)+
    ggplot2::geom_vline(xintercept = RAF_fitness[mt1==driven_mt2&Nham_aa==1,nor_gr],color=colour_scheme[["blue"]],linetype="dashed")+
    ggplot2::annotate("text",x=RAF_fitness[mt1==driven_mt2&Nham_aa==1,nor_gr],y=12,
             label = driven_mt2,
             size=7*0.35,color=colour_scheme[["blue"]], angle=45)+
    ggplot2::geom_vline(xintercept = RAF_fitness[mt1==driven_mt3&Nham_aa==1,nor_gr],color=colour_scheme[["blue"]],linetype="dashed")+
    ggplot2::annotate("text",x=RAF_fitness[mt1==driven_mt3&Nham_aa==1,nor_gr],y=12,
             label = driven_mt3,
             size=7*0.35,color=colour_scheme[["blue"]], angle=45)+
    ggplot2::geom_vline(xintercept = RAF_fitness[mt1==driven_mt4&Nham_aa==1,nor_gr],color=colour_scheme[["blue"]],linetype="dashed")+
    ggplot2::annotate("text",x=RAF_fitness[mt1==driven_mt4&Nham_aa==1,nor_gr],y=12,
             label = driven_mt4,
             size=7*0.35,color=colour_scheme[["blue"]], angle=45)+
    ggplot2::geom_vline(xintercept = RAF_fitness[mt1==driven_mt5&Nham_aa==1,nor_gr],color=colour_scheme[["blue"]],linetype="dashed")+
    ggplot2::annotate("text",x=RAF_fitness[mt1==driven_mt5&Nham_aa==1,nor_gr],y=12,
             label = driven_mt5,
             size=7*0.35,color=colour_scheme[["blue"]], angle=45)+
    ggplot2::geom_vline(xintercept = RAF_fitness[mt1==driven_mt6&Nham_aa==1,nor_gr],color=colour_scheme[["blue"]],linetype="dashed")+
    ggplot2::annotate("text",x=RAF_fitness[mt1==driven_mt6&Nham_aa==1,nor_gr],y=12,
             label = driven_mt6,
             size=7*0.35,color=colour_scheme[["blue"]], angle=45)+
    ggplot2::geom_vline(xintercept = RAF_fitness[WT==TRUE,nor_gr],color=colour_scheme[["blue"]])+
    ggplot2::theme_classic()+
    ggplot2::xlab("growth rates in BindingPCA")+
    ggplot2::ylab("Density")+
    ggplot2::theme(text = ggplot2::element_text(size=7),
          axis.text = ggplot2::element_text(size=7),
          legend.text = ggplot2::element_text(size=7),
          legend.key.size = ggplot2::unit(1, "cm"),
          plot.title = ggplot2::element_text(size=7))
}