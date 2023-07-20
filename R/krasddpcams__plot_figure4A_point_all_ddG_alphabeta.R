
#' all ddG point plot
#' 
#' This function allows you to plot all ddG heatmap.
#' @param ddG1 free energy data1
#' @param assay1 data1's name
#' @param ddG2 free energy data2
#' @param assay2 data2's name
#' @param ddG3 free energy data3
#' @param assay3 data3's name
#' @param ddG4 free energy data4
#' @param assay4 data4's name
#' @param ddG5 free energy data5
#' @param assay5 data5's name
#' @param ddG6 free energy data6
#' @param assay6 data6's name
#' @param ddG7 free energy data7
#' @param assay7 data7's name
#' @param anno annotation df final_distance_dc_anno_for_anno_mt_4
#' @param rect_input rect for beta strands
#' @param rect_alpha rect for alpha 1
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_figure4A_point_all_ddG_alphabeta<-function(
  ddG1=ddG1,
  assay1=assay1,
  ddG2=ddG2,
  assay2=assay2,
  ddG3=ddG3,
  assay3=assay3,
  ddG4=ddG4,
  assay4=assay4,
  ddG5=ddG5,
  assay5=assay5,
  ddG6=ddG6,
  assay6=assay6,
  ddG7=ddG7,
  assay7=assay7,
  anno=anno,
  rect_input=rect_input,
  rect_alpha=rect_alpha,
  colour_scheme
  ){
  ddG1<-krasddpcams__read_ddG(ddG1,assay1)
  ddG2<-krasddpcams__read_ddG(ddG2,assay2)
  ddG3<-krasddpcams__read_ddG(ddG3,assay3)
  ddG4<-krasddpcams__read_ddG(ddG4,assay4)
  ddG5<-krasddpcams__read_ddG(ddG5,assay5)
  ddG6<-krasddpcams__read_ddG(ddG6,assay6)
  ddG7<-krasddpcams__read_ddG(ddG7,assay7)
  rects_dt<-as.data.table(rect_input)
  rects_alphas<-as.data.table(rect_alpha)
  all_ddG<-rbind(ddG1,ddG2,ddG3,ddG4,ddG5,ddG6,ddG7)
  all_data_25_anno_bs<-merge(all_ddG,
                             anno,
                             by.x=c("Pos_real","assay"),
                             by.y=c("Pos","assay"),all=T)
  all_data_25_anno_bs[is.na(binding),binding:="no"]
  all_data_25_anno_bs<-within(all_data_25_anno_bs, 
                              assay <- factor(assay,
                                              levels=c("folding",
                                                       "RAF1",
                                                       "PIK3CG",
                                                       "RALGDS",
                                                       "SOS1",
                                                       "DARPin K27",
                                                       "DARPin K55")))
  all_data_25_anno_bs<-within(all_data_25_anno_bs, 
                              binding <- factor(binding,
                                                levels=c("core",
                                                         "surface",
                                                         "HVR",
                                                         "GTP",
                                                         "GDP",
                                                         "both",
                                                         "RAF",
                                                         "PI3",
                                                         "RAL",
                                                         "SOS",
                                                         "K27",
                                                         "K55",
                                                         "no")))
  ggplot2::ggplot()+
    ggplot2::geom_rect(data=rects_dt,ggplot2::aes(ymin=-2,ymax=3.5,xmin=xstart-0.5,xmax=xend+0.5),fill="black",alpha=0.08)+
    ggplot2::geom_rect(data=rects_alphas,ggplot2::aes(ymin=-2,ymax=3.5,xmin=xstart-0.5,xmax=xend+0.5),fill="black",alpha=0.05)+
    ggplot2::geom_jitter(data=all_data_25_anno_bs[id!="WT"],ggplot2::aes(x=Pos_real,y=`mean_kcal/mol`,color=binding),size=0.1,width=0.4,height = 0)+
    ggplot2::facet_wrap(~assay,nrow = 7)+
    ggplot2::scale_x_continuous(expand=c(1/188,11/188),breaks=seq(0,187,5),
                       minor_breaks = seq(2,187,1))+
    ggplot2::scale_color_manual(values=c(colour_scheme[["dark green"]], colour_scheme[["orange"]],"gray",
                                colour_scheme[["blue"]],colour_scheme[["blue"]],
                                colour_scheme[["purple"]],colour_scheme[["red"]],colour_scheme[["red"]],
                                colour_scheme[["red"]],colour_scheme[["red"]],colour_scheme[["red"]],colour_scheme[["red"]],"gray"))+
    ggplot2::theme_classic()+
    ggplot2::theme(text = ggplot2::element_text(size = 7),
          axis.text = ggplot2::element_text(size = 7),legend.text = ggplot2::element_text(size = 7),
          legend.position="bottom",
          panel.spacing.y = ggplot2::unit(1, "mm"),
          strip.background = ggplot2::element_rect(color="black",fill=NULL,linetype="solid")
    )+
    ggplot2::coord_fixed(ratio=2)
}