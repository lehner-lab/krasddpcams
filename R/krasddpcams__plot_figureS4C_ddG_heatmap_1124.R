
#' all ddG heatmap
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
#' @param wt_aa wt amino acid sequence
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_figureS4C_ddG_heatmap_1124<-function(
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
  wt_aa,
  colour_scheme
  ){
  ddG1<-krasddpcams__read_ddG(ddG1,assay1)
  ddG2<-krasddpcams__read_ddG(ddG2,assay2)
  ddG3<-krasddpcams__read_ddG(ddG3,assay3)
  ddG4<-krasddpcams__read_ddG(ddG4,assay4)
  ddG5<-krasddpcams__read_ddG(ddG5,assay5)
  ddG6<-krasddpcams__read_ddG(ddG6,assay6)
  ddG7<-krasddpcams__read_ddG(ddG7,assay7)
  all_ddG<-rbind(ddG1,ddG2,ddG3,ddG4,ddG5,ddG6,ddG7)
  all_ddG[wt_codon==mt_codon,`mean_kcal/mol`:=0]
  all_ddG<-within(all_ddG,
                  assay <- factor(assay,
                                  levels=c(assay1,
                                           assay2,
                                           assay3,
                                           assay4,
                                           assay5,
                                           assay6,
                                           assay7)))
  all_ddG<-within(all_ddG,
                  mt_codon <- factor(mt_codon, 
                                     levels = c("D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  ggplot2::ggplot(all_ddG[Pos_real>1,],ggplot2::aes(y=mt_codon,x=Pos_real))+
    ggplot2::geom_tile(ggplot2::aes(fill=`mean_kcal/mol`))+
    ggplot2::scale_x_discrete(limits=c(2:188),labels=c(strsplit(wt_aa,"")))+
    ggplot2::scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value = "white")+
    ggplot2::geom_text(data=all_ddG[Pos_real>1&wt_codon==mt_codon,],
              ggplot2::aes(x=Pos_real,y=mt_codon,label="-"),size=5*0.35)+
    
    ggplot2::facet_wrap(~assay,ncol=1,nrow = 7)+
    ggplot2::ggtitle("free energy change")+
    ggpubr::theme_classic2()+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(family="Courier",size = 5, vjust = 3,hjust = 0.5),
          axis.text.y = ggplot2::element_text(family="Courier",size=5,angle=90, vjust = -3,hjust = 0.5),
          text = ggplot2::element_text(size=5),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position="bottom")+
    ggplot2::labs(fill='ddG (kcal/mol)')+
    ggplot2::coord_fixed()
}
