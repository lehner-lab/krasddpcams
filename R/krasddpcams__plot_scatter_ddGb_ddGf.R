
#' A Function to plot binding ddG against folding ddG colored by binding interface
#' 
#' This function allows you to plot binding ddG against folding ddG
#' @param ddG1 free energy data1
#' @param assay1 data1's name
#' @param ddG2 free energy data2
#' @param assay2 data2's name
#' @param anno all_distance_anno
#' @param colour_scheme colour scheme list
#' 
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_scatter_ddGb_ddGf<-function(
  ddG1=ddG1,
  assay1_sele=assay1_sele,
  ddG2=ddG2,
  assay2_sele=assay2_sele,
  anno=anno,
  colour_scheme
  ){
  ddG1<-krasddpcams__read_ddG(ddG = ddG1,
                 assay_sele = assay1_sele)
  ddG2<-krasddpcams__read_ddG(ddG = ddG2,
                 assay_sele = assay2_sele)
  all_ddG<-rbind(ddG1,ddG2)
  all_ddG_dc<-dcast(all_ddG[!is.na(mt),],mt+Pos_real~assay,value.var="mean_kcal/mol")
  all_ddG_dc_anno<-merge(all_ddG_dc,anno,by.x="Pos_real",by.y="Pos")
  all_ddG_dc_anno[,RAF_type_bs:="others"]
  all_ddG_dc_anno[scHAmin_ligand_RAF<=5,RAF_type_bs:="binding interface"]
  ggplot2::ggplot()+
    ggplot2::geom_point(data=all_ddG_dc_anno[Pos_real>1&RAF_type_bs=="others",],
               ggplot2::aes(x=folding,y=RAF1),color="black",alpha=0.6,size=0.1)+
    ggplot2::geom_point(data=all_ddG_dc_anno[Pos_real>1&RAF_type_bs=="binding interface"],
               ggplot2::aes(x=folding,y=RAF1),color=colour_scheme[["red"]],alpha=0.6,size=0.5)+
    ggpubr::theme_classic2()+
    ggplot2::labs(color=NULL)+
    ggplot2::xlab("Folding ddG")+
    ggplot2::ylab("Binding ddG")+
    ggplot2::theme(text = ggplot2::element_text(size=7),
          legend.position="right",
          legend.text = ggplot2::element_text(size=7),
          axis.text.x = ggplot2::element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5,margin=ggplot2::margin(0,0,0,0,"mm")),
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"),
          plot.margin=ggplot2::margin(0,0,0,0))+
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}