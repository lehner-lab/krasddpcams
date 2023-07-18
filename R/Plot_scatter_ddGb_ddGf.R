#' A Function to plot binding ddG against folding ddG colored by binding interface
#' 
#' This function allows you to plot binding ddG against folding ddG
#' @param ddG1 free energy data1
#' @param assay1 data1's name
#' @param ddG2 free energy data2
#' @param assay2 data2's name
#' @param anno all_distance_anno
#' 
#' @return Nothing
#' @export
Plot_scatter_ddGb_ddGf<-function(ddG1=ddG1,
                                      assay1_sele=assay1_sele,
                                      ddG2=ddG2,
                                      assay2_sele=assay2_sele,
                                      anno=anno){
  ddG1<-Read_ddG(ddG = ddG1,
                 assay_sele = assay1_sele)
  ddG2<-Read_ddG(ddG = ddG2,
                 assay_sele = assay2_sele)
  all_ddG<-rbind(ddG1,ddG2)
  all_ddG_dc<-dcast(all_ddG[!is.na(mt),],mt+Pos_real~assay,value.var="mean_kcal/mol")
  all_ddG_dc_anno<-merge(all_ddG_dc,anno,by.x="Pos_real",by.y="Pos")
  all_ddG_dc_anno[,RAF_type_bs:="others"]
  all_ddG_dc_anno[scHAmin_ligand_RAF<=5,RAF_type_bs:="binding interface"]
  ggplot()+
    geom_point(data=all_ddG_dc_anno[Pos_real>1&RAF_type_bs=="others",],
               aes(x=folding,y=RAF1),color="black",alpha=0.6,size=0.1)+
    geom_point(data=all_ddG_dc_anno[Pos_real>1&RAF_type_bs=="binding interface"],
               aes(x=folding,y=RAF1),color=colour_scheme[["red"]],alpha=0.6,size=0.5)+
    theme_classic2()+
    labs(color=NULL)+
    xlab("Folding ddG")+
    ylab("Binding ddG")+
    theme(text = element_text(size=7),
          legend.position="right",
          legend.text = element_text(size=7),
          axis.text.x = element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))+
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}