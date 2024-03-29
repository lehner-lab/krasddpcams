
#' plot chosen position ddG heatmap
#' 
#' This function allows you to plot chosen position ddG heatmap.
#' @param ddG1 free energy data set assay1
#' @param assay1 assay1
#' @param ddG2 free energy data set assay2
#' @param assay2 assay2
#' @param ddG3 free energy data set assay3
#' @param assay3 assay23
#' @param anno anno_final2 are chosen 
#' @param bind1 bind1, scHAmin_ligand_RAF
#' @param bind2 bind2, scHAmin_ligand_PI3
#' @param bind3 bind3, scHAmin_ligand_RAL
#' @param wt_aa wt amino acid sequence
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_ddG_effector_bf_heatmap<-function(
  ddG1=ddG1,
  assay1=assay1,
  ddG2=ddG2,
  assay2=assay2,
  ddG3=ddG3,
  assay3=assay3,
  anno=anno,
  bind1=bind1,
  bind2=bind2,
  bind3=bind3,
  wt_aa,
  colour_scheme
  ){
  
  bind1<-substitute(bind1)
  bind2<-substitute(bind2)
  bind3<-substitute(bind3)
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:188,each=20),
                           mt_codon = unlist(aa_list))
  ddG1<-krasddpcams__read_ddG(ddG1,assay1)
  input1_heatmap<-merge(ddG1,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG2<-krasddpcams__read_ddG(ddG2,assay2)
  input2_heatmap<-merge(ddG2,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG3<-krasddpcams__read_ddG(ddG3,assay3)
  input3_heatmap<-merge(ddG3,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG<-rbind(input1_heatmap,input2_heatmap,input3_heatmap)
  ddG[wt_codon==mt_codon,`mean_kcal/mol`:=0]
  output<-merge(ddG,anno,by.x="Pos_real",by.y="Pos",all=T)
  output<-output[eval(bind1)<5|eval(bind2)<5|eval(bind3)<5,]
  output[,Pos_wt:=paste0(wt_codon,Pos_real)]
  output<-within(output, 
                 mt_codon <- factor(mt_codon, 
                                    levels = c("D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  output<-within(output, 
                 Pos_wt <- factor(Pos_wt,
                                  levels = unique(output[order(Pos_real),Pos_wt])))
  output<-within(output, 
                 assay <- factor(assay, 
                                    levels = c("RAF1","PIK3CG","RALGDS")))
  ggplot2::ggplot()+
    ggplot2::geom_tile(data=output,mapping=ggplot2::aes(x=assay,y=mt_codon,fill=`mean_kcal/mol`))+
    ggplot2::geom_text(data=output[Pos_real>1&wt_codon==mt_codon,],
              mapping=ggplot2::aes(x=assay,y=mt_codon,label="-"),size=5*0.35)+
    ggplot2::scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value = "white")+
    ggplot2::facet_wrap(~Pos_wt,nrow = 1)+
    ggplot2::ylab("Mutant AA")+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size =7, angle=90, vjust=.5, hjust=1),
          axis.text.y = ggplot2::element_text(family="Courier",size=5, vjust = 1,hjust = 1),
          text = ggplot2::element_text(size=7),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position="bottom",
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.text = ggplot2::element_text(size=7),
          strip.text.x = ggplot2::element_text(size=7),
          strip.background = ggplot2::element_rect(colour="black", fill="white"))+
    ggplot2::coord_fixed()
}