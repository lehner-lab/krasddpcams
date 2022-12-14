#' plot chosen position ddG heatmap
#' 
#' This function allows you to plot chosen position ddG heatmap.
#' @param ddG1 free energy data set assay1
#' @param assay1 assay1
#' @param ddG2 free energy data set assay2
#' @param assay2 assay2
#' @param ddG3 free energy data set assay3
#' @param assay3 assay3
#' @param ddG4 free energy data set assay4
#' @param assay4 assay4
#' @param ddG5 free energy data set assay5
#' @param assay5 assay5
#' @param ddG6 free energy data set assay6
#' @param assay6 assay6
#' @param anno anno_final3 
#' @param bind1 bind1
#' @param bind2 bind2
#' @param bind3 bind3
#' @param bind4 bind4
#' @param bind5 bind5
#' @param bind6 bind6
#' @return Nothing
#' @export
Plot_ddG_all_bf_heatmap<-function(ddG1=ddG1,
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
                                       anno=anno,
                                       bind1=bind1,
                                       bind2=bind2,
                                       bind3=bind3,
                                       bind4=bind4,
                                       bind5=bind5,
                                       bind6=bind6){
  
  bind1<-substitute(bind1)
  bind2<-substitute(bind2)
  bind3<-substitute(bind3)
  bind4<-substitute(bind4)
  bind5<-substitute(bind5)
  bind6<-substitute(bind6)
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:188,each=20),
                           mt_codon = unlist(aa_list))
  ddG1<-Read_ddG(ddG1,assay1)
  input1_heatmap<-merge(ddG1,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG2<-Read_ddG(ddG2,assay2)
  input2_heatmap<-merge(ddG2,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG3<-Read_ddG(ddG3,assay3)
  input3_heatmap<-merge(ddG3,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG4<-Read_ddG(ddG4,assay4)
  input4_heatmap<-merge(ddG4,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG5<-Read_ddG(ddG5,assay5)
  input5_heatmap<-merge(ddG5,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  ddG6<-Read_ddG(ddG6,assay6)
  input6_heatmap<-merge(ddG6,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  
  ddG<-rbind(input1_heatmap,input2_heatmap,input3_heatmap,input4_heatmap,input5_heatmap,input6_heatmap)
  ddG[wt_codon==mt_codon,`mean_kcal/mol`:=0]
  output<-merge(ddG,anno,by.x="Pos_real",by.y="Pos",all=T)
  output<-output[Pos_real>1&(eval(bind1)<5|eval(bind2)<5|eval(bind3)<5|eval(bind4)<5|eval(bind5)<5|eval(bind6)<5),]
  output<-within(output, 
                 mt_codon <- factor(mt_codon, 
                                    levels = c("D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  output<-within(output, 
                 assay <- factor(assay, 
                                 levels = c("RAF1","PIK3CG","RALGDS","SOS1","DARPin K27","DARPin K55")))
  output[,Pos_wt:=paste0(wt_codon,Pos_real)]
  output<-within(output, 
                 Pos_wt <- factor(Pos_wt,
                                  levels = unique(output[order(Pos_real),Pos_wt])))
  ggplot2::ggplot()+
    ggplot2::geom_tile(data=output,mapping=ggplot2::aes(x=assay,y=mt_codon,fill=`mean_kcal/mol`))+
    ggplot2::geom_text(data=output[Pos_real>1&wt_codon==mt_codon,],
              mapping=ggplot2::aes(x=assay,y=mt_codon,label="-"),size=5*5/14)+
    ggplot2::scale_fill_gradient2(low=colour_scheme[["blue"]],mid="grey",high=colour_scheme[["red"]],na.value = "white")+
    ggplot2::facet_wrap(~Pos_wt,nrow = 3)+
    ggplot2::ylab("Mutant AA")+
    ggplot2::theme_bw()+
    ggplot2::theme(text = ggplot2::element_text(size=5),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position="bottom",
          legend.text = ggplot2::element_text(size=5),
          strip.text.x = ggplot2::element_text(size=7),
          strip.background = ggplot2::element_rect(colour="black", fill="white"),
          axis.text.x = ggplot2::element_text(size =5, angle=90, vjust=.5, hjust=1),
          axis.text.y = ggplot2::element_text(family="Courier New",size=5, vjust = 1,hjust = -10,margin=ggplot2::margin(0,0,0,0)),
          panel.spacing.x =  ggplot2::unit(1, "mm"),
          panel.spacing.y =  ggplot2::unit(3, "mm"),
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"))+
    ggplot2::coord_fixed()
}