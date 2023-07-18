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
#' @return Nothing
#' @export
Plot_ddG_effector_bf_heatmap<-function(ddG1=ddG1,
                                       assay1=assay1,
                                       ddG2=ddG2,
                                       assay2=assay2,
                                       ddG3=ddG3,
                                       assay3=assay3,
                                       anno=anno,
                                       bind1=bind1,
                                       bind2=bind2,
                                       bind3=bind3){
  
  bind1<-substitute(bind1)
  bind2<-substitute(bind2)
  bind3<-substitute(bind3)
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
  ggplot()+
    geom_tile(data=output,mapping=aes(x=assay,y=mt_codon,fill=`mean_kcal/mol`))+
    geom_text(data=output[Pos_real>1&wt_codon==mt_codon,],
              mapping=aes(x=assay,y=mt_codon,label="-"),size=5*0.35)+
    scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value = "white")+
    facet_wrap(~Pos_wt,nrow = 1)+
    ylab("Mutant AA")+
    theme_bw()+
    theme(axis.text.x = element_text(size =7, angle=90, vjust=.5, hjust=1),
          axis.text.y = element_text(family="Courier New",size=5, vjust = 1,hjust = 1),
          text = element_text(size=7),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.text = element_text(size=7),
          strip.text.x = element_text(size=7),
          strip.background = element_rect(colour="black", fill="white"))+
    coord_fixed()
}