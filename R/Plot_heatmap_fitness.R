#' fitness heatmap Function
#' 
#' This function allows you to plot heatmap of fitness.
#' @param input dimsum data
#' @param assay_name folding or binding
#' 
#' @return Nothing
#' @export
Plot_heatmap_fitness<-function(input=input,assay_name=assay_name,anno_input=anno_input){
  input_assay<-input[assay==assay_name,]
  input_single <- Nor_overlap_single_mt_fitness(input_assay)
  input_single_pos <-input_single
  input_single_pos[,position:=AA_Pos1]
  input_single_pos[,WT_AA:=wtcodon1]
  heatmap_tool_fitness<-data.table(wtcodon1 = rep(unlist(strsplit(wt_aa,"")),each=21),
                                   position = rep(2:188,each=21),
                                   codon1 = c(unlist(aa_list),"*"))
  heatmap_tool_fitness_anno_single<-merge(input_single_pos,heatmap_tool_fitness,by=c("wtcodon1","position","codon1"),all=T)
  heatmap_tool_fitness_anno_single<-within(heatmap_tool_fitness_anno_single, 
                                           codon1 <- factor(codon1, 
                                                            levels = c("*","D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  heatmap_tool_fitness_anno_single[wtcodon1==codon1,nor_fitness_nooverlap:=0]
  heatmap_text_anno<-data.table(Pos=1:188,aa=c("M",unlist(strsplit(wt_aa,""))),
                                RAF_bsite=1,K27_bsite=1,PI3_bsite=1,RAL_bsite=1,SOS_bsite=1,
                                K55_bsite=1,SOS1_bsite=1,SOS2_bsite=1)
  anno<-as.data.table(anno_input)                              
  heatmap_text_anno<-merge(heatmap_text_anno,anno,by="Pos",all=T)
  heatmap_text_anno[scHAmin_ligand_RAF<5,RAF_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_K27<5,K27_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_PI3<5,PI3_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_RAL<5,RAL_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_SOS<5,SOS_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_K55<5,K55_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_SOS1<5,SOS1_bsite:=2]
  heatmap_text_anno[scHAmin_ligand_SOS2<5,SOS2_bsite:=2]
  heatmap_text_anno[,second_structure_hm:=1]
  heatmap_text_anno[second_structure%in%c("a1","a2","a3","a4","a5"),second_structure_hm:=2]
  heatmap_text_anno[second_structure%in%c("b1","b2","b3","b4","b5","b6"),second_structure_hm:=3]
  heatmap_text_anno[second_structure%in%c("HVR"),second_structure_hm:=4]
  ggplot2::ggplot()+
    ggplot2::theme_classic()+
    ggplot2::geom_tile(data=heatmap_tool_fitness_anno_single[position>1,],ggplot2::aes(x=position,y=codon1,fill=nor_fitness_nooverlap))+
    ggplot2::scale_x_discrete(limits=c(2:188),labels=c(2:188))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, vjust = 0.5,hjust = 0.5, 
                                     color = c(NA,NA,NA,rep(c("black",NA,NA,NA,NA),37))))+
    ggplot2::scale_fill_gradient2(low=colour_scheme[["red"]],mid="gray",high=colour_scheme[["blue"]],na.value = "white")+
    ggplot2::ylab("Mutant aa")+
    ggplot2::ggtitle(paste(assay_name,"fitness"))+
    ggplot2::labs(fill=NULL)+
    ggplot2::geom_text(data=heatmap_tool_fitness_anno_single[position>1&wtcodon1==codon1,],ggplot2::aes(x=position,y=codon1),label="-",size=3)+
    ggplot2::theme(text = ggplot2::element_text(size=5),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position="bottom",
          legend.text = ggplot2::element_text(size=5),
          axis.text.x = ggplot2::element_text(size =5, angle=90, vjust=.5, hjust=1),
          axis.text.y = ggplot2::element_text(family="Courier",size=5, vjust = 1,hjust = -10,margin=ggplot2::margin(0,0,0,0)),
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"),
          plot.margin=ggplot2::margin(0,0,0,0))+
    ggplot2::coord_fixed()
}


