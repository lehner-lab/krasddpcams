#' fitness heatmap Function for block1 only
#' 
#' This function allows you to plot heatmap of fitness for block1.
#' @param input dimsum data
#' @param assay_name folding or binding
#' 
#' @return Nothing
#' @export
Plot_heatmap_fitness_block1<-function(input=input){
  input_assay<-input
  input_single <- Nor_overlap_single_mt_fitness(input_assay)
  input_single_pos <-input_single
  input_single_pos[,position:=AA_Pos1]
  input_single_pos[,WT_AA:=wtcodon1]
  heatmap_tool_fitness<-data.table(wtcodon1 = rep(unlist(strsplit(wt_aa,"")),each=21),
                                   position = rep(2:188,each=21),
                                   codon1 = c(unlist(aa_list),"*"))
  heatmap_tool_fitness_block1<-heatmap_tool_fitness[position<65,]
  heatmap_tool_fitness_anno_single<-merge(input_single_pos,heatmap_tool_fitness_block1,by=c("wtcodon1","position","codon1"),all=T)
  heatmap_tool_fitness_anno_single<-within(heatmap_tool_fitness_anno_single, 
                                           codon1 <- factor(codon1, 
                                                            levels = c("*","D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  heatmap_tool_fitness_anno_single[wtcodon1==codon1,nor_fitness_nooverlap:=0]
  ggplot()+
    theme_classic()+
    geom_tile(data=heatmap_tool_fitness_anno_single[position>1,],aes(x=position,y=codon1,fill=nor_fitness_nooverlap))+
    scale_x_discrete(limits=c(2:188),labels=c(2:188))+
    theme(axis.text.x = element_text(size = 8, vjust = 0.5,hjust = 0.5, 
                                     color = c(NA,NA,NA,rep(c("black",NA,NA,NA,NA),37))))+
    scale_fill_gradient2(low=colour_scheme[["red"]],mid="gray",high=colour_scheme[["blue"]],na.value = "white")+
    ylab("Mutant aa")+
    ggtitle(paste("fitness"))+
    labs(fill=NULL)+
    geom_text(data=heatmap_tool_fitness_anno_single[position>1&wtcodon1==codon1,],aes(x=position,y=codon1),label="-",size=3)+
    theme(text = element_text(size=5),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.text = element_text(size=5),
          axis.text.x = element_text(size =5, angle=90, vjust=.5, hjust=1),
          axis.text.y = element_text(family="Courier New",size=5, vjust = 1,hjust = -10,margin=margin(0,0,0,0)),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))+
    coord_fixed()
}
