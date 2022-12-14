#' ddG heatmap Function
#' 
#' This function allows you to plot heatmap of free energy change.
#' @param input mochi data
#' @param assay_name folding or binding
#' 
#' @return Nothing
#' @export
Plot_heatmap_ddG<-function(input=input,assay_name=assay_name){
  input<-Read_ddG(ddG = input,assay_sele = assay_name)
  input[id!="WT",wt_codon:=substr(id,1,1)]
  input[id!="WT",mt_codon:=substr(id,nchar(id),nchar(id))]
  input[,mt:=paste0(wt_codon,Pos_real,mt_codon)]
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:188,each=20),
                           mt_codon = unlist(aa_list))
  input_heatmap<-merge(input,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  input_heatmap<-within(input_heatmap, mt_codon <- factor(mt_codon, levels = c("D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  input_heatmap[wt_codon==mt_codon,`mean_kcal/mol`:=0]
  ggplot2::ggplot()+
    ggplot2::theme_classic2()+
    ggplot2::geom_tile(data=input_heatmap[Pos_real>1,],ggplot2::aes(x=Pos_real,y=mt_codon,fill=`mean_kcal/mol`))+
    ggplot2::scale_x_discrete(limits=c(2:188),labels=c(2:188))+
    ggplot2::scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value ="white")+
    ggplot2::ggtitle(assay_name)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, vjust = 0.5,hjust = 0.5, 
                                     color = c(NA,NA,NA,rep(c("black",NA,NA,NA,NA),37))))+
    ggplot2::geom_text(data=input_heatmap[Pos_real>1&wt_codon==mt_codon,],
              ggplot2::aes(x=Pos_real,y=mt_codon),label="-",size=3)+
    ggplot2::labs(fill=NULL)+
    ggplot2::ylab("Mutant AA")+
    ggplot2::theme(text = ggplot2::element_text(size=5),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position="bottom",
          legend.text = ggplot2::element_text(size=5),
          axis.text.x = ggplot2::element_text(size =5, angle=90, vjust=.5, hjust=1),
          axis.text.y = ggplot2::element_text(family="Courier New",angle=90,size=5, vjust = .5,hjust = .5,margin=ggplot2::margin(0,-0.5,0,0,"mm")),
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"),
          plot.margin=ggplot2::margin(0,0,0,0))+
    ggplot2::coord_fixed()
}
