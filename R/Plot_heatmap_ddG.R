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
  ggplot()+
    theme_classic2()+
    geom_tile(data=input_heatmap[Pos_real>1,],aes(x=Pos_real,y=mt_codon,fill=`mean_kcal/mol`))+
    scale_x_discrete(limits=c(2:188),labels=c(2:188))+
    scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value ="white")+
    ggtitle(assay_name)+
    theme(axis.text.x = element_text(size = 5, vjust = 0.5,hjust = 0.5, 
                                     color = c(NA,NA,NA,rep(c("black",NA,NA,NA,NA),37))))+
    geom_text(data=input_heatmap[Pos_real>1&wt_codon==mt_codon,],
              aes(x=Pos_real,y=mt_codon),label="-",size=3)+
    labs(fill=NULL)+
    ylab("Mutant AA")+
    theme(text = element_text(size=5),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.text = element_text(size=5),
          axis.text.x = element_text(size =5, angle=90, vjust=.5, hjust=1),
          axis.text.y = element_text(family="Courier New",angle=90,size=5, vjust = .5,hjust = .5,margin=margin(0,-0.5,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))+
    coord_fixed()
}
