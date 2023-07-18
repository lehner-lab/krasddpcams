#' Plot binding interface mutation heatmap
#' 
#' This function allows you to plot binding interface mutations heatmap.
#' @param ddG free energy data
#' @param assay_sele assay_sele
#' @return Nothing
#' @export
Plot_binding_interface_mutation_heatmap<-function(ddG=ddG,
                                                  assay_sele=assay_sele){
  ddG<-Read_ddG(ddG,assay_sele)
  ddG<-ddG[id!="WT",]
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:188,each=20),
                           mt_codon = unlist(aa_list))
  input_heatmap<-merge(ddG,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  input_heatmap<-within(input_heatmap, 
              mt_codon <- factor(mt_codon, 
                                 levels = c("D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  ddG_plot<-input_heatmap[Pos_real%in%c(31,33,36,37,38,39,40,24,25,41,3),]
  ddG_plot[,wt_Pos:=paste0(wt_codon,Pos_real)]
  ddG_plot<-within(ddG_plot, 
                   wt_Pos <- factor(wt_Pos, 
                                    levels = c("E31","D33","I36","E37","D38","S39","Y40","I24","Q25","R41","E3")))
  ddG_plot[,assay:="RAF1"]
  ddG_plot[wt_codon==mt_codon,`mean_kcal/mol`:=0]

  ggplot()+
    geom_tile(data=ddG_plot,
              aes(y=assay,x=mt_codon,fill=`mean_kcal/mol`))+
    scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value = "white")+
    geom_text(data=ddG_plot[Pos_real>1&wt_codon==mt_codon,],
              aes(x=mt_codon,y=assay),label="-",size=5*5/14)+
    facet_wrap(~wt_Pos,ncol=1)+
    labs(fill="\u0394\u0394G(kcal/mol)",)+
    xlab("Mutant aa")+
    ylab(NULL)+
    theme_classic2()+
    ggplot2::theme(text = element_text(size=7),
                   axis.ticks.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   legend.position="right",
                   legend.text = element_text(size=7),
                   axis.text.x = element_text(size =5, vjust=.5, hjust=.5),
                   axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
                   legend.key.height= unit(3.1, 'mm'),
                   legend.key.width = unit(3.1, 'mm'),
                   legend.key.size = unit(1,"mm"),
                   plot.margin=margin(0,0,0,0),
                   strip.background = element_rect(color="white",fill=NULL,linetype = NULL))+
    coord_fixed()
}