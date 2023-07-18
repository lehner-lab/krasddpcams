#' fitness distribution
#' 
#' This function allows you to plot fitness distribution.
#' @param input dimsum data
#' @param assay_name folding or binding
#' @param anno_input anno_input
#' 
#' @return Nothing
#' @export
Plot_distribution_fitness<-function(input=input,assay_name=assay_name,anno_input=anno_input){
  input_assay<-input[assay==assay_name,]
  input_single<-Nor_overlap_single_mt_fitness(input_assay)
  input_single_pos <-Pos_id(input_single,wt_aa = wt_aa)
  input_single_pos[,position:=AA_Pos1]
  input_single_pos[,WT_AA:=wtcodon1]
  #input_single_pos[,fitness_fdr:=p.adjust(doubledeepms__pvalue(nor_fitness,nor_fitness_sigma),method="BH")]
  output<-merge(input_single_pos,anno_input,by.x="AA_Pos1",by.y="Pos")
  output[,class:="others"]
  output[SASA<0.25,class:="core"]
  output[get(paste0("scHAmin_ligand_",assay_name))<5,class:="binding interface"]
  output[,type:="no differences"]
  #output[nor_fitness<0&fitness_fdr<0.05,type:="Decrease binding fitness"]
  #output[nor_fitness>0&fitness_fdr<0.05,type:="Increase binding fitness"]
  #output<-within(output,class<-factor(class,levels=c("binding interface","core","others")))
  #output<-within(output,type<-factor(type,levels=c("Increase binding fitness","Decrease binding fitness","no differences")))
  ggplot()+
    geom_density(data=output,aes(x=nor_fitness_nooverlap,color=class))+
    ggplot2::xlab("RAF1 Binding") +
    ggplot2::ylab("Density") +
    ggplot2::theme_classic() +
    scale_color_manual(values = c(colour_scheme[["red"]],"black","gray"),
                       labels=c("binding interface","core","others"))+
    ggplot2::labs(color=NULL)+
    theme(text = element_text(size=7),
          legend.position="right",
          legend.text = element_text(size=7),
          axis.text.x = element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))
  
}