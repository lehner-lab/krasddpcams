#' RAF in vitro correlation Function
#' 
#' This function allows you to plot RAF1 binding free energy Mochi inferred against in vitro data.
#' @param input RAF1 free energy data
#' @param assay_name assay
#' 
#' @return Nothing
#' @export
Plot_RAF_invitro_cor<-function(input=input,
                                    assay_name=assay_name){
  input_ddG<-Read_ddG(ddG = input,assay_sele = assay_name)
  christina_b_raf_energy<-data.table(id.Raf.RBD.christina=c("WT","R59A","N64A","Q66A","R67A","T68A","V69A","K84A","V88A",rep("WT",10)),
                                     id.Rasprotein.christina=c(rep("WT",9),"I21A","H27A","E31A","D33A","I36A","E37A","D38A","S39A","R41A","V45A"),
                                     G_ITC=c(-10,-8.3,-9.9,-8.3,-8.0,-8.5,-9.0,-7.8,-9.9,-9.1,-10.3,-9.5,-8.9,-8.2,-8.2,-6.7,-9.7,-10.6,-10.3)
  )
  christina_mochi<-input_ddG[christina_b_raf_energy,on=.(mt=id.Rasprotein.christina)]
  christina_mochi[mt=="WT"&`id.Raf.RBD.christina`=="WT",`mean_kcal/mol`:=0]
  christina_mochi_m<-lm(`mean_kcal/mol`~G_ITC,christina_mochi)
  christina_mochi[,ddG_ITC:=G_ITC-christina_mochi[id.Raf.RBD.christina=="WT"&mt=="WT",G_ITC]]
  ggplot2::ggplot(christina_mochi,ggplot2::aes(x=ddG_ITC,y=`mean_kcal/mol`,label=mt))+
    ggplot2::geom_smooth(method=lm, se=T,size=0.1,color=colour_scheme[["red"]])+
    ggplot2::geom_point(size=0.1)+
    ggplot2::xlab("Binding \u2206\u2206G to RAF1 (in vitro)\n(kcal/mol)")+
    ggplot2::ylab("Binding \u2206\u2206G to RAF1 (inferred)\n(kcal/mol)")+
    ggplot2::geom_pointrange(ggplot2::aes(ymin=`mean_kcal/mol`-`std_kcal/mol`, ymax=`mean_kcal/mol`+`std_kcal/mol`),size=0.1)+
    ggplot2::geom_text_repel(size=7*0.35)+
    ggplot2::annotate("text",x=0,y=2.0,,label = paste0("r = ",round(sqrt(summary(christina_mochi_m)$r.squared),2)) ,size=7*0.35)+
    ggplot2::theme_classic2()+
    ggplot2::theme(text = ggplot2::element_text(size=7),
          legend.position="right",
          legend.text = ggplot2::element_text(size=7),
          axis.text.x = ggplot2::element_text(size =7,vjust=.5, hjust=1),
          axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5,margin=ggplot2::margin(0,0,0,0,"mm")),
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"),
          plot.margin=ggplot2::margin(0,0,0,0))+
    ggplot2::coord_fixed()
}
