#' RALGDS in vitro cor Function
#' 
#' This function allows you to plot RAF1 binding free energy Mochi inferred against in vitro data.
#' @param ddG RALGDS free energy data
#' @param assay_sele assay_sele
#' @return Nothing
#' @export
Plot_RAL_invitro_cor_GDI<-function(ddG=ddG,
                                    assay_sele=assay_sele){
  ddG<-Read_ddG(ddG,assay_sele)
  christina_b_ralgds_energy<-data.table(id.RalGDS.RBD.christina=c("WT","I14A","R16A","N23A","N25A","Y27F","Y27A","K28A","S29A","K44A","D47A","K48A","H49A","E52A","E53A",rep("WT",12)),
                                        id.Rasprotein.christina=c(rep("WT",15),"Q25A","V29A","E31A","D33A","E37A","D38A","S39A","Y40F","Y40A","R41A","E62A","E63A"),
                                        G_ITC=c(-8.4,-7.4,-6.2,-7.4,-8.1,-7.4,NA,-6.6,-7.6,-8.1,-8.7,-6.0,-7.4,-8.6,-8.5,-7.5,NA,-8.2,-7.3,NA,NA,-9.0,-6.1,NA,-7.5,NA,NA),
                                        G_GDI=c(-8.4,-6.5,-6.0,-7.5,-7.7,-6.7,-4.8,-5.2,-7.3,-7.8,-8.7,-5.2,-6.4,NA,-8.2,NA,-7.9,-8.1,NA,-7.2,-4.5,-9.3,-6.4,-4.7,-7.6,-8.2,-8.1))
  christina_mochi_RAL<-ddG[christina_b_ralgds_energy,on=.(mt=id.Rasprotein.christina)]
  christina_mochi_RAL[mt=="WT"&`id.RalGDS.RBD.christina`=="WT",`mean_kcal/mol`:=0]
  christina_mochi_RAL_m_ITC<-lm(`mean_kcal/mol`~G_ITC,christina_mochi_RAL)
  christina_mochi_RAL_m_GDI<-lm(`mean_kcal/mol`~G_GDI,christina_mochi_RAL)
  christina_mochi_RAL[,ddG_ITC:=G_ITC-christina_mochi_RAL[id.RalGDS.RBD.christina=="WT"&mt=="WT",G_ITC]]
  christina_mochi_RAL[,ddG_GDI:=G_GDI-christina_mochi_RAL[id.RalGDS.RBD.christina=="WT"&mt=="WT",G_GDI]]
  
    ggplot2::ggplot(christina_mochi_RAL,ggplot2::aes(x=ddG_GDI,y=`mean_kcal/mol`,label=mt))+
      ggplot2::geom_smooth(method=lm, se=T,color=colour_scheme[["red"]])+
      ggplot2::geom_point(size=0.1)+
      ggplot2::xlab("Binding ddG to RALGDS (in vitro/GDI)\n(kcal/mol)")+
      ggplot2::ylab("Binding ddG to RALGDS (inferred)\n(kcal/mol)")+
      ggplot2::geom_pointrange(aes(ymin=`mean_kcal/mol`-`std_kcal/mol`, ymax=`mean_kcal/mol`+`std_kcal/mol`),size=0.1)+
      ggplot2::geom_text_repel(hjust=0.5, vjust=2,size=7*0.35)+
      ggplot2::annotate("text",x=0,y=2,
               label = paste0("r = ",round(sqrt(summary(christina_mochi_RAL_m_GDI)$r.squared
                                                ),3
                                           )
                              ),
               size=7*0.35 )+
      ggplot2::theme_classic2()+
      ggplot2::ylim(-1,2.5)+
      ggplot2::theme(text = ggplot2::element_text(size = 7),
            axis.text = ggplot2::element_text(size = 7),legend.text = ggplot2::element_text(size = 7))+
      ggplot2::coord_fixed()
  
}