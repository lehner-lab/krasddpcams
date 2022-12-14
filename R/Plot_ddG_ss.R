#' Plot allosteric mutation in secondary structure
#' 
#' This function allows you to plot allosteric mutation in secondary structure.
#' @param weighted_mean_ddG weighted mean free energy data
#' @param ddG free energy data
#' @param anno structure information annotation data
#' @param assay_sele assay:"RAF","PI3"...
#' @param rect_input rects
#' @param rect_alpha rects_alpha_all
#' @return Nothing
#' @export
#' 
Plot_ddG_ss<-function(
  ddG=ddG,
  anno=anno,
  assay_sele=assay_sele,
  rect_input=rect_input,
  rect_alpha=rect_alpha
){
  ddG1<-Read_ddG(ddG = ddG,
                 assay_sele = assay_sele)
  
  weighted_mean_ddG1<-Get_weighted_mean_abs_ddG(ddG=ddG,assay_sele = assay_sele)
  weighted_mean_ddG1[,Pos:=Pos_real]
  data_plot<-merge(weighted_mean_ddG1,anno,by="Pos",all=T)
  data_plot[,binding_type:="allosteric site"]
  data_plot[get(paste0("scHAmin_ligand_",assay_sele))<5,binding_type:="binding site"]
  data_plot[,binding_type_gtp_included:=binding_type]
  data_plot[get(paste0("GXPMG_scHAmin_ligand_",assay_sele))<5,binding_type_gtp_included:="GTP binding site"]
  reg_threshold<-data_plot[binding_type=="binding site",sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean","sigma")]
  data_plot[,site_type:="Reminder"]
  data_plot[binding_type_gtp_included=="binding site",site_type:="Binding interface site"]
  data_plot[binding_type_gtp_included=="GTP binding site",site_type:="GTP binding interface site"]
  ddG1[,Pos:=Pos_real]
  data_plot_mutation1<-merge(ddG1,anno,by="Pos",all=T)
  data_plot_mutation<-data_plot_mutation1[Pos>1&!is.na(id),]
  data_plot_mutation[,mutation_type:="Reminder"]
  data_plot_mutation[,allosteric_mutation := p.adjust(doubledeepms__pvalue(abs(mean)-reg_threshold, std), method = "BH")<0.05 & (abs(mean)-reg_threshold)>0]
  data_plot_mutation[Pos%in%data_plot[site_type=="Binding interface site",Pos],mutation_type:="Orthosteric site"]
  data_plot_mutation[Pos%in%data_plot[site_type=="GTP binding interface site",Pos]&
                       allosteric_mutation==T,mutation_type:="GTP binding allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="GTP binding interface site",Pos]&
                       allosteric_mutation==F,mutation_type:="GTP binding other mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Reminder",Pos]&
                       allosteric_mutation==T&
                       SASA<=0.2,mutation_type:="Core allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Reminder",Pos]&
                       allosteric_mutation==T&
                       SASA>0.2,mutation_type:="Surface allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Reminder",Pos]&
                       allosteric_mutation==T&
                       is.na(SASA),mutation_type:="HVR allosteric mutation"]
  data_plot_mutation<-within(data_plot_mutation,
                             mutation_type<-factor(mutation_type,
                                                   levels=c("GTP binding allosteric mutation",
                                                            "GTP binding other mutation",
                                                            "Surface allosteric mutation",
                                                            "Core allosteric mutation",
                                                            "HVR allosteric mutation",
                                                            "Orthosteric site",
                                                            "Reminder")))
  rects_dt<-as.data.table(rect_input)
  rects_alphas<-as.data.table(rect_alpha)
  data_plot_mutation[,colors_type:="others"]
  data_plot_mutation[Pos>=rects_dt[col=="b1",xstart]&Pos<=rects_dt[col=="b1",xend],colors_type:="b1"]
  data_plot_mutation[Pos>=rects_dt[col=="b2",xstart]&Pos<=rects_dt[col=="b2",xend],colors_type:="b2"]
  data_plot_mutation[Pos>=rects_dt[col=="b3",xstart]&Pos<=rects_dt[col=="b3",xend],colors_type:="b3"]
  data_plot_mutation[Pos>=rects_dt[col=="b4",xstart]&Pos<=rects_dt[col=="b4",xend],colors_type:="b4"]
  data_plot_mutation[Pos>=rects_dt[col=="b5",xstart]&Pos<=rects_dt[col=="b5",xend],colors_type:="b5"]
  data_plot_mutation[Pos>=rects_dt[col=="b6",xstart]&Pos<=rects_dt[col=="b6",xend],colors_type:="b6"]
  data_plot_mutation[Pos>=rects_alphas[col=="a1",xstart]&Pos<=rects_alphas[col=="a1",xend],colors_type:="a1"]
  data_plot_mutation[Pos>=rects_alphas[col=="a2",xstart]&Pos<=rects_alphas[col=="a2",xend],colors_type:="a2"]
  data_plot_mutation[Pos>=rects_alphas[col=="a3",xstart]&Pos<=rects_alphas[col=="a3",xend],colors_type:="a3"]
  data_plot_mutation[Pos>=rects_alphas[col=="a4",xstart]&Pos<=rects_alphas[col=="a4",xend],colors_type:="a4"]
  data_plot_mutation[Pos>=rects_alphas[col=="a5",xstart]&Pos<=rects_alphas[col=="a5",xend],colors_type:="a5"]
  data_plot_mutation_final<-data_plot_mutation[mutation_type!="Orthosteric site",]
  ggplot2::ggplot()+
    ggplot2::geom_bar(data=data_plot_mutation_final,ggplot2::aes(x=..count..,fill=allosteric_mutation,y=colors_type,group=allosteric_mutation),position="stack",stat="count")+
    ggplot2::geom_text(data=data_plot_mutation_final,ggplot2::aes(x=..count..,label=..count..,y=colors_type,group=allosteric_mutation),position="stack",stat="count",size=7*0.35)+
    ggplot2::xlab("Mutations") +
    ggplot2::ylab("Secondary structure") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 7),
                   axis.text = ggplot2::element_text(size = 7),legend.text = ggplot2::element_text(size = 7))+
    ggplot2::ggtitle(assay_sele) +
    ggplot2::labs(fill = "Allosteric\nmutations")+
    ggplot2::scale_fill_manual(values = c("gray",colour_scheme[["red"]]),
                      labels=c("Others","Allosteric mutations")
    )
}