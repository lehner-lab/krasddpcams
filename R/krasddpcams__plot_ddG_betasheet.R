
#' Plot allosteric mutation vs distance
#' 
#' This function allows you to plot binding interface mutations heatmap.
#' @param weighted_mean_ddG weighted mean free energy data
#' @param ddG free energy data
#' @param anno structure information annotation data
#' @param assay_name assay:"RAF","PI3"...
#' @param rect_input rect_input
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_ddG_betasheet<-function(
  weighted_mean_ddG=weighted_mean_ddG,
  ddG=ddG,
  anno=anno,
  assay_sele=assay_sele,
  rect_input=rect_input
  ){
  weighted_mean_ddG1<-weighted_mean_ddG
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
  data_plot[binding_type_gtp_included=="allosteric site"&mean>reg_threshold,site_type:="Novel allosteric site"]
  ddG<-fread(ddG)
  ddG[,Pos:=Pos_ref+1]
  data_plot_mutation<-merge(ddG,anno,by="Pos",all=T)
  data_plot_mutation[,mutation_type:="Reminder"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Binding interface site",Pos],mutation_type:="Orthosteric site"]
  data_plot_mutation[Pos%in%data_plot[site_type=="GTP binding interface site",Pos],mutation_type:="GTP binding site"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Novel allosteric site",Pos],mutation_type:="Novel allosteric site"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Reminder",Pos]&
                       (`mean_kcal/mol`>reg_threshold|`mean_kcal/mol`+reg_threshold<0)&
                       SASA<=0.2,mutation_type:="Core allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Reminder",Pos]&
                       (`mean_kcal/mol`>reg_threshold|`mean_kcal/mol`+reg_threshold<0)&
                       SASA>0.2,mutation_type:="Surface allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Reminder",Pos]&
                       (`mean_kcal/mol`>reg_threshold|`mean_kcal/mol`+reg_threshold<0)&
                       is.na(SASA),mutation_type:="HVR allosteric mutation"]
  data_plot_mutation<-within(data_plot_mutation,
                             mutation_type<-factor(mutation_type,
                                                   levels=c("GTP binding site",
                                                            "Novel allosteric site",
                                                            "Surface allosteric mutation",
                                                            "Core allosteric mutation",
                                                            "HVR allosteric mutation",
                                                            "Orthosteric site",
                                                            "Reminder")))
  rects_dt<-as.data.table(rect_input)
  data_plot_mutation[,colors_type:="others"]
  data_plot_mutation[Pos>=rects_dt[col=="b1",xstart]&Pos<=rects_dt[col=="b1",xend],colors_type:="b1"]
  data_plot_mutation[Pos>=rects_dt[col=="b2",xstart]&Pos<=rects_dt[col=="b2",xend],colors_type:="b2"]
  data_plot_mutation[Pos>=rects_dt[col=="b3",xstart]&Pos<=rects_dt[col=="b3",xend],colors_type:="b3"]
  data_plot_mutation[Pos>=rects_dt[col=="b4",xstart]&Pos<=rects_dt[col=="b4",xend],colors_type:="b4"]
  data_plot_mutation[Pos>=rects_dt[col=="b5",xstart]&Pos<=rects_dt[col=="b5",xend],colors_type:="b5"]
  data_plot_mutation[Pos>=rects_dt[col=="b6",xstart]&Pos<=rects_dt[col=="b6",xend],colors_type:="b6"]
  data_plot_mutation[,colors_type2:="others"]
  data_plot_mutation[colors_type!="others",colors_type2:="beta sheet"]
  data_plot_mutation[mutation_type=="GTP binding site",colors_type2:="GTP"]
  data_plot_mutation[mutation_type=="GTP binding site"&colors_type!="others",colors_type2:="GTP+beta sheet"]
  #data_plot_mutation
  # ggplot(data_plot_mutation,aes(x=scHAmin_ligand_RAF,y=`mean_kcal/mol`))+
  #   geom_point(aes(color=colors_type2),size=1)+
  #   geom_hline(yintercept = reg_threshold,linetype =2)+
  #   geom_hline(yintercept = -reg_threshold,linetype =2)+
  #   geom_vline(xintercept = 5,linetype =2)+
  #   scale_color_manual(values = c("red","blue","black"))+
  #   ylab("Binding free energy change \n(kcal/mol)")+
  #   xlab(expression(paste("Distance to binding parnter ("*ring(A)*")")))+
  #   theme_bw()+
  #   coord_fixed(ratio = 5)
  data_plot_mutation_beta<-data_plot_mutation[colors_type2=="beta sheet",]
  data_plot_mutation_beta<-within(data_plot_mutation_beta,
                                  colors_type,
                                  levels = c("b2","b3","b1","b4","b5","b6"))
  ggplot2::ggplot(data_plot_mutation_beta,ggplot2::aes(x=factor(colors_type, level = c("b2","b3","b1","b4","b5","b6")),y=`mean_kcal/mol`))+
    ggplot2::geom_violin()+
    ggplot2::geom_jitter(size=1*0.35,height = 0)+
    ggplot2::ylab("Binding free energy change \n(kcal/mol)")+
    ggplot2::xlab("beta sheet")+
    ggplot2::theme_classic()+
    ggplot2::theme(text = ggplot2::element_text(size = 7),
          axis.text = ggplot2::element_text(size = 7),legend.text = ggplot2::element_text(size = 7))
}