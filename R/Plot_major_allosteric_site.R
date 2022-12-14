#' Plot major allosteric sites
#' 
#' This function allows you to plot major allosteric sites.
#' @param ddG free energy data
#' @param anno structure information annotation data
#' @param assay_name assay:"RAF","PI3"...
#' @param allosteric_list output data with allosteric sites and ggplot
#' @param rect_input rect_input
#' @return list of allosteric sites information
#' @export
#' 
Plot_major_allosteric_site<-function(ddG=ddG,
                                          anno=anno,
                                          assay_name=assay_name,
                                          allosteric_list=allosteric_list,
                                          rect_input=rect_input){
  weighted_mean_ddG<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG,
                                                    assay_sele = assay_name)
  data_plot<-merge(weighted_mean_ddG,anno,by="Pos",all=T)
  data_plot[,binding_type:="allosteric site"]
  data_plot[get(paste0("scHAmin_ligand_",assay_name))<5,binding_type:="binding site"]
  data_plot[,binding_type_gtp_included:=binding_type]
  reg_threshold<-data_plot[binding_type=="binding site",sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean","sigma")]
  data_plot[get(paste0("GXPMG_scHAmin_ligand_",assay_name))<5,binding_type_gtp_included:="GTP binding site"]
  data_plot[,site_type:="Reminder"]
  data_plot[binding_type_gtp_included=="binding site",site_type:="Binding interface site"]
  data_plot[binding_type_gtp_included=="GTP binding site",site_type:="Other GTP pocket site"]
  data_plot[binding_type_gtp_included=="GTP binding site"&mean>reg_threshold&get(paste0("scHAmin_ligand_",assay_name))>=5,site_type:="Allosteric GTP pocket site"]
  data_plot[binding_type_gtp_included=="allosteric site"&mean>reg_threshold,site_type:="Major allosteric site"]
  data_plot[,colors_type:="others"]
  rects_dt<-as.data.table(rect_input)
  data_plot[Pos>=rects_dt[col=="b1",xstart]&Pos<=rects_dt[col=="b1",xend],colors_type:="b1"]
  data_plot[Pos>=rects_dt[col=="b2",xstart]&Pos<=rects_dt[col=="b2",xend],colors_type:="b2"]
  data_plot[Pos>=rects_dt[col=="b3",xstart]&Pos<=rects_dt[col=="b3",xend],colors_type:="b3"]
  data_plot[Pos>=rects_dt[col=="b4",xstart]&Pos<=rects_dt[col=="b4",xend],colors_type:="b4"]
  data_plot[Pos>=rects_dt[col=="b5",xstart]&Pos<=rects_dt[col=="b5",xend],colors_type:="b5"]
  data_plot[Pos>=rects_dt[col=="b6",xstart]&Pos<=rects_dt[col=="b6",xend],colors_type:="b6"]
  data_plot[,shape_beta:="others"]
  data_plot[colors_type!="others",shape_beta:="beta strand"]
  data_plot<-data_plot[Pos_real>1&count>9.5,]
  data_plot<-within(data_plot,
                    site_type<-factor(site_type,
                                      levels = c("Binding interface site",
                                                 "Allosteric GTP pocket site",
                                                 "Other GTP pocket site",
                                                 "Major allosteric site",
                                                 "Reminder")))
  allosteric_list<-list()
  allosteric_list[["Binding interface site"]]<-data_plot[get(paste0("scHAmin_ligand_",assay_name))<5,Pos]
  allosteric_list[["Allosteric GTP pocket site"]]<-data_plot[site_type=="Allosteric GTP pocket site",Pos]
  allosteric_list[["Other GTP pocket site"]]<-data_plot[site_type=="Other GTP pocket site",Pos]
  allosteric_list[["Major allosteric site"]]<-data_plot[site_type=="Major allosteric site",Pos]
  allosteric_list[["p"]]<-ggplot2::ggplot()+
    ggplot2::geom_point(data=data_plot,
               mapping = ggplot2::aes(x=get(paste0("scHAmin_ligand_",assay_name)),
                             y=mean,
                             color=site_type,
                             shape=as.factor(shape_beta)),
               size=1,alpha=1/3)+
    ggplot2::scale_color_manual(values=c(colour_scheme[["red"]], 
                                colour_scheme[["blue"]],
                                colour_scheme[["light blue"]],
                                colour_scheme[["green"]],
                                "gray"),
                       labels=c(paste("Binding interface",
                                      "=",length(data_plot[get(paste0("scHAmin_ligand_",assay_name))<5,Pos])),
                                paste("Allosteric GTP pocket",
                                      "=",length(data_plot[site_type=="Allosteric GTP pocket site",Pos])),
                                paste("Other GTP pocket",
                                      "=",length(data_plot[site_type=="Other GTP pocket site",Pos])),
                                paste("Major allosteric",
                                      "=",length(data_plot[site_type=="Major allosteric site",Pos])),
                                paste("Others",
                                      "=",length(data_plot[site_type=="Reminder",Pos]))
                       ))+
    ggplot2::geom_pointrange(data=data_plot,ggplot2::aes(x=get(paste0("scHAmin_ligand_",assay_name)),
                                       y=mean,
                                       color=site_type,ymin=mean-sigma, ymax=mean+sigma,
                                       shape=as.factor(shape_beta)),size=0.35)+
    ggplot2::geom_hline(yintercept = reg_threshold,linetype =2,size=0.1)+
    ggplot2::geom_vline(xintercept = 5,linetype =2,size=0.1)+
    ggplot2::geom_text_repel(data=data_plot[site_type=="Major allosteric site",],ggplot2::aes(x=get(paste0("scHAmin_ligand_",assay_name)),
                                                                                                              y=mean,label=Pos),nudge_y=0.05,color=colour_scheme[["green"]],size=7*0.35)+
    ggplot2::geom_text_repel(data=data_plot[site_type=="Allosteric GTP pocket site",],ggplot2::aes(x=get(paste0("scHAmin_ligand_",assay_name)),
                                                                                                                    y=mean,label=Pos),nudge_y=0.05,color=colour_scheme[["blue"]],size=7*0.35)+
    ylab(paste("Weighted mean |ddG|(kcal/mol)"))+
    xlab(expression(paste("Distance to binding partner ("*ring(A)*")")))+
    ggplot2::labs(color=NULL)+
    ggplot2::theme_classic2()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust =.5),
          text = ggplot2::element_text(size=7),
          legend.position="right",
          #legend.key.height= unit(3.1, 'mm'),
          #legend.key.width = unit(5, 'mm'),
          legend.text = ggplot2::element_text(size=7),plot.margin=ggplot2::margin(0,1,0,1,"mm"),
          legend.margin=ggplot2::margin(0,0,0,l=-2,"mm"),
          legend.spacing.y = ggplot2::unit(0, 'mm'),
          legend.key.height=ggplot2::unit(4,"mm"))+
    ggplot2::labs(shape=NULL)+
    ggplot2::coord_cartesian(xlim = c(0,34),ylim=c(0,2.3),)
  return(allosteric_list)
}