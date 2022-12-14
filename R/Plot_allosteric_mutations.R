#' Plot allosteric mutations
#' 
#' This function allows you to plot allosteric mutations.
#' @param ddG free energy data
#' @param anno structure information annotation data
#' @param assay_name assay:"RAF","PI3"...
#' @param rect_input rect_input
#' @param rect_alpha rect_alpha
#' @return Nothing
#' @export
#' 
Plot_allosteric_mutations<-function(
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
  data_plot_mutation[Pos%in%data_plot[site_type=="Binding interface site",Pos]&
                       allosteric_mutation==T,mutation_type:="Orthosteric site huge differences"]
  data_plot_mutation[Pos%in%data_plot[site_type=="Binding interface site",Pos]&
                       allosteric_mutation==F,mutation_type:="Orthosteric site small differences"]
  data_plot_mutation[Pos%in%data_plot[site_type=="GTP binding interface site",Pos]&
                       allosteric_mutation==T,mutation_type:="GTP binding allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[site_type=="GTP binding interface site",Pos]&
                       allosteric_mutation==F,mutation_type:="GTP binding other mutation"]
  data_plot_mutation[Pos%in%data_plot[!site_type%in%c("GTP binding interface site","Binding interface site"),Pos]&
                       allosteric_mutation==T,mutation_type:="Allosteric mutation"]
  data_plot_mutation[Pos%in%data_plot[!site_type%in%c("GTP binding interface site","Binding interface site"),Pos]&
                       allosteric_mutation==F,mutation_type:="Other mutation"]
  
  data_plot_mutation<-within(data_plot_mutation,
                             mutation_type<-factor(mutation_type,
                                                   levels=c("Orthosteric site huge differences",
                                                            "Orthosteric site small differences",
                                                            "GTP binding allosteric mutation",
                                                            "GTP binding other mutation",
                                                            "Allosteric mutation",
                                                            "Other mutation")))
  ggplot2::ggplot()+
    ggplot2::geom_rect(data=rect_input,ggplot2::aes(ymin=-1.5,ymax=3.5,xmin=xstart-0.5,xmax=xend+0.5),fill="black",alpha=0.08)+
    ggplot2::geom_rect(data=rect_alpha,ggplot2::aes(ymin=-1.5,ymax=3.5,xmin=xstart-0.5,xmax=xend+0.5),fill="black",alpha=0.05)+
    ggplot2::geom_point(data=data_plot_mutation,ggplot2::aes(x=Pos,y=`mean_kcal/mol`,color=mutation_type),size=0.1)+
      ggplot2::scale_color_manual(values=c(colour_scheme[["red"]],
                                  colour_scheme[["pink"]],
                                  colour_scheme[["blue"]],
                                  colour_scheme[["light blue"]],
                                  colour_scheme[["orange"]],
                                  "gray"),
                         labels=c(paste("Ortho site huge dif",
                                        "=",nrow(data_plot_mutation[mutation_type=="Orthosteric site huge differences",]),
                                        "pos=",length(unique(data_plot_mutation[mutation_type=="Orthosteric site huge differences",Pos]))),
                                  paste("Ortho site small dif",
                                        "=",nrow(data_plot_mutation[mutation_type=="Orthosteric site small differences",]),
                                        "pos=",length(unique(data_plot_mutation[mutation_type=="Orthosteric site small differences",Pos]))),
                                  paste("GTP binding allosteric",
                                        "=",nrow(data_plot_mutation[mutation_type=="GTP binding allosteric mutation",]),
                                        "pos=",length(unique(data_plot_mutation[mutation_type=="GTP binding allosteric mutation",Pos]))),
                                  paste("GTP binding other",
                                        "=",nrow(data_plot_mutation[mutation_type=="GTP binding other mutation",]),
                                        "pos=",length(unique(data_plot_mutation[mutation_type=="GTP binding other mutation",Pos]))),
                                  paste("Allosteric",
                                        "=",nrow(data_plot_mutation[mutation_type=="Allosteric mutation",]),
                                        "pos=",length(unique(data_plot_mutation[mutation_type=="Allosteric mutation",Pos]))),
                                  paste("Other",
                                        "=",length(unique(data_plot_mutation[mutation_type=="Other mutation",Pos])),
                                        "muts = ",nrow(data_plot_mutation[mutation_type=="Other mutation",]))
                         ))+
      ggplot2::geom_hline(yintercept = 0,linetype =2)+
      ggplot2::scale_x_continuous(expand=c(1/188,11/188))+
      ggplot2::ylab("Binding G(kcal/mol)")+
      ggplot2::xlab("Amino acid position")+
      ggplot2::labs(color=NULL)+
      ggplot2::theme_classic2()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size =7, vjust=.5, hjust=.5),
            axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5),
            text = ggplot2::element_text(size=7),
            legend.position="top",
            #legend.key.height= unit(3.1, 'mm'),
            #legend.key.width = unit(5, 'mm'),
            legend.text = ggplot2::element_text(size=7),
            strip.background = ggplot2::element_rect(colour="black", fill="white"))+
      ggplot2::coord_fixed(ratio=10)
  
}