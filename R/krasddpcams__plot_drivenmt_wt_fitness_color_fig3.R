
#' Plot fitness correlation between wt background and mt background
#' 
#' This function allows you to plot fitness correlation between wt background and mt background
#' @param input dimsum output
#' @param driven_mt driven mutation
#' @param ddG_RAFRBD free energy data ddG_RAFRBD from six binding partner model
#' @param anno structure information annotation data
#' @param assay_name assay:"RAF","PI3"...
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_drivenmt_wt_fitness_color_fig3<-function(
  input = input,
  driven_mt = driven_mt,
  ddG_RAFRBD=ddG_RAFRBD,
  anno=anno,
  assay_sele=assay_sele,
  colour_scheme
  ){
  # annotation
  ddG1<-krasddpcams__read_ddG(ddG = ddG_RAFRBD,
                 assay_sele = assay_sele)
  weighted_mean_ddG1<-krasddpcams__get_weighted_mean_abs_ddG(ddG=ddG_RAFRBD,assay_sele = assay_sele)
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
  data_plot_mutation[,allosteric_mutation := p.adjust(krasddpcams__pvalue(abs(mean)-reg_threshold, std), method = "BH")<0.05 & (abs(mean)-reg_threshold)>0]
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
  # data loading
  input_WT<-input[Nham_aa<=1,]
  input_WT[,mt_sec:=mt1]
  input_WT[Nham_aa==0,mt_sec:="base"]
  input_MT<-input[Nham_aa>=1&(mt1==driven_mt|mt2==driven_mt),]
  input_MT[mt2==driven_mt,mt_sec:=mt1]
  input_MT[mt1==driven_mt,mt_sec:=mt2]
  input_MT[mt1==driven_mt&Nham_aa==1,mt_sec:="base"]
  input_MT[mt2==driven_mt&Nham_aa==1,mt_sec:="base"]
  input_melt<-input_WT[input_MT,on=.(mt_sec),nomatch = NULL]
  mix_anno<-merge(input_melt,data_plot_mutation,by.x = c("mt_sec"), by.y = c("mt"))
  mix_anno<-mix_anno[STOP!=TRUE&i.STOP!=TRUE&!is.na(nor_fitness)&!is.na(i.nor_fitness),]
  
  # correlation
  r_model<-lm(nor_fitness~i.nor_fitness,mix_anno)
  r_model_allosteric<-lm(nor_fitness~i.nor_fitness,mix_anno[mutation_type=="Allosteric mutation",])
  # plot
  ggplot2::ggplot(mix_anno[AA_Pos1<65,],ggplot2::aes(x=nor_fitness,y=i.nor_fitness,color=mutation_type))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::annotate("text",x=-0.5,y=0.2,label = paste0("r = ",round(sqrt(summary(r_model)$r.squared),2)) ,size=7*0.35,color="black")+
    ggplot2::geom_smooth(method=lm, se=T,size=0.1,color="black")+
    ggplot2::annotate("text",x=-0.5,y=0.1,label = paste0("r = ",round(sqrt(summary(r_model_allosteric)$r.squared),2)) ,size=7*0.35,color=colour_scheme[["orange"]])+
    ggplot2::geom_smooth(method=lm, se=T,size=0.1,color=colour_scheme[["orange"]])+
    ggplot2::xlab("Binding fitness in the background of WT")+
    ggplot2::ylab(paste("Binding fitness in the background of",driven_mt))+
    ggplot2::geom_vline(xintercept = input_WT[WT==TRUE,nor_fitness])+
    ggplot2::geom_hline(yintercept = input_WT[mt1==driven_mt,nor_fitness])+
    ggpubr::theme_classic2()+
    ggplot2::scale_color_manual(values=c(colour_scheme[["red"]],
                                colour_scheme[["pink"]],
                                colour_scheme[["blue"]],
                                colour_scheme[["light blue"]],
                                colour_scheme[["orange"]],
                                "gray"),
                       labels=c(paste("Ortho site huge dif",
                                      "=",nrow(mix_anno[mutation_type=="Orthosteric site huge differences",]),
                                      "pos=",length(unique(mix_anno[mutation_type=="Orthosteric site huge differences",Pos]))),
                                paste("Ortho site small dif",
                                      "=",nrow(mix_anno[mutation_type=="Orthosteric site small differences",]),
                                      "pos=",length(unique(mix_anno[mutation_type=="Orthosteric site small differences",Pos]))),
                                paste("GTP binding allosteric",
                                      "=",nrow(mix_anno[mutation_type=="GTP binding allosteric mutation",]),
                                      "pos=",length(unique(mix_anno[mutation_type=="GTP binding allosteric mutation",Pos]))),
                                paste("GTP binding other",
                                      "=",nrow(mix_anno[mutation_type=="GTP binding other mutation",]),
                                      "pos=",length(unique(mix_anno[mutation_type=="GTP binding other mutation",Pos]))),
                                paste("Allosteric",
                                      "=",nrow(mix_anno[mutation_type=="Allosteric mutation",]),
                                      "pos=",length(unique(mix_anno[mutation_type=="Allosteric mutation",Pos]))),
                                paste("Other",
                                      "=",length(unique(mix_anno[mutation_type=="Other mutation",Pos])),
                                      "muts = ",nrow(mix_anno[mutation_type=="Other mutation",]))
                       ))+
    ggplot2::theme(text = ggplot2::element_text(size = 5),
          strip.text.x = ggplot2::element_text(size = 5),
          strip.text.y = ggplot2::element_text(size = 5),
          legend.text = ggplot2::element_text(size = 5))+
    #scale_x_continuous(limits = c(-1.5, highest_gr)) +
    #scale_y_continuous(limits = c(-1.5, highest_gr)) +
    ggplot2::geom_abline(linetype="dashed")+
    ggplot2::coord_fixed()
}