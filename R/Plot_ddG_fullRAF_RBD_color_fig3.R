#' Plot binding free energy changes correlation between full RAF1 and RAF1RBD bindingPCA
#' 
#' This function allows you to plot binding free energy changes correlation between full RAF1 and RAF1RBD bindingPCA.
#' @param ddG1 mochi output
#' @param ddG2 mochi output
#' @param ddG_RAFRBD free energy data ddG_RAFRBD from six binding partner model
#' @param anno structure information annotation data
#' @param assay_name assay:"RAF","PI3"...
#' @return Nothing
#' @export
Plot_ddG_fullRAF_RBD_color_fig3<-function(ddG1=ddG1,
                               ddG2=ddG2,
                               ddG_RAFRBD=ddG_RAFRBD,
                               anno=anno,
                               assay_sele=assay_sele){
  # merge data
  ddG1=Read_ddG(ddG=ddG1,assay_sele="a")
  ddG2=Read_ddG(ddG=ddG2,assay_sele="b")
  all_ddG<-rbind(ddG1,ddG2)
  all_ddG_dc<-dcast(all_ddG[!is.na(mt),],mt+Pos_real~assay,value.var=c("mean_kcal/mol","std_kcal/mol","ci95_kcal/mol"))
  # structure information
  ddG_RAFRBD1<-Read_ddG(ddG = ddG_RAFRBD,
                 assay_sele = assay_sele)
  # ddG_RAFRBD
  weighted_mean_ddG1<-Get_weighted_mean_abs_ddG(ddG=ddG_RAFRBD,assay_sele = "RAF")
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
  ddG_RAFRBD1[,Pos:=Pos_real]
  data_plot_mutation1<-merge(ddG_RAFRBD1,anno,by="Pos",all=T)
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
  #
  all_ddG_dc_anno<-merge(all_ddG_dc,data_plot_mutation,by=c("mt"))
  all_ddG_dc_anno<-all_ddG_dc_anno[!is.na(`mean_kcal/mol_a`)&!is.na(`mean_kcal/mol_b`),]
  ddG_lm<-lm(`mean_kcal/mol_a`~`mean_kcal/mol_b`,all_ddG_dc_anno)
  ddG_lm_allosteric<-lm(`mean_kcal/mol_a`~`mean_kcal/mol_b`,all_ddG_dc_anno[mutation_type=="Allosteric mutation",])
  #plot
  ggplot()+
    geom_point(data=all_ddG_dc_anno,aes(`mean_kcal/mol_a`,`mean_kcal/mol_b`,color=mutation_type),size=0.1)+
    annotate("text",x=0,y=2,,label = paste0("r = ",round(sqrt(summary(ddG_lm)$r.squared),2)) ,size=7*0.35)+
    annotate("text",x=0,y=1.5,,label = paste0("r = ",round(sqrt(summary(ddG_lm_allosteric)$r.squared),2)) ,size=7*0.35,color=colour_scheme[["orange"]])+
    geom_smooth(data=all_ddG_dc_anno,aes(`mean_kcal/mol_a`,`mean_kcal/mol_b`),method=lm, se=T,size=0.1,color="black")+
    geom_smooth(data=all_ddG_dc_anno[mutation_type=="Allosteric mutation",],aes(`mean_kcal/mol_a`,`mean_kcal/mol_b`),method=lm, se=T,size=0.1,color=colour_scheme[["orange"]])+
    theme_classic2()+
    scale_color_manual(values=c(colour_scheme[["red"]],
                                colour_scheme[["pink"]],
                                colour_scheme[["blue"]],
                                colour_scheme[["light blue"]],
                                colour_scheme[["orange"]],
                                "gray"),
                       labels=c(paste("Ortho site huge dif",
                                      "=",nrow(all_ddG_dc_anno[mutation_type=="Orthosteric site huge differences",]),
                                      "pos=",length(unique(all_ddG_dc_anno[mutation_type=="Orthosteric site huge differences",Pos]))),
                                paste("Ortho site small dif",
                                      "=",nrow(all_ddG_dc_anno[mutation_type=="Orthosteric site small differences",]),
                                      "pos=",length(unique(all_ddG_dc_anno[mutation_type=="Orthosteric site small differences",Pos]))),
                                paste("GTP binding allosteric",
                                      "=",nrow(all_ddG_dc_anno[mutation_type=="GTP binding allosteric mutation",]),
                                      "pos=",length(unique(all_ddG_dc_anno[mutation_type=="GTP binding allosteric mutation",Pos]))),
                                paste("GTP binding other",
                                      "=",nrow(all_ddG_dc_anno[mutation_type=="GTP binding other mutation",]),
                                      "pos=",length(unique(all_ddG_dc_anno[mutation_type=="GTP binding other mutation",Pos]))),
                                paste("Allosteric",
                                      "=",nrow(all_ddG_dc_anno[mutation_type=="Allosteric mutation",]),
                                      "pos=",length(unique(all_ddG_dc_anno[mutation_type=="Allosteric mutation",Pos]))),
                                paste("Other",
                                      "=",length(unique(all_ddG_dc_anno[mutation_type=="Other mutation",Pos])),
                                      "muts = ",nrow(all_ddG_dc_anno[mutation_type=="Other mutation",]))
                       ))+
    labs(color=NULL)+
    theme_classic2()+
    xlab("Binding ddG of RAF1")+
    ylab("Binding ddG of RAF1RBD")+
    theme(text = element_text(size = 5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5))+
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}

