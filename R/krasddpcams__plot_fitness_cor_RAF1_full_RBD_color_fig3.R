
#' Plot fitness correlation between full RAF1 and RAF1RBD bindingPCA
#' 
#' This function allows you to plot fitness correlation between full RAF1 and RAF1RBD bindingPCA.
#' @param fitness1 dimsum output
#' @param fitness2 dimsum output
#' @param ddG_RAFRBD free energy data ddG_RAFRBD from six binding partner model
#' @param anno structure information annotation data
#' @param assay_name assay:"RAF","PI3"...
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_fitness_cor_RAF1_full_RBD_color_fig3<-function(
  fitness1=fitness1,
  fitness2=fitness2,
  ddG_RAFRBD=ddG_RAFRBD,
  anno=anno,
  assay_sele=assay_sele,
  colour_scheme
  ){
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
  #return(data_plot_mutation)
  # merge data
  mix_1<-merge(fitness1[Nham_aa==1&STOP!=TRUE,],fitness2[Nham_aa==1&STOP!=TRUE&block=="block1",],by = c("mt1"))
  mix_anno<-merge(mix_1,data_plot_mutation,by.x = c("mt1"), by.y = c("mt"))
  mix_anno<-mix_anno[STOP.x!=TRUE&!is.na(nor_fitness.y)&!is.na(nor_fitness.x),]
  r_model<-lm(nor_fitness.y~nor_fitness.x,mix_anno)
  ggplot2::ggplot(data=mix_anno,ggplot2::aes(nor_fitness.x,nor_fitness.y,color=mutation_type))+
    ggplot2::geom_point(size=0.1)+
    ggplot2::annotate("text",x=-1,y=0,,label = paste0("r = ",round(sqrt(summary(r_model)$r.squared),2)) ,size=7*0.35)+
    ggplot2::geom_smooth(method=lm, se=T,size=0.1,color=colour_scheme[["red"]])+
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
    ggplot2::labs(color=NULL)+
    ggpubr::theme_classic2()+
    ggplot2::xlab("fitness of RAF1")+
    ggplot2::ylab("fitness of RAF1RBD")+
    ggplot2::theme(text = ggplot2::element_text(size = 5),
          strip.text.x = ggplot2::element_text(size = 5),
          strip.text.y = ggplot2::element_text(size = 5),
          legend.text = ggplot2::element_text(size = 5))+
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
  # # structure information
  # RAF_CRD<-krasddpcams__minimum_interchain_distances_from_PDB(
  #   "./DATA/pdb6xgu.ent",
  #   chain_query = "A",
  #   chain_target = "B"
  # )
  # RAF_RBD<-krasddpcams__minimum_interchain_distances_from_PDB(
  #   "/Users/cweng/6vjj.pdb",
  #   chain_query = "A",
  #   chain_target = "B"
  # )
  # mix_1[,bind:="others"]
  # mix_1[AA_Pos1.x%in%RAF_CRD[scHAmin_ligand<5,Pos],bind:="CRD"]
  # mix_1[AA_Pos1.x%in%RAF_RBD[scHAmin_ligand<5,Pos],bind:="RBD+CRD"]
  #plot
  #return(mix_anno)
  # ggplot()+
  #   geom_point(data=mix_anno[STOP.x!=TRUE&mutation_type=="Others",],
  #              aes(x=nor_fitness.x,y=nor_fitness.y),color="gray",alpha=0.6,size=0.1)+
  #   geom_point(data=mix_anno[STOP.x!=TRUE&mutation_type=="Allosteric",],
  #              aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["orange"]],alpha=0.6,size=0.1)+
  #   geom_point(data=mix_anno[STOP.x!=TRUE&mutation_type=="GTP binding other",],
  #              aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["light blue"]],alpha=0.6,size=0.1)+
  #   geom_point(data=mix_anno[STOP.x!=TRUE&mutation_type=="GTP binding allosteric",],
  #              aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["blue"]],alpha=0.6,size=0.1)+
  #   geom_point(data=mix_anno[STOP.x!=TRUE&mutation_type=="Ortho site small dif",],
  #              aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["pink"]],alpha=0.6,size=0.1)+
  #   geom_point(data=mix_anno[STOP.x!=TRUE&mutation_type=="Ortho site huge dif",],
  #              aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["red"]],alpha=0.6,size=0.1)+
  #   theme_classic2()+
  #   labs(color=NULL)+
  #   theme_classic2()+
  #   xlab("fitness of RAF1")+
  #   ylab("fitness of RAF1RBD")+
  #   theme(text = element_text(size = 5),
  #         strip.text.x = element_text(size = 5),
  #         strip.text.y = element_text(size = 5),
  #         legend.text = element_text(size = 5))+
  #   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}

