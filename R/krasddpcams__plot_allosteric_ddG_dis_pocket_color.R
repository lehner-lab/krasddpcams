
#' Plot allosteric mutation vs distance colored by pockets
#' 
#' This function allows you to plot binding interface mutations vs distance.
#' @param ddG free energy data
#' @param anno structure information annotation data
#' @param pocket_anno final_distance_dc_anno_for_anno_mt_4
#' @param assay_sele assay:"RAF","PI3"...
#' @param pocket pocket name
#' @param pocket2 pocket2 name, Sotorasib_dis
#' @param colour_scheme colour scheme list
#' 
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_allosteric_ddG_dis_pocket_color<-function(
  ddG=ddG,
  anno=anno,
  pocket_anno=pocket_anno,
  assay_sele=assay_sele,
  pocket=pocket,
  pocket2=pocket2,
  colour_scheme
  ){
  ddG1<-krasddpcams__read_ddG(ddG = ddG,
                 assay_sele = assay_sele)
  
  weighted_mean_ddG1<-krasddpcams__get_weighted_mean_abs_ddG(ddG=ddG,assay_sele = assay_sele)
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
  # data_plot_mutation[,mutation_type:="Reminder"]
  data_plot_mutation[,allosteric_mutation := p.adjust(krasddpcams__pvalue(abs(mean)-reg_threshold, std), method = "BH")<0.05 & (abs(mean)-reg_threshold)>0]
  # data_plot_mutation[,colors_type:="others"]
  # data_plot_mutation[mutation_type=="GTP binding site",colors_type:="GTP binding site"]
  if(pocket=="pocket1"){
    data_plot_mutation[Pos%in%pocket_anno[contact==pocket,Pos],colors_type:="pocket"]
  }else if(pocket=="pocket2"){
    data_plot_mutation[Pos%in%pocket2[drug_scHAmin_ligand<5,Pos],colors_type:="pocket"]
  }else if(pocket=="pocket3"){
    data_plot_mutation[Pos%in%pocket_anno[contact==pocket,Pos],colors_type:="pocket"]
  }else if(pocket=="pocket4"){
    data_plot_mutation[Pos%in%pocket_anno[contact==pocket,Pos],colors_type:="pocket"]
  }
  data_plot_mutation[,colors_type2:="Other"]
  # data_plot_mutation[Pos%in%data_plot[get(paste0("scHAmin_ligand_",assay_sele))<5,Pos]&colors_type=="pocket",colors_type2:="Binding interface mutation & pocket"]
  #data_plot_mutation[!Pos%in%data_plot[get(paste0("scHAmin_ligand_",assay_sele))<5,Pos]&colors_type!="pocket"&allosteric_mutation==T,colors_type2:="Other allosteric mutation"]
  data_plot_mutation[colors_type=="pocket"&allosteric_mutation==T,colors_type2:="Pocket allosteric mutation"]
  #data_plot_mutation[!Pos%in%data_plot[get(paste0("scHAmin_ligand_",assay_sele))<5,Pos]&colors_type!="pocket"&allosteric_mutation==F,colors_type2:="Other mutation"]
  data_plot_mutation[colors_type=="pocket"&allosteric_mutation==F,colors_type2:="Other pocket mutation"]
  data_plot_mutation<-within(data_plot_mutation,
                             factor(colors_type2,
                                    levels = c("Pocket allosteric mutation","Other pocket mutation","Other")))
  ggplot2::ggplot(data_plot_mutation,ggplot2::aes(x=scHAmin_ligand_RAF,y=`mean_kcal/mol`))+
    ggplot2::geom_point(ggplot2::aes(color=colors_type2),size=0.1)+
    ggplot2::geom_vline(xintercept = 5,linetype =2)+
    ggplot2::scale_color_manual(values = c("gray",colour_scheme[["pink"]],colour_scheme[["red"]]),
                       labels = c(paste("Pocket allosteric mutation","=",
                                        nrow(data_plot_mutation[colors_type2=="Pocket allosteric mutation",]),
                                        "site=",
                                        length(unique(data_plot_mutation[colors_type2=="Pocket allosteric mutation",Pos_real]))),
                                  paste("Other pocket mutation","=",
                                        nrow(data_plot_mutation[colors_type2=="Other pocket mutation",])),
                                  paste("Other","=",
                                        nrow(data_plot_mutation[colors_type2=="Other",]))
                                  )
                       
                        )+
  ggplot2::ylab(paste("Binding free energy change \n",assay_sele,"(kcal/mol)"))+
  ggplot2::xlab(expression(paste("Distance to binding partner ("*ring(A)*")")))+
  ggplot2::theme_classic()+
  ggplot2::theme(text = ggplot2::element_text(size = 7),
        axis.text = ggplot2::element_text(size = 7),legend.text = ggplot2::element_text(size = 7))+
  ggplot2::ggtitle(pocket)+
  ggplot2::coord_fixed(ratio = 5)
    

  
  # data_pl
}