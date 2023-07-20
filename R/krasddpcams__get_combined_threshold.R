
#' Plot major allosteric sites
#' 
#' This function allows you to plot major allosteric sites.
#' @param ddG1 free energy data1
#' @param assay1 data1's name
#' @param ddG2 free energy data2
#' @param assay2 data2's name
#' @param ddG3 free energy data3
#' @param assay3 data3's name
#' @param ddG4 free energy data4
#' @param assay4 data4's name
#' @param ddG5 free energy data5
#' @param assay5 data5's name
#' @param ddG6 free energy data6
#' @param assay6 data6's name
#' @param anno anno_final3
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__get_combined_threshold<-function(
  ddG1=ddG1,
  assay1=assay1,
  ddG2=ddG2,
  assay2=assay2,
  ddG3=ddG3,
  assay3=assay3,
  ddG4=ddG4,
  assay4=assay4,
  ddG5=ddG5,
  assay5=assay5,
  ddG6=ddG6,
  assay6=assay6,
  anno=anno
  ){
  weighted_mean_ddG<-list()
  weighted_mean_ddG[[assay1]]<-krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG=ddG1,
                                                              assay_sele = assay1)
  weighted_mean_ddG[[assay2]]<-krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG=ddG2,
                                                              assay_sele = assay2)
  weighted_mean_ddG[[assay3]]<-krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG=ddG3,
                                                              assay_sele = assay3)
  weighted_mean_ddG[[assay4]]<-krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG=ddG4,
                                                              assay_sele = assay4)
  weighted_mean_ddG[[assay5]]<-krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG=ddG5,
                                                              assay_sele = assay5)
  weighted_mean_ddG[[assay6]]<-krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG=ddG6,
                                                              assay_sele = assay6)
  
  assay_list<-list(assay1,assay2,assay3,assay4,assay5,assay6)
  
  data_plot<-data.table()
  for(assayi in assay_list){
    data_plot_assayi<-merge(weighted_mean_ddG[[assayi]],anno,by="Pos",all=T)
    data_plot_assayi[,binding_type:="allosteric site"]
    data_plot_assayi[get(paste0("scHAmin_ligand_",assayi))<5,binding_type:="binding site"]
    data_plot_assayi[,binding_type_gtp_included:=binding_type]
    data_plot_assayi[get(paste0("GXPMG_scHAmin_ligand_",assayi))<5,binding_type_gtp_included:="GTP binding site"]
    data_plot<-rbind(data_plot,data_plot_assayi)
  }
  reg_threshold<-data_plot[binding_type=="binding site",sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean","sigma")]
  reg_threshold
  
  
}