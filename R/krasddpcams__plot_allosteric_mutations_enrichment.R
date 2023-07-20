
#' Plot allosteric mutations enrichment
#' 
#' This function allows you to plot allosteric mutations.
#' @param ddG free energy data
#' @param anno structure information annotation data
#' @param assay_sele assay:"RAF","PI3"...
#' @param threshold threshold
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_allosteric_mutations_enrichment<-function(
  ddG=ddG,
  anno=anno,
  assay_sele=assay_sele,
  threshold=threshold,
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
  # reg_threshold<-data_plot[binding_type=="binding site",sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean","sigma")]
  reg_threshold<-threshold
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
  data_plot_mutation
  result_list <- list()
  bset_list <- list(
    "Charged" = c("R", "H", "D", "E", "K"), 
    "Hydrophobic" = c("A", "V", "I", "L", "M", "F", "Y", "W"))
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  names(aa_list) <- unlist(aa_list)
  bset_list <- c(bset_list, aa_list)
  data_plot_mutation[,WT_AA:=wt_codon]
  data_plot_mutation[,Mut:=mt_codon]
  for(lset_name in names(bset_list)){
    lset <- bset_list[[lset_name]]
    in_lset_allo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][WT_AA %in% lset & allosteric_mutation==T,.N]
    out_lset_allo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][!WT_AA %in% lset & allosteric_mutation==T,.N]
    in_lset_nallo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][WT_AA %in% lset & allosteric_mutation==F,.N]
    out_lset_nallo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][!WT_AA %in% lset & allosteric_mutation==F,.N]
    temp_test <- fisher.test(matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2))
    result_list <- c(result_list, list(data.table(assay = assay_sele, set_name = lset_name, mutation = "WT", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
    in_lset_allo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][Mut %in% lset & allosteric_mutation==T,.N]
    out_lset_allo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][!Mut %in% lset & allosteric_mutation==T,.N]
    in_lset_nallo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][Mut %in% lset & allosteric_mutation==F,.N]
    out_lset_nallo <- data_plot_mutation[!mutation_type%in%c("Orthosteric site huge differences","Orthosteric site small differences")][!Mut %in% lset & allosteric_mutation==F,.N]
    temp_test <- fisher.test(matrix(c(in_lset_allo, out_lset_allo, in_lset_nallo, out_lset_nallo), nrow = 2))
    result_list <- c(result_list, list(data.table(assay = assay_sele, set_name = lset_name, mutation = "Mutant", odds_ratio = temp_test$estimate, p_value = temp_test$p.value)))
  }
  result<-rbindlist(result_list)
  result<-within(result,
                 set_name<-factor(set_name,
                                  levels = c("P","G","D","E","R","H","K","Charged",
                                             "A","V","I","L","M","F","Y","W","Hydrophobic",
                                             "C","S","T","N","Q")))
  ggplot2::ggplot(result,ggplot2::aes(set_name, log2(odds_ratio), fill = mutation, alpha = p_value<0.05)) +
    ggplot2::scale_fill_manual(values=c(colour_scheme[["red"]],
                                         colour_scheme[["blue"]]))+
    ggplot2::geom_col(position = "dodge") +
    ggplot2::xlab("Binding partners") +
    ggplot2::ylab("log2(odds ratio)") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size=7),
                   legend.position="right",
                   legend.text = ggplot2::element_text(size=7),
                   axis.text.x = ggplot2::element_text(size =7, vjust=.5, hjust=.5),
                   axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5,margin=ggplot2::margin(0,0,0,0,"mm")),
                   legend.key.height= ggplot2::unit(3.1, 'mm'),
                   legend.key.width = ggplot2::unit(3.1, 'mm'),
                   legend.key.size = ggplot2::unit(1,"mm"),
                   plot.margin=ggplot2::margin(0,0,0,0))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}