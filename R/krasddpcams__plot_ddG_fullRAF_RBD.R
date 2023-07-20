
#' Plot binding free energy changes correlation between full RAF1 and RAF1RBD bindingPCA
#' 
#' This function allows you to plot binding free energy changes correlation between full RAF1 and RAF1RBD bindingPCA.
#' @param ddG1 mochi output
#' @param ddG2 mochi output
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_ddG_fullRAF_RBD<-function(
  ddG1=ddG1,
  ddG2=ddG2,
  ci95_1=ci95_1,
  ci95_2=ci95_2,
  colour_scheme
  ){
  # merge data
  ddG1=krasddpcams__read_ddG(ddG=ddG1,assay_sele="a")
  ddG2=krasddpcams__read_ddG(ddG=ddG2,assay_sele="b")
  all_ddG<-rbind(ddG1,ddG2)
  all_ddG_dc<-dcast(all_ddG[!is.na(mt),],mt+Pos_real~assay,value.var=c("mean_kcal/mol","std_kcal/mol","ci95_kcal/mol"))
  # structure information
  RAF_CRD<-krasddpcams__minimum_interchain_distances_from_PDB(
    "./DATA/pdb6xgu.ent",
    chain_query = "A",
    chain_target = "B"
  )
  RAF_RBD<-krasddpcams__minimum_interchain_distances_from_PDB(
    "./DATA/6vjj.pdb",
    chain_query = "A",
    chain_target = "B"
  )
  all_ddG_dc[,bind:="others"]
  all_ddG_dc[Pos_real%in%RAF_CRD[scHAmin_ligand<5,Pos],bind:="CRD"]
  all_ddG_dc[Pos_real%in%RAF_RBD[scHAmin_ligand<5,Pos],bind:="RBD+CRD"]
  # bind_lm_ddG<-lm(`mean_kcal/mol_a`~`mean_kcal/mol_b`,all_ddG_dc[`std_kcal/mol_a`<0.1&`std_kcal/mol_b`<0.1,])
  # #plot
  # ggplot()+
  #   geom_point(data=all_ddG_dc[bind=="others"&`std_kcal/mol_a`<0.1&`std_kcal/mol_b`<0.1,],
  #              ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),color="black",alpha=0.6,size=0.1)+
  #   geom_point(data=all_ddG_dc[bind=="CRD"&`std_kcal/mol_a`<0.1&`std_kcal/mol_b`<0.1,],
  #              ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),color=colour_scheme[["yellow"]],alpha=0.6,size=0.1)+
  #   geom_point(data=all_ddG_dc[bind=="RBD+CRD"&`std_kcal/mol_a`<0.1&`std_kcal/mol_b`<0.1,],
  #              ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),color=colour_scheme[["red"]],alpha=0.6,size=0.1)+
  #   geom_smooth(data=all_ddG_dc[`std_kcal/mol_a`<0.1&`std_kcal/mol_b`<0.1,],ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),
  #               method=lm, se=T,size=0.1,color=colour_scheme[["red"]])+
  #   annotate("text",x=-0.5,y=2.00,label = paste0("r = ",round(sqrt(summary(bind_lm_ddG)$r.squared),2)) ,size=7*0.35)+
  #   theme_classic2()+
  #   labs(color=NULL)+
  #   geom_hline(yintercept = 0) +
  #   geom_vline(xintercept = 0) +
  #   theme_classic()+
  #   geom_abline(linetype="dashed")+
  #   xlab("Binding ddG of RAF1")+
  #   ylab("Binding ddG of RAF1RBD")+
  #   theme(text = element_text(size = 5),
  #         strip.text.x = element_text(size = 5),
  #         strip.text.y = element_text(size = 5),
  #         legend.text = element_text(size = 5))+
  #   coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
  bind_lm_ddG<-lm(`mean_kcal/mol_a`~`mean_kcal/mol_b`,all_ddG_dc[`ci95_kcal/mol_a`<ci95_1&`ci95_kcal/mol_b`<ci95_2,])
  #plot
  ggplot2::ggplot()+
    ggplot2::geom_point(data=all_ddG_dc[bind=="others"&`ci95_kcal/mol_a`<ci95_1&`ci95_kcal/mol_b`<ci95_2,],
               ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),color="black",alpha=0.6,size=2)+
    ggplot2::geom_point(data=all_ddG_dc[bind=="CRD"&`ci95_kcal/mol_a`<ci95_1&`ci95_kcal/mol_b`<ci95_2,],
               ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),color=colour_scheme[["yellow"]],alpha=0.6,size=2)+
    ggplot2::geom_point(data=all_ddG_dc[bind=="RBD+CRD"&`ci95_kcal/mol_a`<ci95_1&`ci95_kcal/mol_b`<ci95_2,],
               ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),color=colour_scheme[["red"]],alpha=0.6,size=2)+
    ggplot2::geom_smooth(data=all_ddG_dc[`ci95_kcal/mol_a`<ci95_1&`ci95_kcal/mol_b`<ci95_2,],ggplot2::aes(x=`mean_kcal/mol_a`,y=`mean_kcal/mol_b`),
                method=lm, se=T,size=0.1,color=colour_scheme[["red"]])+
    ggplot2::annotate("text",x=-0.5,y=2.00,label = paste0("r = ",round(sqrt(summary(bind_lm_ddG)$r.squared),2)) ,size=7*0.35)+
    ggpubr::theme_classic2()+
    ggplot2::labs(color=NULL)+
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()+
    ggplot2::geom_abline(linetype="dashed")+
    ggplot2::xlab("Binding ddG of RAF1")+
    ggplot2::ylab("Binding ddG of RAF1RBD")+
    ggplot2::theme(text = ggplot2::element_text(size = 5),
          strip.text.x = ggplot2::element_text(size = 5),
          strip.text.y = ggplot2::element_text(size = 5),
          legend.text = ggplot2::element_text(size = 5))+
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}

