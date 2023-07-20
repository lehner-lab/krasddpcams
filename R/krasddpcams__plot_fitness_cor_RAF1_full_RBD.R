
#' Plot fitness correlation between full RAF1 and RAF1RBD bindingPCA
#' 
#' This function allows you to plot fitness correlation between full RAF1 and RAF1RBD bindingPCA.
#' @param fitness1 dimsum output
#' @param fitness2 dimsum output
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_fitness_cor_RAF1_full_RBD<-function(
  fitness1=fitness1,
  fitness2=fitness2,
  colour_scheme
  ){
  # merge data
  mix_1<-merge(fullRAF_nor_df_pos[Nham_aa==1,],RAF_nor_df_pos[Nham_aa==1,],by = c("mt1"))
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
  mix_1[,bind:="others"]
  mix_1[AA_Pos1.x%in%RAF_CRD[scHAmin_ligand<5,Pos],bind:="CRD"]
  mix_1[AA_Pos1.x%in%RAF_RBD[scHAmin_ligand<5,Pos],bind:="RBD+CRD"]
  #plot
  ggplot2::ggplot()+
    ggplot2::geom_point(data=mix_1[STOP.x!=TRUE&bind=="others",],
               ggplot2::aes(x=nor_fitness.x,y=nor_fitness.y),color="black",alpha=0.6,size=0.1)+
    ggplot2::geom_point(data=mix_1[STOP.x!=TRUE&bind=="CRD",],
               ggplot2::aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["yellow"]],alpha=0.6,size=0.1)+
    ggplot2::geom_point(data=mix_1[STOP.x!=TRUE&bind=="RBD+CRD",],
               ggplot2::aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["red"]],alpha=0.6,size=0.1)+
    ggpubr::theme_classic2()+
    ggplot2::labs(color=NULL)+
    ggpubr::theme_classic2()+
    ggplot2::xlab("fitness of RAF1")+
    ggplot2::ylab("fitness of RAF1RBD")+
    ggplot2::theme(text = ggplot2::element_text(size = 5),
          strip.text.x = ggplot2::element_text(size = 5),
          strip.text.y = ggplot2::element_text(size = 5),
          legend.text = ggplot2::element_text(size = 5))+
    ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}

