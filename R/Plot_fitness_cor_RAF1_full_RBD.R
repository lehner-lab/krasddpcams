#' Plot fitness correlation between full RAF1 and RAF1RBD bindingPCA
#' 
#' This function allows you to plot fitness correlation between full RAF1 and RAF1RBD bindingPCA.
#' @param fitness1 dimsum output
#' @param fitness2 dimsum output
#' @return Nothing
#' @export
Plot_fitness_cor_RAF1_full_RBD<-function(fitness1=fitness1,
                                         fitness2=fitness2){
  # merge data
  mix_1<-merge(fullRAF_nor_df_pos[Nham_aa==1,],RAF_nor_df_pos[Nham_aa==1,],by = c("mt1"))
  # structure information
  RAF_CRD<-doubledeepms__minimum_interchain_distances_from_PDB(
    "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/pdb6xgu.ent",
    chain_query = "A",
    chain_target = "B"
  )
  RAF_RBD<-doubledeepms__minimum_interchain_distances_from_PDB(
    "/Users/cweng/6vjj.pdb",
    chain_query = "A",
    chain_target = "B"
  )
  mix_1[,bind:="others"]
  mix_1[AA_Pos1.x%in%RAF_CRD[scHAmin_ligand<5,Pos],bind:="CRD"]
  mix_1[AA_Pos1.x%in%RAF_RBD[scHAmin_ligand<5,Pos],bind:="RBD+CRD"]
  #plot
  ggplot()+
    geom_point(data=mix_1[STOP.x!=TRUE&bind=="others",],
               aes(x=nor_fitness.x,y=nor_fitness.y),color="black",alpha=0.6,size=0.1)+
    geom_point(data=mix_1[STOP.x!=TRUE&bind=="CRD",],
               aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["yellow"]],alpha=0.6,size=0.1)+
    geom_point(data=mix_1[STOP.x!=TRUE&bind=="RBD+CRD",],
               aes(x=nor_fitness.x,y=nor_fitness.y),color=colour_scheme[["red"]],alpha=0.6,size=0.1)+
    theme_classic2()+
    labs(color=NULL)+
    theme_classic2()+
    xlab("fitness of RAF1")+
    ylab("fitness of RAF1RBD")+
    theme(text = element_text(size = 5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5))+
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}

