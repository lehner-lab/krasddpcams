#' Scatter plot fitness of binding vs abundance
#' 
#' This function allows you to scatter plot binding fitness vs abundance fitness.
#' @param input all_data_pos
#' @param assay_sele assay_sele
#' @param anno anno_final3
#' 
#' @return Nothing
#' @export
Scatterplot_fitness_bs<-function(input=input,
                                      assay_sele=assay_sele,
                                      anno=anno){
  input_abundance<-input[assay=="stab",]
  input_abundance_single <- Nor_overlap_single_mt_fitness(input_abundance)
  input_binding<-input[assay==assay_sele,]
  input_binding_single <- Nor_overlap_single_mt_fitness(input_binding)
  input_long<-rbind(input_abundance_single,input_binding_single)
  
  input_dc <- dcast(input_long, nt_seq+aa_seq+Nham_aa+AA_Pos1+wtcodon1~assay, value.var=c("nor_fitness_nooverlap", "nor_fitness_nooverlap_sigma","nor_gr_nooverlap","nor_gr_nooverlap_sigma"),drop=TRUE)
  
  input_single_pos <- input_dc
  input_single_pos[,position:=AA_Pos1]
  input_single_pos[,WT_AA:=wtcodon1]
  anno_single<-merge(input_single_pos,anno,by.x=c("position","WT_AA"),by.y =c("Pos","WT_AA"),all=T )
  # anno_single[RAF_scHAmin_ligand<=5,RAF_type:="binding interface"]
  # anno_single[RAF_scHAmin_ligand>5&SASA<0.25,RAF_type:="core"]
  # anno_single[RAF_scHAmin_ligand>5&SASA>=0.25,RAF_type:="surface"]
  # anno_single[as.numeric(position)>166,RAF_type:="HVR"]
  anno_single[,RAF_type_bs:="others"]
  anno_single[get(paste0("scHAmin_ligand_",assay_sele))<5,RAF_type_bs:="binding interface"]
  ggplot()+
    geom_point(data=anno_single[position>1&RAF_type_bs=="others",],
               aes(x=nor_fitness_nooverlap_stab,y=get(paste0("nor_fitness_nooverlap_",assay_sele))),color="black",alpha=0.6,size=0.1)+
    geom_point(data=anno_single[position>1&RAF_type_bs=="binding interface"],
               aes(x=nor_fitness_nooverlap_stab,y=get(paste0("nor_fitness_nooverlap_",assay_sele))),color=colour_scheme[["red"]],alpha=0.6,size=0.5)+
    theme_classic2()+
    labs(color=NULL)+
    xlab("Abundance")+
    ylab("RAF1 Binding")+
    theme(text = element_text(size=7),
          legend.position="right",
          legend.text = element_text(size=7),
          axis.text.x = element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))+
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
}