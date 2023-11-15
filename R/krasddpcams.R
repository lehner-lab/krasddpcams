
#' krasddpcams
#'
#' Main analysis script.
#'
#' @param base_dir Base directory for all output files (default:NB private CRG server path; change accordingly)
#'
#' @return Nothing
#' @export
#' @import data.table
krasddpcams <- function(
  base_dir = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_updating_20231115/"
  ){

  #Set working directory and create output folders
  setwd(base_dir)
  dir.create("FIGURES")
  dir.create("Figures_revision")
  dir.create("RESULTS")
  
  ######
  #basic setting
  ######

  options(ggrepel.max.overlaps = Inf)
  colour_scheme<- list(
    "blue"="#1B38A6",#rgb(27, 56, 166)
    "red"="#F4270C",#rgb(244, 39, 12)
    "orange"="#F4AD0C",#rgb(244, 173, 12)
    "green"="#09B636",#rgb(9, 182, 54)
    "yellow"="#F1DD10",#rgb(241, 221, 16)
    "purple"="#6D17A0",#rgb(109, 23, 160)
    "pink"="#FFB0A5",#rgb(255, 176, 165)
    "light orange"="#FFE4A5",#rgb(255, 228, 165)
    "light blue"="#9DACE3",#rgb(157, 172, 227)
    "light green"="#97E9AD",#rgb(151, 233, 173)
    "light red"="#FF6A56",#rgb(255, 106, 86)
    "dark red"="#A31300",#rgb(163, 19, 0)
    "dark blue"="#0C226F",#rgb(12, 34, 111)
    "dark green"="#007A20"#rgb(0, 122, 32)
  )
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))

  ######
  #KRAS information
  ######

  wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
  ###beta strands
  rects<-data.frame(xstart=c(2,37,49,77,111,141),
                    xend=c(10,46,58,83,116,143), 
                    col=c("b1","b2","b3","b4","b5","b6"))
  ###alpha helix 1
  rects_alpha<-data.frame(xstart=c(16),
                          xend=c(25
                          ), 
                          col=c("a1"))
  ###all alpha helix
  rects_alpha_all<-data.frame(xstart=c(16,62,87,127, 152),
                              xend=c(25,74,104, 137,166 ), 
                              col=c("a1","a2","a3","a4","a5"))
  ##
  all_distance_anno<-fread("./DATA/all_distance_anno.csv")
  all_long_distance<-fread("./DATA/all_long_distance.csv")
  all_long_distance_GXPMG<-fread("./DATA/all_long_distance_GXPMG.csv")
  final_distance_dc_anno_for_anno_mt_5<-fread("./DATA/final_distance_dc_anno_for_anno_mt_5.csv")
  anno_final3<-fread("./DATA/anno_final3.csv")

  ######
  #Load data
  ######

  stability_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                   block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                   block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData")
  RAF_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                             block2_dimsum_df = "./DATA/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                             block3_dimsum_df = "./DATA/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")
  PI3_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                             block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                             block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData")
  RAL_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                                             block2_dimsum_df = "./DATA/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                                             block3_dimsum_df = "./DATA/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData")
  SOS_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                                             block2_dimsum_df = "./DATA/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                                             block3_dimsum_df = "./DATA/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData")
  K27_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                                             block2_dimsum_df = "./DATA/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                                             block3_dimsum_df = "./DATA/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData")
  K55_nor_df <- krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                                             block2_dimsum_df = "./DATA/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                                             block3_dimsum_df = "./DATA/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData")

  stab <- stability_nor_df
  RAF <- RAF_nor_df
  K27 <- K27_nor_df
  PI3 <- PI3_nor_df
  RAL <- RAL_nor_df
  SOS <- SOS_nor_df
  K55 <- K55_nor_df
  all_data <- krasddpcams__merge_dimsum_df(stab,RAF,K27,PI3,RAL,SOS,K55)
  all_data_pos<-krasddpcams__pos_id(all_data,wt_aa)
  fullRAF_nor_df <- krasddpcams__normalize_growthrate_fitness_block1(block1_dimsum_df = "./DATA/CW_RAS_binding_FULLRAF_1_fitness_replicates.RData")
  fullRAF_nor_df_pos <- krasddpcams__pos_id(fullRAF_nor_df,wt_aa)
  GAP_nor_df <- krasddpcams__normalize_growthrate_fitness_block1(block1_dimsum_df = "./DATA/CW_RAS_binding_RAFGAP_1_fitness_replicates_fullseq.RData")
  GAP_nor_df_pos <- krasddpcams__pos_id(GAP_nor_df,wt_aa)
  RAF_nor_df_pos <- krasddpcams__pos_id(RAF_nor_df,wt_aa)

  ######
  #Figure 1
  ######

  #1e
  krasddpcams__plot_heatmap_fitness(all_data_pos,"stab",anno_final3,wt_aa, colour_scheme)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(family="Courier",angle=90,size=5, vjust =.5,hjust =.5,margin=ggplot2::margin(0,-0.5,0,0,"mm")))
  ggplot2::ggsave("./FIGURES/figure1e_fitness_heatmap.pdf", device = cairo_pdf,height = 4,width=20)

  #1f
  krasddpcams__plot_heatmap_fitness(all_data_pos,"RAF",anno_final3,wt_aa, colour_scheme)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(family="Courier",angle=90,size=5, vjust =.5,hjust =.5,margin=ggplot2::margin(0,-0.5,0,0,"mm")))
  ggplot2::ggsave("./FIGURES/figure1f_fitness_heatmap.pdf", device = cairo_pdf,height = 4,width=20)

  #1g
  krasddpcams__plot_heatmap_ddG(input="./DATA/weights_Folding.txt",
                        assay_name = "folding", wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure1g_ddg_heatmap.pdf", device = cairo_pdf,height = 4,width=20)

  #1h
  krasddpcams__plot_heatmap_ddG(input="./DATA/weights_Binding_RAF.txt",
                        assay_name = "RAF1", wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure1h_ddg_heatmap.pdf", device = cairo_pdf,height = 4,width=20)

  #1j
  krasddpcams__scatterplot_fitness_bs(input = all_data_pos,
                              assay_sele="RAF",
                              anno = anno_final3, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure1j_fitness_bs.pdf", device = cairo_pdf,height = 4,width=4)

  #1k
  krasddpcams__plot_RAF_invitro_cor(input="./DATA/weights_Binding_RAF.txt",
                            assay_name = "RAF1", colour_scheme)
  ggplot2::ggsave("./FIGURES/figure1k_RAF1_invitro_ddG_cor.pdf", device = cairo_pdf,height = 50,width=50,units = "mm")

  #1l
  krasddpcams__weighted_mean_ddG_structure(ddG="./DATA/weights_Folding.txt",
                              input_PDB = "./DATA/pdb6vjj.ent",
                              chain_KRAS = "A",
                              output_PDB_file="./FIGURES/figure1l_folding_ddG_pdb6vjj.pdb")
  #to make figure 1h using chimeraX

  ######
  #Figure 2
  ######

  #2A
  krasddpcams__plot_scatter_ddGb_ddGf(ddG1="./DATA/weights_Folding.txt",
                              assay1="folding",
                              ddG2="./DATA/weights_Binding_RAF.txt",
                              assay2="RAF1",
                              anno=anno_final3, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure2a_scatter_ddGb_ddGf.pdf", device = cairo_pdf,height = 4,width=4)

  #2B
  weighted_mean_abs_ddG_fold<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Folding.txt",
                                                        assay_sele = "folding")
  weighted_mean_abs_ddG_RAF<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_RAF.txt",
                                                       assay_sele = "RAF")
  krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                 weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                 weighted_meab_abs_ddGb=weighted_mean_abs_ddG_RAF,
                                                 anno_input=anno_final3,
                                                 bind="RAF",
                                                 bind_ligand = "scHAmin_ligand_RAF")+ggplot2::ggtitle("RAF1")
  ggplot2::ggsave("./FIGURES/figure2b_ROC.pdf", device = cairo_pdf,height = 4,width=4)

  #2C
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                                   input_PDB = "./DATA/pdb6vjj.ent",
                                   chain_KRAS = "A",
                                   output_PDB_file="./FIGURES/figure2c_RAF_binding_ab_ddG_6vjj.pdb")

  #to make figure 2C using chimeraX, key command: color bfactor /A range -2,2
  #2D
  krasddpcams__plot_binding_interface_mutation_heatmap("./DATA/weights_Binding_RAF.txt",
                                               assay_sele = "RAF1", wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure2d_binding_interface_mutation_heatmap.pdf", device = cairo_pdf,height = 10,width=4)

  ######
  #Figure 3
  ######

  #3A
  krasddpcams__plot_allosteric_mutations(ddG="./DATA/weights_Binding_RAF.txt",
                            anno=anno_final3,
                            assay_sele="RAF",
                            rect_input=rects,
                            rect_alpha=rects_alpha, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure3a_allosteric_mutations.pdf", device = cairo_pdf,height = 4,width=10)

  #3B
  RAF_list<-krasddpcams__plot_major_allosteric_site(ddG="./DATA/weights_Binding_RAF.txt",
                                       anno=anno_final3,
                                       assay_name="RAF",
                                       allosteric_list = RAF_list,
                                       rect_input=rects, colour_scheme)
  RAF_list[["p"]]
  ggplot2::ggsave("./FIGURES/figure3b_major_allosteric_site.pdf", device = cairo_pdf,height = 4,width=6)

  #3c
  #pdb 6vjj, chimeraX,
  #paste0(RAF_list[["Binding interface site"]],collapse=",")
  #paste0(RAF_list[["Allosteric GTP pocket site"]],collapse=",")
  #paste0(RAF_list[["Major allosteric site"]],collapse=",")

  #3d
  krasddpcams__maximum_ddG_structure(ddG="./DATA/weights_Binding_RAF.txt",
                        input_PDB = "./DATA/pdb6vjj.ent",
                        chain_KRAS = "A",
                        output_PDB_file="./FIGURES/figure3d_max_RAF_binding_ddG_6vjj.pdb")
  #to make figure3d using chimeraX, key: color bfactor /A range -2,2

  #3e
  krasddpcams__plot_ddG_betasheet(ddG="./DATA/weights_Binding_RAF.txt",
                     weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                     anno=anno_final3,
                     assay_sele="RAF",
                     rect_input=rects)
  ggplot2::ggsave("./FIGURES/figure3e_ddG_betasheet.pdf", device = cairo_pdf,height = 4,width=4)

  #3f
  #6vjj and 6oim, using chimeraX, key: merge #1/A to #2/B

  #3g
  krasddpcams__maximum_ddG_structure(ddG="./DATA/weights_Binding_RAF.txt",
                        input_PDB = "./DATA/pdb6oim.ent",
                        chain_KRAS = "A",
                        output_PDB_file="./FIGURES/figure3g_max_RAF_binding_ddG_6oim.pdb")
  #to make figure3g using chimeraX, key:color bfactor sel range -2,2

  ######
  #Figure 4
  ######

  #4a
  krasddpcams__plot_figure4A_point_all_ddG_alphabeta(ddG1="./DATA/weights_Folding.txt",
                                             assay1="folding",
                                             ddG2="./DATA/weights_Binding_RAF.txt",
                                             assay2="RAF1",
                                             ddG3="./DATA/weights_Binding_PI3.txt",
                                             assay3="PIK3CG",
                                             ddG4="./DATA/weights_Binding_RAL.txt",
                                             assay4="RALGDS",
                                             ddG5="./DATA/weights_Binding_SOS.txt",
                                             assay5="SOS1",
                                             ddG6="./DATA/weights_Binding_K27.txt",
                                             assay6="DARPin K27",
                                             ddG7="./DATA/weights_Binding_K55.txt",
                                             assay7="DARPin K55",
                                             anno=final_distance_dc_anno_for_anno_mt_5,
                                             rect_input=rects,
                                             rect_alpha = rects_alpha, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure4a_point_all_ddG_alphabeta.pdf", device = cairo_pdf,height = 10,width=10)

  #4b
  krasddpcams__plot_RAL_invitro_cor_ITC(ddG="./DATA/weights_Binding_RAL.txt",
                                assay_sele = "RAL", colour_scheme)
  ggplot2::ggsave("./FIGURES/figure4b_RALGDS_ITC_invitro_ddG_cor.pdf", device = cairo_pdf,height = 40,width=40,units = "mm")

  #4c
  krasddpcams__plot_RAL_invitro_cor_GDI(ddG="./DATA/weights_Binding_RAL.txt",
                                assay_sele = "RAL", colour_scheme)
  ggplot2::ggsave("./FIGURES/figure4c_RALGDS_GDI_invitro_ddG_cor.pdf", device = cairo_pdf,height = 45,width=45,units = "mm")

  ######
  #Figure 5
  ######

  #5a
  figure5A<-krasddpcams__get_common_interface(anno_input=anno_final3)
  figure5A[["p"]]
  ggplot2::ggsave("./FIGURES/figure5a_common_interface.pdf", device = cairo_pdf,height = 4, width=7)

  #5b
  #top
  #to make figure5b, use chimeraX color binding interface residues
  #binding interface residues of 5 binding partners
  # paste0(figure5A[["matrix"]][assay=="K55"&binding=="yes",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="K55"&binding=="common",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="K27"&binding=="yes",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="K27"&binding=="common",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="RAL"&binding=="yes",Pos]+200,collapse=",")
  # paste0(figure5A[["matrix"]][assay=="RAL"&binding=="common",Pos]+200,collapse=",")
  # paste0(figure5A[["matrix"]][assay=="SOS"&binding=="yes",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="SOS"&binding=="common",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="PI3"&binding=="yes",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="PI3"&binding=="common",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="RAF"&binding=="yes",Pos],collapse=",")
  # paste0(figure5A[["matrix"]][assay=="RAF"&binding=="common",Pos],collapse=",")

  #5b
  #bottom
  #RAF
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                                   input_PDB = "./DATA/pdb6vjj.ent",
                                   chain_KRAS = "A",
                                   output_PDB_file="./FIGURES/figure5b_RAF_binding_ab_ddG_6vjj.pdb")

  #PI3
  weighted_mean_abs_ddG_PI3<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_PI3.txt",
                                                       assay_sele = "PI3")
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_PI3,
                                   input_PDB = "./DATA/pdb1he8.ent",
                                   chain_KRAS = "B",
                                   output_PDB_file="./FIGURES/figure5c_PI3_binding_ab_ddG_1he8.pdb")

  #RAL
  weighted_mean_abs_ddG_RAL<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_RAL.txt",
                                                       assay_sele = "RAL")
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_RAL,
                                   input_PDB = "./DATA/pdb1lfd.ent",
                                   chain_KRAS = "B",
                                   Pos_correction=200,
                                   output_PDB_file="./FIGURES/figure5c_RAL_binding_ab_ddG_1lfd.pdb")

  #K55
  weighted_mean_abs_ddG_K55<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_K55.txt",
                                                       assay_sele = "K55")
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_K55,
                                   input_PDB = "./DATA/pdb5o2t.ent",
                                   chain_KRAS = "A",
                                   output_PDB_file="./FIGURES/figure5c_K55_binding_ab_ddG_5o2t.pdb")

  #SOS
  weighted_mean_abs_ddG_SOS<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_SOS.txt",
                                                       assay_sele = "SOS")
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_SOS,
                                   input_PDB = "./DATA/1nvw.pdb",
                                   chain_KRAS = "R",
                                   output_PDB_file="./FIGURES/figure5c_SOS_binding_ab_ddG_1nvw_R.pdb")
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_SOS,
                                   input_PDB = "./FIGURES/figure5c_SOS_binding_ab_ddG_1nvw_R.pdb",
                                   chain_KRAS = "Q",
                                   output_PDB_file="./FIGURES/figure5c_SOS_binding_ab_ddG_1nvv_RQ.pdb")

  #K27
  weighted_mean_abs_ddG_K27<-krasddpcams__get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_K27.txt",
                                                       assay_sele = "K27")
  krasddpcams__weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_K27,
                                   input_PDB = "./DATA/pdb5o2s.ent",
                                   chain_KRAS = "A",
                                   output_PDB_file="./FIGURES/figure5c_K27_binding_ab_ddG_5o2s.pdb")
  #to make figure5b, use chimeraX. key:color bfactor sel range -2,2

  #5c
  #to make figure5c, use chimeraX. key:color sel gray
  #binding interface of binding partners
  pdb_list <- list(
    RAF = list(pdb = "./DATA/pdb6vjj.ent", chain_query = "B", chain_target = "A"),
    K27 = list(pdb = "./DATA/pdb5o2s.ent", chain_query = "B", chain_target = "A"),
    PI3 = list(pdb = "./DATA/pdb1he8.ent", chain_query = "A", chain_target = "B"),
    RAL = list(pdb = "./DATA/pdb1lfd.ent", chain_query = "A", chain_target = "B"),
    K55 = list(pdb = "./DATA/pdb5o2t.ent", chain_query = "B", chain_target = "A"),
    SOS1 = list(pdb = "./DATA/1nvw.pdb", chain_query = "S", chain_target = "R"),
    SOS2 = list(pdb = "./DATA/1nvw.pdb", chain_query = "S", chain_target = "Q"))
  figure5c <- krasddpcams__get_bindingresidue_bindingpartners(pdb_list)
  # paste0(figure5c[binding_partner=="RAF"&scHAmin_ligand<5,Pos],collapse=",")
  # paste0(figure5c[binding_partner=="PI3"&scHAmin_ligand<5,Pos],collapse=",")
  # paste0(figure5c[binding_partner=="RAL"&scHAmin_ligand<5,Pos],collapse=",")
  # #SOS1 chain R
  # paste0(figure5c[binding_partner=="SOS1"&scHAmin_ligand<5,Pos],collapse=",")
  # #SOS1 chain Q
  # paste0(figure5c[binding_partner=="SOS2"&scHAmin_ligand<5,Pos],collapse=",")
  # paste0(figure5c[binding_partner=="K27"&scHAmin_ligand<5,Pos],collapse=",")
  # paste0(figure5c[binding_partner=="K55"&scHAmin_ligand<5,Pos],collapse=",")

  #5d
  krasddpcams__plot_ddG_effector_bf_heatmap(ddG1="./DATA/weights_Binding_RAF.txt",
                                    assay1="RAF1",
                                    ddG2="./DATA/weights_Binding_PI3.txt",
                                    assay2="PIK3CG",
                                    ddG3="./DATA/weights_Binding_RAL.txt",
                                    assay3="RALGDS",
                                    anno=anno_final3,
                                    bind1=scHAmin_ligand_RAF,
                                    bind2=scHAmin_ligand_PI3,
                                    bind3=scHAmin_ligand_RAL, wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure5d_ddG_effector_bf_heatmap.pdf", device = cairo_pdf,height = 4, width=7)

  ######
  #Figure 6
  ######

  #6a
  all_allosteric_site_list<-krasddpcams__plot_major_allosteric_site_combined_threshold(ddG1="./DATA/weights_Binding_RAF.txt",
                                                                               assay1="RAF",
                                                                               ddG2="./DATA/weights_Binding_PI3.txt",
                                                                               assay2="PI3",
                                                                               ddG3="./DATA/weights_Binding_RAL.txt",
                                                                               assay3="RAL",
                                                                               ddG4="./DATA/weights_Binding_SOS.txt",
                                                                               assay4="SOS",
                                                                               ddG5="./DATA/weights_Binding_K27.txt",
                                                                               assay5="K27",
                                                                               ddG6="./DATA/weights_Binding_K55.txt",
                                                                               assay6="K55",
                                                                               anno=anno_final3,
                                                                               allosteric_list = all_allosteric_site_list,
                                                                               rect_input=rects, colour_scheme)
  all_allosteric_site_list[["p"]]
  ggplot2::ggsave("./FIGURES/figure6a_major_allosteric_site_combined_threshold.pdf", device = cairo_pdf,height = 4, width=7)

  #6b
  krasddpcams__plot_ddG_all_allosteric_heatmap(ddG1="./DATA/weights_Binding_RAF.txt",
                                       assay1="RAF1",
                                       ddG2="./DATA/weights_Binding_PI3.txt",
                                       assay2="PIK3CG",
                                       ddG3="./DATA/weights_Binding_RAL.txt",
                                       assay3="RALGDS",
                                       ddG4="./DATA/weights_Binding_SOS.txt",
                                       assay4="SOS1",
                                       ddG5="./DATA/weights_Binding_K27.txt",
                                       assay5="DARPin K27",
                                       ddG6="./DATA/weights_Binding_K55.txt",
                                       assay6="DARPin K55",
                                       anno=anno_final3,
                                       list1 = all_allosteric_site_list[["RAF"]],
                                       list2 = all_allosteric_site_list[["PI3"]],
                                       list3 = all_allosteric_site_list[["RAL"]],
                                       list4 = all_allosteric_site_list[["SOS"]],
                                       list5 = all_allosteric_site_list[["K27"]],
                                       list6 = all_allosteric_site_list[["K55"]], wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figure6b_ddG_all_allosteric_heatmap.pdf", device = cairo_pdf,height = 4, width=7)

  ######
  #ED figure 1
  ######

  #S1a, b
  d <- krasddpcams__plot_fitness_correlation_blocks(stability_nor_df,"Abundance", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS1a_fitness_correlation_blocks_abundance.pdf", d, device = cairo_pdf,height = 4, width=4)
  d <- krasddpcams__plot_fitness_correlation_blocks(RAF_nor_df,"RAF1 binding", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS1b_fitness_correlation_blocks_binding.pdf", d, device = cairo_pdf,height = 4, width=4)

  #S1c
  krasddpcams__plot_kuriyan_cor(kuriyan_list="./DATA/kuriyan_list.RData", colour_scheme, all_data_pos)
  ggplot2::ggsave("./FIGURES/figureS1c_kuriyan_fitness_cor.pdf", device = cairo_pdf,height = 40,width=40,units = "mm")

  #S1d
  krasddpcams__plot_cor_gr_tecan_binding(tecandata = "./DATA/20230616_tecanA_1.txt",
                            plasmidid = "./DATA/plasmid_plate_0615.csv",
                            single = all_data_pos[Nham_aa<=1,], colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS1d_tecantest_fitness_cor.pdf", device = cairo_pdf,height = 40,width=40,units = "mm")

  #S1e
  krasddpcams__plot_distribution_fitness(input=all_data_pos,assay_name = "RAF",anno_input = anno_final3, wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS1e_distribution_fitness_cor.pdf", device = cairo_pdf,height = 20,width=80,units = "mm")

  #S1f
  Folding_pre_ob_fitness<-krasddpcams__merge_ddGf_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                                              folding="./DATA/weights_Folding.txt",
                                              block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                              block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                              block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                              wt_aa_input=wt_aa)
  krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness,
                           fold_n=1,
                           mochi_parameters = "./DATA/linears_weights_Abundance1.txt",
                           phenotypen = 1,RT=0.001987*(273+30),bin_input = 50)
  ggplot2::ggsave("./FIGURES/figureS1f1_block1_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness,
                           fold_n=1,
                           mochi_parameters = "./DATA/linears_weights_Abundance2.txt",
                           phenotypen = 2,RT=0.001987*(273+30),bin_input = 50)
  ggplot2::ggsave("./FIGURES/figureS1f2_block2_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness,
                           fold_n=1,
                           mochi_parameters = "./DATA/linears_weights_Abundance3.txt",
                           phenotypen = 3,RT=0.001987*(273+30),bin_input = 50)
  ggplot2::ggsave("./FIGURES/figureS1f3_block3_ddGf_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

  #S1g
  krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = Folding_pre_ob_fitness,phenotypen = 1)
  ggplot2::ggsave("./FIGURES/figureS1g1_block1_pre_ob_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = Folding_pre_ob_fitness,phenotypen = 2)
  ggplot2::ggsave("./FIGURES/figureS1g2_block2_pre_ob_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = Folding_pre_ob_fitness,phenotypen = 3)
  ggplot2::ggsave("./FIGURES/figureS1g3_block3_pre_ob_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

  #S1h
  krasddpcams__plot_pre_12_ob_3_fitness_heldout_cor(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                                       observed="./DATA/predicted_phenotypes_all.txt",
                                       phenotypen=1)
  ggplot2::ggsave("./FIGURES/figureS1h1_block1_pre12_ob3_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_heldout_cor(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                               observed="./DATA/predicted_phenotypes_all.txt",
                               phenotypen=2)
  ggplot2::ggsave("./FIGURES/figureS1h2_block2_pre12_ob3_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_heldout_cor(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                               observed="./DATA/predicted_phenotypes_all.txt",
                               phenotypen=3)
  ggplot2::ggsave("./FIGURES/figureS1h3_block3_pre12_ob3_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

  #S1i
  RAF_pre_ob_fitness<-krasddpcams__merge_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                                                 block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                 block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                 block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                                 assay_sele="RAF1",
                                                 wt_aa_input=wt_aa)

  Cairo::CairoPDF(file = "./FIGURES/figureS1i1_block1_ddGf_ddGb_ob_pre_fitness.pdf")
  krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input=RAF_pre_ob_fitness,
                                       folding_assay = Abundance1,
                                       binding_assay = Binding1_RAF,
                                       RT=0.001987*(273+30),
                                       mochi_parameters = "./DATA/linears_weights_Binding1_RAF.txt", colour_scheme)
  dev.off()

  Cairo::CairoPDF(file = "./FIGURES/figureS1i2_block2_ddGf_ddGb_ob_pre_fitness.pdf")
  krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input=RAF_pre_ob_fitness,
                                       folding_assay = Abundance2,
                                       binding_assay = Binding2_RAF,
                                       RT=0.001987*(273+30),
                                       mochi_parameters = "./DATA/linears_weights_Binding2_RAF.txt", colour_scheme)
  dev.off()

  Cairo::CairoPDF(file = "./FIGURES/figureS1i3_block3_ddGf_ddGb_ob_pre_fitness.pdf")
  ggplot2::ggsave("./FIGURES/figureS1i3_block2_ddGf_ddGb_ob_pre_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness(binding_input=RAF_pre_ob_fitness,
                                       folding_assay = Abundance3,
                                       binding_assay = Binding3_RAF,
                                       RT=0.001987*(273+30),
                                       mochi_parameters = "./DATA/linears_weights_Binding3_RAF.txt", colour_scheme)
  dev.off()

  #S1j
  krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF_pre_ob_fitness,phenotypen = 13)
  ggplot2::ggsave("./FIGURES/figureS1j1_block1_pre_ob_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF_pre_ob_fitness,phenotypen = 14)
  ggplot2::ggsave("./FIGURES/figureS1j2_block2_pre_ob_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF_pre_ob_fitness,phenotypen = 15)
  ggplot2::ggsave("./FIGURES/figureS1j3_block3_pre_ob_fitness.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  #S1k
  krasddpcams__plot_pre_12_ob_3_fitness_heldout_cor(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                               observed="./DATA/predicted_phenotypes_all.txt",
                               phenotypen=13)
  ggplot2::ggsave("./FIGURES/figureS1k1_block1_pre12_ob3_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_heldout_cor(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                               observed="./DATA/predicted_phenotypes_all.txt",
                               phenotypen=14)
  ggplot2::ggsave("./FIGURES/figureS1k2_block2_pre12_ob3_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=70,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_heldout_cor(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                               observed="./DATA/predicted_phenotypes_all.txt",
                               phenotypen=15)
  ggplot2::ggsave("./FIGURES/figureS1k3_block3_pre12_ob3_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

  ######
  #ED figure 2
  ######

  #S2a
  krasddpcams__plot_beta_strand_propagation_heatmap(ddG="./DATA/weights_Binding_RAF.txt",
                                             assay_sele = "RAF", wt_aa, colour_scheme)+ggplot2::theme(strip.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.4,0,b = 0.2,0, "mm")),
                                                                       strip.text = ggplot2::element_text(size=5))+ggplot2::coord_fixed(ratio=2)
  ggplot2::ggsave("./FIGURES/figureS2a_beta_strand_propagation_heatmap.pdf", device = cairo_pdf,height = 4,width=10)

  #S2b
  krasddpcams__plot_ddG_ss(ddG="./DATA/weights_Binding_RAF.txt",
                   anno=anno_final3,
                   assay_sele="RAF",
                   rect_input=rects,
                   rect_alpha=rects_alpha_all, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS2b_ddG_ss.pdf", device = cairo_pdf,height = 3,width=7)

  #S2b - test
  krasddpcams__plot_fisher_test(target_group1="b1",
                   other_group2="others",
                   condition1="allosteric",
                   condition2="others",
                   g1c1 = 29,
                   g1c2 = 142,
                   g2c1 = 363-29,
                   g2c2 = 2712-142,
                   x_y = c("ddG", "structure"))
  #test secondary structures one by one
  ggplot2::ggsave("./FIGURES/figureS2b_fisher_test.pdf", device = cairo_pdf,height = 4,width=4)

  #S2c
  krasddpcams__plot_allosteric_mut_distance_2color(ddG="./DATA/weights_Binding_RAF.txt",
                                            weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                                            anno=anno_final3,
                                            assay_sele="RAF",
                                            rect_input=rects, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS2c_allosteric_mut_distance_2color.pdf", device = cairo_pdf,height = 4,width=4)

  #S2defg
  Sotorasib_dis<-krasddpcams__get_sotorasib_pocket(input_file = "./DATA/pdb6oim.ent")
  ggpubr::ggarrange(
    krasddpcams__plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                              anno=anno_final3,
                                              pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                              pocket = "pocket1",
                                              assay_sele="RAF",
                                              pocket2=Sotorasib_dis, colour_scheme = colour_scheme),
    krasddpcams__plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                              anno=anno_final3,
                                              pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                              pocket = "pocket2",
                                              assay_sele="RAF",
                                              pocket2=Sotorasib_dis, colour_scheme = colour_scheme),
    krasddpcams__plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                              anno=anno_final3,
                                              pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                              pocket = "pocket3",
                                              assay_sele="RAF",
                                              pocket2=Sotorasib_dis, colour_scheme = colour_scheme),
    krasddpcams__plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                              anno=anno_final3,
                                              pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                              pocket = "pocket4",
                                              assay_sele="RAF",
                                              pocket2=Sotorasib_dis, colour_scheme = colour_scheme),
    ncol=1,align="v"
  )
  ggplot2::ggsave("./FIGURES/figureS2def_allosteric_ddG_dis_pocket_color.pdf", device = cairo_pdf,height = 7,width=4)

  #S2h
  krasddpcams__plot_cor_gr_tecan_binding_MoCHi_prediction_double(tecandata = "./DATA/20230616_tecanA_1.txt",
                                                    plasmidid = "./DATA/plasmid_plate_0615.csv",
                                                    ddG1 = "./DATA/weights_Folding.txt",
                                                    ddG2 = "./DATA/weights_Binding_RAF.txt", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS2h_tecantest_doublemutations_cor.pdf", device = cairo_pdf,height = 45,width=45,units = "mm")

  ######
  #ED figure 3
  ######

  #S3a
  d<-krasddpcams__plot_fitness_correlation_blocks(PI3_nor_df,"PIK3CG binding", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS3a_fitness_correlation_blocks_PIK3CG.pdf", d, device = cairo_pdf,height = 4, width=4)
  d<-krasddpcams__plot_fitness_correlation_blocks(RAL_nor_df,"RALGDS binding", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS3a_fitness_correlation_blocks_RALGDS.pdf", d, device = cairo_pdf,height = 4, width=4)
  d<-krasddpcams__plot_fitness_correlation_blocks(SOS_nor_df,"SOS1 binding", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS3a_fitness_correlation_blocks_SOS1.pdf", d, device = cairo_pdf,height = 4, width=4)
  d<-krasddpcams__plot_fitness_correlation_blocks(K27_nor_df,"DARPin K27 binding", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS3a_fitness_correlation_blocks_K27.pdf", d, device = cairo_pdf,height = 4, width=4)
  d<-krasddpcams__plot_fitness_correlation_blocks(K55_nor_df,"DARPin K55 binding", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS3a_fitness_correlation_blocks_K55.pdf", d, device = cairo_pdf,height = 4, width=4)

  #S3b
  krasddpcams__plot2d_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                                  block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                  block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                  block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData",
                                  assay_sele = "PIK3CG",
                                  wt_aa_input=wt_aa)
  ggplot2::ggsave("./FIGURES/figureS3b1_PIK3CG_pre_ob_fitness_R2.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                             block1_dimsum_df = "./DATA/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                             block2_dimsum_df = "./DATA/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                             block3_dimsum_df = "./DATA/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData",
                             assay_sele = "RALGDS",
                             wt_aa_input=wt_aa)
  ggplot2::ggsave("./FIGURES/figureS3b2_RALGDS_pre_ob_fitness_R2.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                             block1_dimsum_df = "./DATA/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                             block2_dimsum_df = "./DATA/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                             block3_dimsum_df = "./DATA/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData",
                             assay_sele = "SOS1",
                             wt_aa_input=wt_aa)
  ggplot2::ggsave("./FIGURES/figureS3b3_SOS1_pre_ob_fitness_R2.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                             block1_dimsum_df = "./DATA/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                             block2_dimsum_df = "./DATA/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                             block3_dimsum_df = "./DATA/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
                             assay_sele = "DARPin K27",
                             wt_aa_input=wt_aa)
  ggplot2::ggsave("./FIGURES/figureS3b4_K27_pre_ob_fitness_R2.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot2d_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                             block1_dimsum_df = "./DATA/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                             block2_dimsum_df = "./DATA/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                             block3_dimsum_df = "./DATA/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
                             assay_sele = "DARPin K55",
                             wt_aa_input=wt_aa)
  ggplot2::ggsave("./FIGURES/figureS3b5_K55_pre_ob_fitness_R2.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

  #S3c
  krasddpcams__plot_pre_12_ob_3_fitness_cor_all_heldout(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                                   observed="./DATA/predicted_phenotypes_all.txt",
                                   assay_sele = "PIK3CG",
                                   r_y=0.2)
  ggplot2::ggsave("./FIGURES/figureS3c1_PI3_pre_ob_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_cor_all_heldout(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                                   observed="./DATA/predicted_phenotypes_all.txt",
                                   assay_sele = "RALGDS",
                                   r_y=0.2)
  ggplot2::ggsave("./FIGURES/figureS3c2_RALGDS_pre_ob_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_cor_all_heldout(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                                   observed="./DATA/predicted_phenotypes_all.txt",
                                   assay_sele = "SOS1",
                                   r_y=0.2)
  ggplot2::ggsave("./FIGURES/figureS3c3_SOS1_pre_ob_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_cor_all_heldout(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                                   observed="./DATA/predicted_phenotypes_all.txt",
                                   assay_sele = "DARPin K27",
                                   r_y=0.2)
  ggplot2::ggsave("./FIGURES/figureS3c4_K27_pre_ob_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")
  krasddpcams__plot_pre_12_ob_3_fitness_cor_all_heldout(prediction="./DATA/Revision/norep3/Mochi/predicted_phenotypes_all.txt",
                                   observed="./DATA/predicted_phenotypes_all.txt",
                                   assay_sele = "DARPin K55",
                                   r_y=0.2)
  ggplot2::ggsave("./FIGURES/figureS3c5_K55_pre_ob_fitness_heldout.pdf", device = cairo_pdf,height = 35,width=60,units = "mm")

  ######
  #ED figure 4
  ######

  #s4a
  krasddpcams__plot_figureS4C_ddG_heatmap_1124(ddG1="./DATA/weights_Folding.txt",
                                  assay1="folding",
                                  ddG2="./DATA/weights_Binding_RAF.txt",
                                  assay2="RAF1",
                                  ddG3="./DATA/weights_Binding_PI3.txt",
                                  assay3="PIK3CG",
                                  ddG4="./DATA/weights_Binding_RAL.txt",
                                  assay4="RALGDS",
                                  ddG5="./DATA/weights_Binding_SOS.txt",
                                  assay5="SOS1",
                                  ddG6="./DATA/weights_Binding_K27.txt",
                                  assay6="DARPin K27",
                                  ddG7="./DATA/weights_Binding_K55.txt",
                                  assay7="DARPin K55", wt_aa, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS4a_ddG_heatmap_1124.pdf", device = cairo_pdf,height = 10,width=10)

  #s4b
  ggpubr::ggarrange(krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                           weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                           weighted_meab_abs_ddGb=weighted_mean_abs_ddG_PI3,
                                                           anno_input=anno_final3,
                                                           bind="PI3",
                                                           bind_ligand = "scHAmin_ligand_PI3")+ggplot2::ggtitle("PIK3CG"),
            krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                           weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                           weighted_meab_abs_ddGb=weighted_mean_abs_ddG_RAL,
                                                           anno_input=anno_final3,
                                                           bind="RAL",
                                                           bind_ligand = "scHAmin_ligand_RAL")+ggplot2::ggtitle("RALGDS"),
            krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                           weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                           weighted_meab_abs_ddGb=weighted_mean_abs_ddG_K55,
                                                           anno_input=anno_final3,
                                                           bind="SOS",
                                                           bind_ligand = "scHAmin_ligand_SOS")+ggplot2::ggtitle("SOS1"),
            krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                           weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                           weighted_meab_abs_ddGb=weighted_mean_abs_ddG_K27,
                                                           anno_input=anno_final3,
                                                           bind="K27",
                                                           bind_ligand = "scHAmin_ligand_K27")+ggplot2::ggtitle("DARPin K27"),
            krasddpcams__plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                           weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                           weighted_meab_abs_ddGb=weighted_mean_abs_ddG_K55,
                                                           anno_input=anno_final3,
                                                           bind="K55",
                                                           bind_ligand = "scHAmin_ligand_K55")+ggplot2::ggtitle("DARPin K55"),
            ncol=6
  )
  ggplot2::ggsave("./FIGURES/figureS4a_ddG_heatmap_1124.pdf", device = cairo_pdf,height = 4,width=10)

  ######
  #ED figure 5
  ######

  krasddpcams__plot_ddG_all_bf_heatmap(ddG1="./DATA/weights_Binding_RAF.txt",
                               assay1="RAF1",
                               ddG2="./DATA/weights_Binding_PI3.txt",
                               assay2="PIK3CG",
                               ddG3="./DATA/weights_Binding_RAL.txt",
                               assay3="RALGDS",
                               ddG4="./DATA/weights_Binding_SOS.txt",
                               assay4="SOS1",
                               ddG5="./DATA/weights_Binding_K27.txt",
                               assay5="DARPin K27",
                               ddG6="./DATA/weights_Binding_K55.txt",
                               assay6="DARPin K55",
                               anno=anno_final3,
                               bind1=scHAmin_ligand_RAF,
                               bind2=scHAmin_ligand_PI3,
                               bind3=scHAmin_ligand_RAL,
                               bind4=scHAmin_ligand_SOS,
                               bind5=scHAmin_ligand_K27,
                               bind6=scHAmin_ligand_K55, wt_aa, colour_scheme
  )+ggplot2::theme(strip.background = ggplot2::element_rect(colour="white",fill=NULL,linetype = NULL))
  ggplot2::ggsave("./FIGURES/figureS5_ddG_all_bf_heatmap.pdf", device = cairo_pdf,height = 7,width=7)

  ######
  #ED figure 6
  ######

  combined_threshold<-krasddpcams__get_combined_threshold(ddG1="./DATA/weights_Binding_RAF.txt",
                                             assay1="RAF",
                                             ddG2="./DATA/weights_Binding_PI3.txt",
                                             assay2="PI3",
                                             ddG3="./DATA/weights_Binding_RAL.txt",
                                             assay3="RAL",
                                             ddG4="./DATA/weights_Binding_SOS.txt",
                                             assay4="SOS",
                                             ddG5="./DATA/weights_Binding_K27.txt",
                                             assay5="K27",
                                             ddG6="./DATA/weights_Binding_K55.txt",
                                             assay6="K55",
                                             anno=anno_final3)
  krasddpcams__plot_allosteric_mutations_enrichment(ddG="./DATA/weights_Binding_RAF.txt",
                                            anno=anno_final3,
                                            assay_sele ="RAF",
                                            threshold=combined_threshold, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS6b_allosteric_mutations_enrichment_RAF.pdf", device = cairo_pdf,height = 3,width=5)

  krasddpcams__plot_allosteric_mutations_enrichment(ddG="./DATA/weights_Binding_PI3.txt",
                                       anno=anno_final3,
                                       assay_sele ="PI3",
                                       threshold=combined_threshold, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS6b_allosteric_mutations_enrichment_PI3.pdf", device = cairo_pdf,height = 3,width=5)

  krasddpcams__plot_allosteric_mutations_enrichment(ddG="./DATA/weights_Binding_RAL.txt",
                                       anno=anno_final3,
                                       assay_sele ="RAL",
                                       threshold=combined_threshold, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS6b_allosteric_mutations_enrichment_RAL.pdf", device = cairo_pdf,height = 3,width=5)

  krasddpcams__plot_allosteric_mutations_enrichment(ddG="./DATA/weights_Binding_SOS.txt",
                                       anno=anno_final3,
                                       assay_sele ="SOS",
                                       threshold=combined_threshold, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS6b_allosteric_mutations_enrichment_SOS.pdf", device = cairo_pdf,height = 3,width=5)

  krasddpcams__plot_allosteric_mutations_enrichment(ddG="./DATA/weights_Binding_K27.txt",
                                       anno=anno_final3,
                                       assay_sele ="K27",
                                       threshold=combined_threshold, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS6b_allosteric_mutations_enrichment_K27.pdf", device = cairo_pdf,height = 3,width=5)

  krasddpcams__plot_allosteric_mutations_enrichment(ddG="./DATA/weights_Binding_K55.txt",
                                       anno=anno_final3,
                                       assay_sele ="K55",
                                       threshold=combined_threshold, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS6b_allosteric_mutations_enrichment_K55.pdf", device = cairo_pdf,height = 3,width=5)

  ######
  #ED figure 7
  ######

  krasddpcams__plot_allosteric_mutations_all(ddG="./DATA/weights_Binding_RAF.txt",
                                     anno=anno_final3,
                                     assay_sele ="RAF",
                                     threshold=combined_threshold,
                                     rect_input=rects,
                                     rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS7_allosteric_mutations_all_RAF.pdf", device = cairo_pdf,height = 4,width=10)

  krasddpcams__plot_allosteric_mutations_all(ddG="./DATA/weights_Binding_PI3.txt",
                                     anno=anno_final3,
                                     assay_sele ="PI3",
                                     threshold=combined_threshold,
                                     rect_input=rects,
                                     rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS7_allosteric_mutations_all_PI3.pdf", device = cairo_pdf,height = 4,width=10)

  krasddpcams__plot_allosteric_mutations_all(ddG="./DATA/weights_Binding_RAL.txt",
                                     anno=anno_final3,
                                     assay_sele ="RAL",
                                     threshold=combined_threshold,
                                     rect_input=rects,
                                     rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS7_allosteric_mutations_all_RAL.pdf", device = cairo_pdf,height = 4,width=10)

  krasddpcams__plot_allosteric_mutations_all(ddG="./DATA/weights_Binding_SOS.txt",
                                     anno=anno_final3,
                                     assay_sele ="SOS",
                                     threshold=combined_threshold,
                                     rect_input=rects,
                                     rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS7_allosteric_mutations_all_SOS.pdf", device = cairo_pdf,height = 4,width=10)

  krasddpcams__plot_allosteric_mutations_all(ddG="./DATA/weights_Binding_K27.txt",
                                     anno=anno_final3,
                                     assay_sele ="K27",
                                     threshold=combined_threshold,
                                     rect_input=rects,
                                     rect_alpha=rects_alpha,y_min=-2.5,y_max=2.5, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS7_allosteric_mutations_all_K27.pdf", device = cairo_pdf,height = 4,width=10)

  krasddpcams__plot_allosteric_mutations_all(ddG="./DATA/weights_Binding_K55.txt",
                                     anno=anno_final3,
                                     assay_sele ="K55",
                                     threshold=combined_threshold,
                                     rect_input=rects,
                                     rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2, colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS7_allosteric_mutations_all_K55.pdf", device = cairo_pdf,height = 4,width=10)

  ######
  #ED figure 8
  ######

  #a heatmap
  krasddpcams__plot_heatmap_fitness_block1(GAP_nor_df_pos, wt_aa, colour_scheme)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(family="Courier",angle=90,size=5, vjust =.5,hjust =.5,margin=ggplot2::margin(0,-0.5,0,0,"mm")))
  ggplot2::ggsave("./FIGURES/figureS8a_GAP_block1_fitness.pdf", device = cairo_pdf,height = 120,width=120,units = "mm")

  #b fitness correlation
  d <- krasddpcams__plot_fitness_correlation_block1(GAP_nor_df_pos)
  ggplot2::ggsave("./FIGURES/figureS8b_GAP_block1_fitness_cor.pdf", d, device = cairo_pdf,height = 60,width=60,units = "mm")

  #c growthrate of GAP and RAF
  krasddpcams__plot_density_GAP_RAF_single_gr_6mut(GAP_fitness = GAP_nor_df_pos,
                                 RAF_fitness = RAF_nor_df_pos,
                                 driven_mt1 = "G12C",
                                 driven_mt2 = "G12D",
                                 driven_mt3 = "G12V",
                                 driven_mt4 = "Q61L",
                                 driven_mt5 = "Q61H",
                                 driven_mt6 = "G13D", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8c_GAP_RAF_block1_single_growthrate_density_plot_6mut.pdf", device = cairo_pdf,height = 40,width=150,units = "mm")

  # d fitness driven mut vs wt
  krasddpcams__plot_drivenmt_wt_fitness_color_fig3(input = RAF_nor_df_pos,
                                      driven_mt = "G12C",
                                      ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                      anno=anno_final3,
                                      assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8d_RAF_block1_double_fitness_G12CWT_cor.pdf", device = cairo_pdf,height = 50,width=120,units = "mm")
  krasddpcams__plot_drivenmt_wt_fitness_color_fig3(input = RAF_nor_df_pos,
                                      driven_mt = "G12D",
                                      ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                      anno=anno_final3,
                                      assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8e_RAF_block1_double_fitness_G12DWT_cor.pdf", device = cairo_pdf,height = 50,width=120,units = "mm")
  krasddpcams__plot_drivenmt_wt_fitness_color_fig3(input = RAF_nor_df_pos,
                                      driven_mt = "G12V",
                                      ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                      anno=anno_final3,
                                      assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8f_RAF_block1_double_fitness_G12VWT_cor.pdf", device = cairo_pdf,height = 50,width=120,units = "mm")
  krasddpcams__plot_drivenmt_wt_fitness_color_fig3(input = GAP_nor_df_pos,
                                      driven_mt = "G12C",
                                      ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                      anno=anno_final3,
                                      assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8g_GAP_block1_double_fitness_G12CWT_cor.pdf", device = cairo_pdf,height = 50,width=120,units = "mm")
  krasddpcams__plot_drivenmt_wt_fitness_color_fig3(input = GAP_nor_df_pos,
                                      driven_mt = "G12D",
                                      ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                      anno=anno_final3,
                                      assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8h_GAP_block1_double_fitness_G12DWT_cor.pdf", device = cairo_pdf,height = 50,width=120,units = "mm")
  krasddpcams__plot_drivenmt_wt_fitness_color_fig3(input = GAP_nor_df_pos,
                                      driven_mt = "G12V",
                                      ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                      anno=anno_final3,
                                      assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS8i_GAP_block1_double_fitness_G12VWT_cor.pdf", device = cairo_pdf,height = 50,width=120,units = "mm")

  ######
  #ED figure 9
  ######

  #a fitness heatmap of block1 of fullRAF1 bindingPCA
  krasddpcams__plot_heatmap_fitness_block1(fullRAF_nor_df_pos, wt_aa, colour_scheme)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(family="Courier",angle=90,size=5, vjust =.5,hjust =.5,margin=ggplot2::margin(0,-0.5,0,0,"mm")))
  ggplot2::ggsave("./FIGURES/figureS9a_fullRAF1_block1_fitness.pdf", device = cairo_pdf,height = 120,width=120,units = "mm")

  #b fitness correlation across triplicates
  d <- krasddpcams__plot_fitness_correlation_block1(fullRAF_nor_df_pos)
  ggplot2::ggsave("./FIGURES/figureS9b_fullRAF1_block1_fitness_cor.pdf", d, device = cairo_pdf,height = 60,width=60,units = "mm")

  # #c fitness full RAF1 vs RAF1RBD
  # krasddpcams__plot_fitness_cor_RAF1_full_RBD(fullRAF_nor_df_pos,RAF_nor_df_pos, colour_scheme)
  # ggplot2::ggsave("./FIGURES/figureS9c_fullRAF1_RBD_block1_fitness_cor.pdf", device = cairo_pdf,height = 40,width=60,units = "mm")

  # #d mochi fitting Andre's codes
  # krasddpcams__plot_ddG_fullRAF_RBD("./DATA/Revision/fullRAF/Mochi/weights_Binding_RAF.txt",
  #                      "./DATA/Revision/fullRAF/Mochi/weights_Binding_FULLRAF.txt",
  #                      ci95_1 = 1,
  #                      ci95_2 =1, colour_scheme)
  # krasddpcams__plot_ddG_fullRAF_RBD("./DATA/Revision/fullRAF/Mochi/weights_Binding_RAF.txt",
  #                      "./DATA/Revision/fullRAF/Mochi/weights_Binding_FULLRAF.txt",
  #                      ci95_1 = 10000,
  #                      ci95_2 =10000, colour_scheme)
  # ggplot2::ggsave("./FIGURES/figureS9d_fullRAF1_RBD_block1_ddG_cor.pdf", device = cairo_pdf,height = 40,width=60,units = "mm")

  #c
  krasddpcams__plot_fitness_cor_RAF1_full_RBD_color_fig3(fullRAF_nor_df_pos,RAF_nor_df_pos,
                                            ddG_RAFRBD="./DATA/weights_Binding_RAF.txt",
                                            anno=anno_final3,
                                            assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS9c_fullRAF1_RBD_block1_fitness_cor_color.pdf", device = cairo_pdf,height = 40,width=120,units = "mm")

  #d
  krasddpcams__plot_ddG_fullRAF_RBD_color_fig3(ddG1 = "./DATA/Revision/fullRAF/Mochi/weights_Binding_RAF.txt",
                                  ddG2 = "./DATA/Revision/fullRAF/Mochi/weights_Binding_FULLRAF.txt",
                                  ddG_RAFRBD = "./DATA/weights_Binding_RAF.txt",
                                  anno=anno_final3,
                                  assay_sele="RAF", colour_scheme)
  ggplot2::ggsave("./FIGURES/figureS9d_fullRAF1_RBD_block1_ddG_cor_color.pdf", device = cairo_pdf,height = 40,width=120,units = "mm")


  ######
  #sup table 4
  ######

  abundancePCA <- stability_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  abundancePCA[,assay:="AbundancePCA"]
  rafPCA <- RAF_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  rafPCA[,assay:="BindingPCA RAF1RBD"]
  pi3PCA <- PI3_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  pi3PCA[,assay:="BindingPCA PIK3CGRBD"]
  ralPCA <- RAL_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  ralPCA[,assay:="BindingPCA RALGDSRBD"]
  sos1PCA <- SOS_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  sos1PCA[,assay:="BindingPCA SOS1"]
  k27PCA <- SOS_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  k27PCA[,assay:="BindingPCA DARPin K27"]
  k55PCA <- SOS_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  k55PCA[,assay:="BindingPCA DARPin K55"]

  fullrafPCA <- fullRAF_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  fullrafPCA[,assay:="BindingPCA full length RAF1"]

  rafgapPCA <- GAP_nor_df[,.(block, aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma, nor_gr, nor_gr_sigma, nor_fitness, nor_fitness_sigma)]
  rafgapPCA[,assay:="BindingPCA RAF1RBD coexpression GAP"]

  sup_table4 <- rbind(abundancePCA,
                      rafPCA,
                      pi3PCA,
                      ralPCA,
                      sos1PCA,
                      k27PCA,
                      k55PCA,
                      fullrafPCA,
                      rafgapPCA)
  sup_table4 <- as.data.frame(sup_table4)
  write.table(sup_table4,"./RESULTS/sup_table4_fitness_dimsum_updated.csv",sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)
  #copy to excel

  ######
  #sup table 5
  ######

  krasddpcams__write_tableS5(ddG1="./DATA/weights_Folding.txt",
                assay1="folding",
                ddG2="./DATA/weights_Binding_RAF.txt",
                assay2="RAF1",
                ddG3="./DATA/weights_Binding_PI3.txt",
                assay3="PIK3CG",
                ddG4="./DATA/weights_Binding_RAL.txt",
                assay4="RALGDS",
                ddG5="./DATA/weights_Binding_SOS.txt",
                assay5="SOS1",
                ddG6="./DATA/weights_Binding_K27.txt",
                assay6="DARPin K27",
                ddG7="./DATA/weights_Binding_K55.txt",
                assay7="DARPin K55",
                ddG8 = "./DATA/Revision/fullRAF/Mochi/weights_Binding_FULLRAF.txt",
                assay8 = "full length RAF1",
                output="./RESULTS/sup_table5_mochi_output.xlsx")

  ######
  #MAVEdb update
  ######

  stab_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "stab")
  write.csv(stab_MAVE[block=="block1"],"./RESULTS/abundance_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(stab_MAVE[block=="block2"],"./RESULTS/abundance_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(stab_MAVE[block=="block3"],"./RESULTS/abundance_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  RAF1RBD_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "RAF")
  write.csv(RAF1RBD_MAVE[block=="block1"],"./RESULTS/binding_RAF1RBD_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(RAF1RBD_MAVE[block=="block2"],"./RESULTS/binding_RAF1RBD_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(RAF1RBD_MAVE[block=="block3"],"./RESULTS/binding_RAF1RBD_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  RALGDS_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "RAL")
  write.csv(RALGDS_MAVE[block=="block1"],"./RESULTS/binding_RALGDSRBD_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(RALGDS_MAVE[block=="block2"],"./RESULTS/binding_RALGDSRBD_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(RALGDS_MAVE[block=="block3"],"./RESULTS/binding_RALGDSRBD_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  PI3_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "PI3")
  write.csv(PI3_MAVE[block=="block1"],"./RESULTS/binding_PI3KCGRBD_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(PI3_MAVE[block=="block2"],"./RESULTS/binding_PI3KCGRBD_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(PI3_MAVE[block=="block3"],"./RESULTS/binding_PI3KCGRBD_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  SOS1_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "SOS")
  write.csv(SOS1_MAVE[block=="block1"],"./RESULTS/binding_SOS1_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(SOS1_MAVE[block=="block2"],"./RESULTS/binding_SOS1_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(SOS1_MAVE[block=="block3"],"./RESULTS/binding_SOS1_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  K27_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "K27")
  write.csv(K27_MAVE[block=="block1"],"./RESULTS/binding_K27_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(K27_MAVE[block=="block2"],"./RESULTS/binding_K27_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(K27_MAVE[block=="block3"],"./RESULTS/binding_K27_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  K55_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=all_data_pos,assay_sele = "K55")
  write.csv(K55_MAVE[block=="block1"],"./RESULTS/binding_K55_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(K55_MAVE[block=="block2"],"./RESULTS/binding_K55_block2_fitness_MAVEdb.csv", quote = FALSE,row.names = F)
  write.csv(K55_MAVE[block=="block3"],"./RESULTS/binding_K55_block3_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  GAP_nor_df[,assay:="GAP"]
  GAP_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=GAP_nor_df,assay_sele = "GAP")
  write.csv(GAP_MAVE[block=="block1"],"./RESULTS/binding_GAP_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  fullRAF_nor_df[,assay:="full RAF"]
  fullRAF_MAVE<-krasddpcams__write_MAVEdb_scoretable(input=fullRAF_nor_df,assay_sele = "full RAF")
  write.csv(fullRAF_MAVE[block=="block1"],"./RESULTS/binding_fullRAF1_block1_fitness_MAVEdb.csv", quote = FALSE,row.names = F)

  ddG_f_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Folding.txt",
                                          assay_name = "folding")
  write.csv(ddG_f_MAVE,"./RESULTS/folding_ddG_MAVEdb.csv",quote = FALSE,row.names = F)
  ddG_b_RAF1RBD_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Binding_RAF.txt",
                                                  assay_name = "RAF1")
  write.csv(ddG_b_RAF1RBD_MAVE,"./RESULTS/binding_ddG_RAF1RBD_MAVEdb.csv",quote = FALSE,row.names = F)

  ddG_b_PI3_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Binding_PI3.txt",
                                              assay_name = "PI3")
  write.csv(ddG_b_PI3_MAVE,"./RESULTS/binding_ddG_PI3_MAVEdb.csv",quote = FALSE,row.names = F)

  ddG_b_RALGDS_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Binding_RAL.txt",
                                                 assay_name = "RAL")
  write.csv(ddG_b_RALGDS_MAVE,"./RESULTS/binding_ddG_RAL_MAVEdb.csv",quote = FALSE,row.names = F)

  ddG_b_SOS_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Binding_SOS.txt",
                                              assay_name = "SOS")
  write.csv(ddG_b_SOS_MAVE,"./RESULTS/binding_ddG_SOS_MAVEdb.csv",quote = FALSE,row.names = F)

  ddG_b_K27_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Binding_K27.txt",
                                              assay_name = "K27")
  write.csv(ddG_b_K27_MAVE,"./RESULTS/binding_ddG_K27_MAVEdb.csv",quote = FALSE,row.names = F)

  ddG_b_K55_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/weights_Binding_K55.txt",
                                              assay_name = "K55")
  write.csv(ddG_b_K55_MAVE,"./RESULTS/binding_ddG_K55_MAVEdb.csv",quote = FALSE,row.names = F)

  ddG_b_full_RAF1_MAVE<-krasddpcams__write_MAVEdb_scoretable_ddG(input="./DATA/Revision/fullRAF/Mochi/weights_Binding_FULLRAF.txt",
                                              assay_name = "full length RAF1")
  write.csv(ddG_b_full_RAF1_MAVE,"./RESULTS/binding_ddG_fullRAF1_MAVEdb.csv",quote = FALSE,row.names = F)

  ######
  #train dataset
  ######
  # new
  # folding
  model_folding_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                           block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                           block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                           block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                                           assay_sele = "folding",
                                                           wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_folding_pre_ob_fitness_R2_training.pdf",model_folding_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_folding_pre_ob_fitness_R2_test.pdf",model_folding_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_folding_new[["R_squared"]], file = "./RESULTS/model_folding_new.txt", sep = "\t")
  rm(model_folding_new)
  # RAF1
  model_RAF1_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                        block1_dimsum_df = "./DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                                        block2_dimsum_df = "./DATA/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                                        block3_dimsum_df = "./DATA/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
                                                        assay_sele = "RAF1",
                                                        wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_RAF1_pre_ob_fitness_R2_training.pdf",model_RAF1_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_RAF1_pre_ob_fitness_R2_test.pdf",model_RAF1_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_RAF1_new[["R_squared"]], file = "./RESULTS/model_RAF1_new.txt", sep = "\t")
  rm(model_RAF1_new)
  # PIK3CG
  model_PIK3CG_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                          block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                                          block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                                          block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData",
                                                          assay_sele = "PIK3CG",
                                                          wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_PIK3CG_pre_ob_fitness_R2_training.pdf",model_PIK3CG_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_PIK3CG_pre_ob_fitness_R2_test.pdf",model_PIK3CG_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_PIK3CG_new[["R_squared"]], file = "./RESULTS/model_PIK3CG_new.txt", sep = "\t")
  rm(model_PIK3CG_new)
  # RALGDS
  model_RALGDS_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                          block1_dimsum_df = "./DATA/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                                                          block2_dimsum_df = "./DATA/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                                                          block3_dimsum_df = "./DATA/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData",
                                                          assay_sele = "RALGDS",
                                                          wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_RALGDS_pre_ob_fitness_R2_training.pdf",model_RALGDS_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_RALGDS_pre_ob_fitness_R2_test.pdf",model_RALGDS_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_RALGDS_new[["R_squared"]], file = "./RESULTS/model_RALGDS_new.txt", sep = "\t")
  rm(model_RALGDS_new)
  # SOS1
  model_SOS1_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                        block1_dimsum_df = "./DATA/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                                                        block2_dimsum_df = "./DATA/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                                                        block3_dimsum_df = "./DATA/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData",
                                                        assay_sele = "SOS1",
                                                        wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_SOS1_pre_ob_fitness_R2_training.pdf",model_SOS1_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_SOS1_pre_ob_fitness_R2_test.pdf",model_SOS1_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_SOS1_new[["R_squared"]], file = "./RESULTS/model_SOS1_new.txt", sep = "\t")
  rm(model_SOS1_new)
  # DARPinK27
  model_DARPinK27_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                             block1_dimsum_df = "./DATA/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                                                             block2_dimsum_df = "./DATA/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                                                             block3_dimsum_df = "./DATA/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
                                                             assay_sele = "DARPin K27",
                                                             wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_DARPinK27_pre_ob_fitness_R2_training.pdf",model_DARPinK27_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_DARPinK27_pre_ob_fitness_R2_test.pdf",model_DARPinK27_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_DARPinK27_new[["R_squared"]], file = "./RESULTS/model_DARPinK27_new.txt", sep = "\t")
  rm(model_DARPinK27_new)
  # DARPinK55
  model_DARPinK55_new<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA/Revision/DATA_new_model//predicted_phenotypes_all.txt",
                                                             block1_dimsum_df = "./DATA/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                                                             block2_dimsum_df = "./DATA/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                                                             block3_dimsum_df = "./DATA/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
                                                             assay_sele = "DARPin K55",
                                                             wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_DARPinK55_pre_ob_fitness_R2_training.pdf",model_DARPinK55_new[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_newmodel_DARPinK55_pre_ob_fitness_R2_test.pdf",model_DARPinK55_new[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_DARPinK55_new[["R_squared"]], file = "./RESULTS/model_DARPinK55_new.txt", sep = "\t")
  rm(model_DARPinK55_new)
  # old
  # folding
  model_folding_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                           block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                           block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                           block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                                           assay_sele = "folding",
                                                           wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_folding_pre_ob_fitness_R2_test.pdf",model_folding_old[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_folding_pre_ob_fitness_R2_training.pdf",model_folding_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_folding_old[["R_squared"]], file = "./RESULTS/model_folding_old.txt", sep = "\t")
  rm(model_folding_old)
  # RAF1
  model_RAF1_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                        block1_dimsum_df = "./DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                                        block2_dimsum_df = "./DATA/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                                        block3_dimsum_df = "./DATA/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
                                                        assay_sele = "RAF1",
                                                        wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_RAF1_pre_ob_fitness_R2_test.pdf",model_RAF1_old[["test"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_RAF1_pre_ob_fitness_R2_training.pdf",model_RAF1_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_RAF1_old[["R_squared"]], file = "./RESULTS/model_RAF1_old.txt", sep = "\t")
  rm(model_RAF1_old)
  # PIK3CG
  model_PIK3CG_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                          block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                                          block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                                          block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData",
                                                          assay_sele = "PIK3CG",
                                                          wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_PIK3CG_pre_ob_fitness_R2_training.pdf",model_PIK3CG_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_PIK3CG_old[["R_squared"]], file = "./RESULTS/model_PIK3CG_old.txt", sep = "\t")
  rm(model_PIK3CG_old)
  # RALGDS
  model_RALGDS_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                          block1_dimsum_df = "./DATA/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                                                          block2_dimsum_df = "./DATA/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                                                          block3_dimsum_df = "./DATA/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData",
                                                          assay_sele = "RALGDS",
                                                          wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_RALGDS_pre_ob_fitness_R2_training.pdf",model_RALGDS_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_RALGDS_old[["R_squared"]], file = "./RESULTS/model_RALGDS_old.txt", sep = "\t")
  rm(model_RALGDS_old)
  # SOS1
  model_SOS1_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                        block1_dimsum_df = "./DATA/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                                                        block2_dimsum_df = "./DATA/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                                                        block3_dimsum_df = "./DATA/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData",
                                                        assay_sele = "SOS1",
                                                        wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_SOS1_pre_ob_fitness_R2_training.pdf",model_SOS1_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_SOS1_old[["R_squared"]], file = "./RESULTS/model_SOS1_old.txt", sep = "\t")
  rm(model_SOS1_old)
  # DARPinK27
  model_DARPinK27_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                             block1_dimsum_df = "./DATA/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                                                             block2_dimsum_df = "./DATA/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                                                             block3_dimsum_df = "./DATA/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
                                                             assay_sele = "DARPin K27",
                                                             wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_DARPinK27_pre_ob_fitness_R2_training.pdf",model_DARPinK27_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_DARPinK27_old[["R_squared"]], file = "./RESULTS/model_DARPinK27_old.txt", sep = "\t")
  rm(model_DARPinK27_old)
  # DARPinK55
  model_DARPinK55_old<-krasddpcams__plot2d_ddGb_ob_pre_fitness_train_test(prediction="./DATA//predicted_phenotypes_all.txt",
                                                             block1_dimsum_df = "./DATA/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                                                             block2_dimsum_df = "./DATA/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                                                             block3_dimsum_df = "./DATA/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
                                                             assay_sele = "DARPin K55",
                                                             wt_aa_input=wt_aa)
  #ggplot2::ggsave("./Figures_revision/figure_oldmodel_DARPinK55_pre_ob_fitness_R2_training.pdf",model_DARPinK55_old[["training"]], device = cairo_pdf,height = 35,width=60,units = "mm")
  fwrite(model_DARPinK55_old[["R_squared"]], file = "./RESULTS/model_DARPinK55_old.txt", sep = "\t")
  rm(model_DARPinK55_old)
  ############
  ### plot r2 of two models
  ############
  krasddpcams__plot_model_performance_comparison_R2(phenotype=c("folding","RAF1","PIK3CG","RALGDS","SOS1","DARPinK27","DARPinK55"),
                                                    model_type=c("old","new"),
                                                    colour_scheme)
  ggsave("./Figures_revision/figure_model_R2_comparison.pdf", device = cairo_pdf,height = 90,width=160,units = "mm")
}


