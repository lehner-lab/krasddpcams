###beginning
library(ggpubr)
library(GGally)
library(bio3d)
library(ggplot2)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(seqinr)
library(caTools)
library(pROC)
library(readr)
library(knitr)
library(PRROC)
library(rgl)
library(ggVennDiagram)
library(ggstatsplot)
library(plot3D)
library("devtools")
library(roxygen2)
library(data.table)
library(ROCR)
library(openxlsx)
library(ggrepel)
library(Cairo)
######
setwd("/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/")
#
library("devtools")
library(roxygen2)
#create("krasddpcams")
setwd("./krasddpcams")
document()
setwd("..")
library("krasddpcams")
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
#loading data
######
stability_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                 block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                 block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData")
RAF_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")
PI3_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData")
RAL_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData")
SOS_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData")
K27_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData")
K55_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData")
stab <- stability_nor_df
RAF <- RAF_nor_df
K27 <- K27_nor_df
PI3 <- PI3_nor_df
RAL <- RAL_nor_df
SOS <- SOS_nor_df
K55 <- K55_nor_df
all_data <- Merge_dimsum_df(stab,RAF,K27,PI3,RAL,SOS,K55)
all_data_pos<-Pos_id(all_data,wt_aa)
######
#figure 1
######
#1e
Plot_heatmap_fitness(all_data_pos,"stab",anno_final3)+
  theme(axis.text.y = element_text(family="Courier New",angle=90,size=5, vjust =.5,hjust =.5,margin=margin(0,-0.5,0,0,"mm")))
#1f
Plot_heatmap_fitness(all_data_pos,"RAF",anno_final3)+
  theme(axis.text.y = element_text(family="Courier New",angle=90,size=5, vjust =.5,hjust =.5,margin=margin(0,-0.5,0,0,"mm")))
#1g
Plot_heatmap_ddG(input="./DATA/weights_Folding.txt",
                      assay_name = "folding")
#1h
Plot_heatmap_ddG(input="./DATA/weights_Binding_RAF.txt",
                      assay_name = "RAF1")
#1j
Scatterplot_fitness_bs(input = all_data_pos,
                            assay_sele="RAF",
                            anno = anno_final3)
#1k
Plot_RAF_invitro_cor(input="./DATA/weights_Binding_RAF.txt",
                          assay_name = "RAF1")
#1h
Weighted_mean_ddG_structure(ddG="./DATA/weights_Folding.txt",
                            input_PDB = "./DATA/pdb6vjj.ent",
                            chain_KRAS = "A",
                            output_PDB_file="./results/figure1H_folding_ddG_pdb6vjj.pdb")
#to make figure 1h using chimeraX
######
#figure 2
######
#2A
Plot_scatter_ddGb_ddGf(ddG1="./DATA/weights_Folding.txt",
                            assay1="folding",
                            ddG2="./DATA/weights_Binding_RAF.txt",
                            assay2="RAF1",
                            anno=anno_final3)
#2B
weighted_mean_abs_ddG_fold<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Folding.txt",
                                                      assay_sele = "folding")
weighted_mean_abs_ddG_RAF<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_RAF.txt",
                                                     assay_sele = "RAF")
Plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                               weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                               weighted_meab_abs_ddGb=weighted_mean_abs_ddG_RAF,
                                               anno_input=anno_final3,
                                               bind="RAF",
                                               bind_ligand = "scHAmin_ligand_RAF")+ggtitle("RAF1")
#2C
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                                 input_PDB = "./DATA/pdb6vjj.ent",
                                 chain_KRAS = "A",
                                 output_PDB_file="./DATA/figure2C_RAF_binding_ab_ddG_6vjj.pdb")
#to make figure 2C using chimeraX, key command: color bfactor /A range -2,2
#2D
Plot_binding_interface_mutation_heatmap("./DATA/weights_Binding_RAF.txt",
                                             assay_sele = "RAF1")
######
#figure 3
######
#3a
Plot_allosteric_mutations(ddG="./DATA/weights_Binding_RAF.txt",
                          anno=anno_final3,
                          assay_sele="RAF",
                          rect_input=rects,
                          rect_alpha=rects_alpha)
#3b
RAF_list<-Plot_major_allosteric_site(ddG="./DATA/weights_Binding_RAF.txt",
                                     anno=anno_final3,
                                     assay_name="RAF",
                                     allosteric_list = RAF_list,
                                     rect_input=rects)
RAF_list[["p"]]
#3c
#pdb 6vjj, chimeraX,
#paste0(RAF_list[["Binding interface site"]],collapse=",")
#paste0(RAF_list[["Allosteric GTP pocket site"]],collapse=",")
#paste0(RAF_list[["Major allosteric site"]],collapse=",")
#3d
Maximum_ddG_structure(ddG="./DATA/weights_Binding_RAF.txt",
                      input_PDB = "./DATA/pdb6vjj.ent",
                      chain_KRAS = "A",
                      output_PDB_file="./DATA/figure3D_max_RAF_binding_ddG_6vjj.pdb")
#to make figure3d using chimeraX, key: color bfactor /A range -2,2
#3e
Plot_ddG_betasheet(ddG="./DATA/weights_Binding_RAF.txt",
                   weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                   anno=anno_final3,
                   assay_sele="RAF",
                   rect_input=rects)
#3f
#6vjj and 6oim, using chimeraX, key: merge #1/A to #2/B
#3g
Maximum_ddG_structure(ddG="./DATA/weights_Binding_RAF.txt",
                      input_PDB = "./DATA/pdb6oim.ent",
                      chain_KRAS = "A",
                      output_PDB_file="./DATA/figure3D_max_RAF_binding_ddG_6oim.pdb")
#to make figure3g using chimeraX, key:color bfactor sel range -2,2
######
#figure 4
######
#4a
Plot_figure4A_point_all_ddG_alphabeta(ddG1="./DATA/weights_Folding.txt",
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
                                           rect_alpha = rects_alpha)
#4b
Plot_RAL_invitro_cor_ITC(ddG="./DATA/weights_Binding_RAL.txt",
                              assay_sele = "RAL")
#4c
Plot_RAL_invitro_cor_GDI(ddG="./DATA/weights_Binding_RAL.txt",
                              assay_sele = "RAL")
######
#figure 5
######
#5a
figure5A<-Get_common_interface(anno_input=anno_final3)
figure5A[["p"]]
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
#bottom
#RAF
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                                 input_PDB = "./DATA/pdb6vjj.ent",
                                 chain_KRAS = "A",
                                 output_PDB_file="./DATA/figure5C_RAF_binding_ab_ddG_6vjj.pdb")
#PI3
weighted_mean_abs_ddG_PI3<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_PI3.txt",
                                                     assay_sele = "PI3")
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_PI3,
                                 input_PDB = "./DATA/pdb1he8.ent",
                                 chain_KRAS = "B",
                                 output_PDB_file="./DATA/figure5C_PI3_binding_ab_ddG_1he8.pdb")
#RAL
weighted_mean_abs_ddG_RAL<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_RAL.txt",
                                                     assay_sele = "RAL")
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_RAL,
                                 input_PDB = "./DATA/pdb1lfd.ent",
                                 chain_KRAS = "B",
                                 Pos_correction=200,
                                 output_PDB_file="./DATA/figure5C_RAL_binding_ab_ddG_1lfd.pdb")
#SOS
weighted_mean_abs_ddG_SOS<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_SOS.txt",
                                                     assay_sele = "SOS")
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_SOS,
                                 input_PDB = "./DATA/1nvw.pdb",
                                 chain_KRAS = "R",
                                 output_PDB_file="./DATA/figure5C_SOS_binding_ab_ddG_1nvw_R.pdb")
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_SOS,
                                 input_PDB = "./DATA/figure5C_SOS_binding_ab_ddG_1nvw_R.pdb",
                                 chain_KRAS = "Q",
                                 output_PDB_file="./DATA/figure5C_SOS_binding_ab_ddG_1nvv_RQ.pdb")
#K27
weighted_mean_abs_ddG_K27<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_K27.txt",
                                                     assay_sele = "K27")
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_K27,
                                 input_PDB = "./DATA/pdb5o2s.ent",
                                 chain_KRAS = "A",
                                 output_PDB_file="./DATA/figure5C_K27_binding_ab_ddG_5o2s.pdb")
#K55 
weighted_mean_abs_ddG_K55<-Get_weighted_mean_abs_ddG(ddG="./DATA/weights_Binding_K55.txt",
                                                     assay_sele = "K55")
Weighted_mean_ddG_structure_fromweighteddata(weighted_mean_ddG=weighted_mean_abs_ddG_K55,
                                 input_PDB = "./DATA/pdb5o2t.ent",
                                 chain_KRAS = "A",
                                 output_PDB_file="./DATA/figure5C_K55_binding_ab_ddG_5o2t.pdb")
#to make figure5b, use chimeraX. key:color bfactor sel range -2,2
#5c
#to make figure5c, use chimeraX. key:color sel gray
#binding interface of binding partners
figure5B<-Get_bindingresidue_bindingpartners()
# paste0(figure5B[binding_partner=="RAF"&scHAmin_ligand<5,Pos],collapse=",")
# paste0(figure5B[binding_partner=="PI3"&scHAmin_ligand<5,Pos],collapse=",")
# paste0(figure5B[binding_partner=="RAL"&scHAmin_ligand<5,Pos],collapse=",")
# #SOS1 chain R
# paste0(figure5B[binding_partner=="SOS1"&scHAmin_ligand<5,Pos],collapse=",")
# #SOS1 chain Q
# paste0(figure5B[binding_partner=="SOS2"&scHAmin_ligand<5,Pos],collapse=",")
# paste0(figure5B[binding_partner=="K27"&scHAmin_ligand<5,Pos],collapse=",")
# paste0(figure5B[binding_partner=="K55"&scHAmin_ligand<5,Pos],collapse=",")
#5d
Plot_ddG_effector_bf_heatmap(ddG1="./DATA/weights_Binding_RAF.txt",
                                  assay1="RAF1",
                                  ddG2="./DATA/weights_Binding_PI3.txt",
                                  assay2="PIK3CG",
                                  ddG3="./DATA/weights_Binding_RAL.txt",
                                  assay3="RALGDS",
                                  anno=anno_final3,
                                  bind1=scHAmin_ligand_RAF,
                                  bind2=scHAmin_ligand_PI3,
                                  bind3=scHAmin_ligand_RAL)
######
#figure 6
######
#6a
all_allosteric_site_list<-Plot_major_allosteric_site_combined_threshold(ddG1="./DATA/weights_Binding_RAF.txt",
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
                                                                             rect_input=rects)
all_allosteric_site_list[["p"]]
#6b
Plot_ddG_all_allosteric_heatmap(ddG1="./DATA/weights_Binding_RAF.txt",
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
                                     list6 = all_allosteric_site_list[["K55"]]
)
######
#ED figure 1
######
#S1a, b
Plot_fitness_correlation_blocks(stability_nor_df,"Abundance")
Plot_fitness_correlation_blocks(RAF_nor_df,"RAF1 binding")
#S1c
Plot_kuriyan_cor(kuriyan_list="./DATA/kuriyan_list.RData")
#S1d
Plot_distribution_fitness(input=all_data_pos,assay_name = "RAF",anno_input = anno_final3)
#S1e
Folding_pre_ob_fitness<-Merge_ddGf_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                                            folding="./DATA/weights_Folding.txt",
                                            block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                            block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                            block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                            wt_aa_input=wt_aa)
Plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness,
                         fold_n=1,
                         mochi_parameters = "./DATA/linears_weights_Abundance1.txt",
                         phenotypen = 1,RT=0.001987*(273+30),bin_input = 50)
Plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness,
                         fold_n=1,
                         mochi_parameters = "./DATA/linears_weights_Abundance2.txt",
                         phenotypen = 2,RT=0.001987*(273+30),bin_input = 50)
Plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness,
                         fold_n=1,
                         mochi_parameters = "./DATA/linears_weights_Abundance3.txt",
                         phenotypen = 3,RT=0.001987*(273+30),bin_input = 50)
#S1f
Plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = Folding_pre_ob_fitness,phenotypen = 1)
Plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = Folding_pre_ob_fitness,phenotypen = 2)
Plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = Folding_pre_ob_fitness,phenotypen = 3)
#S1g
RAF_pre_ob_fitness<-Merge_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                                               block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                               block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                               block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
                                               assay_sele="RAF1",
                                               wt_aa_input=wt_aa)
Plot3d_ddGf_ddGb_ob_pre_fitness(binding_input=RAF_pre_ob_fitness,
                                     folding_assay = Abundance1,
                                     binding_assay = Binding1_RAF,
                                     RT=0.001987*(273+30),
                                     mochi_parameters = "./DATA/linears_weights_Binding1_RAF.txt")
Plot3d_ddGf_ddGb_ob_pre_fitness(binding_input=RAF_pre_ob_fitness,
                                     folding_assay = Abundance2,
                                     binding_assay = Binding2_RAF,
                                     RT=0.001987*(273+30),
                                     mochi_parameters = "./DATA/linears_weights_Binding2_RAF.txt")
Plot3d_ddGf_ddGb_ob_pre_fitness(binding_input=RAF_pre_ob_fitness,
                                     folding_assay = Abundance3,
                                     binding_assay = Binding3_RAF,
                                     RT=0.001987*(273+30),
                                     mochi_parameters = "./DATA/linears_weights_Binding3_RAF.txt")

#S1h
Plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF_pre_ob_fitness,phenotypen = 13)
Plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF_pre_ob_fitness,phenotypen = 14)
Plot2d_ddGf_ob_pre_fitness_perblock(pre_nor = RAF_pre_ob_fitness,phenotypen = 15)
######
#ED figure 2
######
#S2a
Plot_beta_strand_propagation_heattmap(ddG="./DATA/weights_Binding_RAF.txt",
                                           assay_sele = "RAF")+theme(strip.text.x = element_text(margin = margin(t = 0.4,0,b = 0.2,0, "mm")),
                                                                     strip.text = element_text(size=5))+
  coord_fixed(ratio=2)
#S2b
Plot_ddG_ss(ddG="./DATA/weights_Binding_RAF.txt",
                 anno=anno_final3,
                 assay_sele="RAF",
                 rect_input=rects,
                 rect_alpha=rects_alpha_all)
#test
Plot_fisher_test(target_group1="b1",
                 other_group2="others",
                 condition1="allosteric",
                 condition2="others",
                 g1c1 = 29,
                 g1c2 = 142,
                 g2c1 = 363-29,
                 g2c2 = 2712-142,
                 x_y = c("ddG", "structure"))
#test secondary structures one by one
#S2c
Plot_allosteric_mut_distance_2colors(ddG="./DATA/weights_Binding_RAF.txt",
                                          weighted_mean_ddG=weighted_mean_abs_ddG_RAF,
                                          anno=anno_final3,
                                          assay_sele="RAF",
                                          rect_input=rects)
#S2defg
Sotorasib_dis<-Get_sotorasib_pocket(input_file = "./DATA/pdb6oim.ent")
ggarrange(
  Plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                            anno=anno_final3,
                                            pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                            pocket = "pocket1",
                                            assay_sele="RAF",
                                            pocket2=Sotorasib_dis),
  Plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                            anno=anno_final3,
                                            pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                            pocket = "pocket2",
                                            assay_sele="RAF",
                                            pocket2=Sotorasib_dis),
  Plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                            anno=anno_final3,
                                            pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                            pocket = "pocket3",
                                            assay_sele="RAF",
                                            pocket2=Sotorasib_dis),
  Plot_allosteric_ddG_dis_pocket_color(ddG="./DATA/weights_Binding_RAF.txt",
                                            anno=anno_final3,
                                            pocket_anno=final_distance_dc_anno_for_anno_mt_5,
                                            pocket = "pocket4",
                                            assay_sele="RAF",
                                            pocket2=Sotorasib_dis),
  ncol=1,align="v"
)
######
#ED figure 3
######
#S3a
Plot_fitness_correlation_blocks(PI3_nor_df,"PIK3CG binding")
Plot_fitness_correlation_blocks(RAL_nor_df,"RALGDS binding")
Plot_fitness_correlation_blocks(SOS_nor_df,"SOS1 binding")
Plot_fitness_correlation_blocks(K27_nor_df,"DARPin K27 binding")
Plot_fitness_correlation_blocks(K55_nor_df,"DARPin K55 binding")
#S3b
Plot2d_ddGb_ob_pre_fitness(prediction="./DATA/predicted_phenotypes_all.txt",
                                block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData",
                                assay_sele = "PIK3CG",
                                wt_aa_input=wt_aa)
######
#ED figure 4
######
#s4a
Plot_figureS4C_ddG_heatmap_1124(ddG1="./DATA/weights_Folding.txt",
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
                                assay7="DARPin K55")
#s4b
ggarrange(Plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                         weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                         weighted_meab_abs_ddGb=weighted_mean_abs_ddG_PI3,
                                                         anno_input=anno_final3,
                                                         bind="PI3",
                                                         bind_ligand = "scHAmin_ligand_PI3")+ggtitle("PIK3CG"),
          Plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                         weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                         weighted_meab_abs_ddGb=weighted_mean_abs_ddG_RAL,
                                                         anno_input=anno_final3,
                                                         bind="RAL",
                                                         bind_ligand = "scHAmin_ligand_RAL")+ggtitle("RALGDS"),
          Plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                         weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                         weighted_meab_abs_ddGb=weighted_mean_abs_ddG_K55,
                                                         anno_input=anno_final3,
                                                         bind="SOS",
                                                         bind_ligand = "scHAmin_ligand_SOS")+ggtitle("SOS1"),
          Plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                         weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                         weighted_meab_abs_ddGb=weighted_mean_abs_ddG_K27,
                                                         anno_input=anno_final3,
                                                         bind="K27",
                                                         bind_ligand = "scHAmin_ligand_K27")+ggtitle("DARPin K27"),
          Plot_ROCAUC_ddG_fitness_binding_interface(fitness_all=all_data_pos,
                                                         weighted_meab_abs_ddGf=weighted_mean_abs_ddG_fold,
                                                         weighted_meab_abs_ddGb=weighted_mean_abs_ddG_K55,
                                                         anno_input=anno_final3,
                                                         bind="K55",
                                                         bind_ligand = "scHAmin_ligand_K55")+ggtitle("DARPin K55"),
          ncol=6
)
######
#ED figure 5
######
Plot_ddG_all_bf_heatmap(ddG1="./DATA/weights_Binding_RAF.txt",
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
                             bind6=scHAmin_ligand_K55
)+
  theme(strip.background = element_rect(colour="white",fill=NULL,linetype = NULL))
######
#ED figure 6
######
combined_threshold<-Get_combined_threshold(ddG1="./Data/weights_Binding_RAF.txt",
                                           assay1="RAF",
                                           ddG2="./Data/weights_Binding_PI3.txt",
                                           assay2="PI3",
                                           ddG3="./Data/weights_Binding_RAL.txt",
                                           assay3="RAL",
                                           ddG4="./Data/weights_Binding_SOS.txt",
                                           assay4="SOS",
                                           ddG5="./Data/weights_Binding_K27.txt",
                                           assay5="K27",
                                           ddG6="./Data/weights_Binding_K55.txt",
                                           assay6="K55",
                                           anno=anno_final3)
Plot_allosteric_mutations_all(ddG="./Data/weights_Binding_RAF.txt",
                                   anno=anno_final3,
                                   assay_sele ="RAF",
                                   threshold=combined_threshold,
                                   rect_input=rects,
                                   rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2)
Plot_allosteric_mutations_all(ddG="./Data/weights_Binding_PI3.txt",
                                   anno=anno_final3,
                                   assay_sele ="PI3",
                                   threshold=combined_threshold,
                                   rect_input=rects,
                                   rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2)
Plot_allosteric_mutations_all(ddG="./Data/weights_Binding_RAL.txt",
                                   anno=anno_final3,
                                   assay_sele ="RAL",
                                   threshold=combined_threshold,
                                   rect_input=rects,
                                   rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2)
Plot_allosteric_mutations_all(ddG="./Data/weights_Binding_SOS.txt",
                                   anno=anno_final3,
                                   assay_sele ="SOS",
                                   threshold=combined_threshold,
                                   rect_input=rects,
                                   rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2)
Plot_allosteric_mutations_all(ddG="./Data/weights_Binding_K27.txt",
                                   anno=anno_final3,
                                   assay_sele ="K27",
                                   threshold=combined_threshold,
                                   rect_input=rects,
                                   rect_alpha=rects_alpha,y_min=-2.5,y_max=2.5)
Plot_allosteric_mutations_all(ddG="./Data/weights_Binding_K55.txt",
                                   anno=anno_final3,
                                   assay_sele ="K55",
                                   threshold=combined_threshold,
                                   rect_input=rects,
                                   rect_alpha=rects_alpha,y_min=-1.5,y_max=3.2)
######
#ED figure 7
######
Plot_allosteric_mutations_enrichment(ddG="./Data/weights_Binding_RAF.txt",
                                          anno=anno_final3,
                                          assay_sele ="RAF",
                                          threshold=combined_threshold)
Plot_allosteric_mutations_enrichment(ddG="./Data/weights_Binding_PI3.txt",
                                     anno=anno_final3,
                                     assay_sele ="PI3",
                                     threshold=combined_threshold)
Plot_allosteric_mutations_enrichment(ddG="./Data/weights_Binding_RAL.txt",
                                     anno=anno_final3,
                                     assay_sele ="RAL",
                                     threshold=combined_threshold)
Plot_allosteric_mutations_enrichment(ddG="./Data/weights_Binding_SOS.txt",
                                     anno=anno_final3,
                                     assay_sele ="SOS",
                                     threshold=combined_threshold)
Plot_allosteric_mutations_enrichment(ddG="./Data/weights_Binding_K27.txt",
                                     anno=anno_final3,
                                     assay_sele ="K27",
                                     threshold=combined_threshold)
Plot_allosteric_mutations_enrichment(ddG="./Data/weights_Binding_K55.txt",
                                     anno=anno_final3,
                                     assay_sele ="K55",
                                     threshold=combined_threshold)

