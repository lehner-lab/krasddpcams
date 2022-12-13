#' Get binding residues of binding partners
#' 
#' This function allows you to get binding residues of binding partners.
#' @return Nothing
#' @export
Get_bindingresidue_bindingpartners<-function(){
  raf_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/pdb6vjj.ent",
                                                                    chain_query = "B",chain_target = "A")
  K27_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/pdb5o2s.ent",
                                                                    chain_query = "B",chain_target = "A")
  PI3_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/pdb1he8.ent",
                                                                    chain_query = "A",chain_target = "B")
  RAL_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/pdb1lfd.ent",
                                                                    chain_query = "A",chain_target = "B")
  K55_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/pdb5o2t.ent",
                                                                    chain_query = "B",chain_target = "A")
  SOS1_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/1nvw.pdb",
                                                                     chain_query = "S",chain_target = "R")
  SOS2_distance_bp<-doubledeepms__minimum_interchain_distances_from_PDB("/Users/cweng/Desktop/Ras_structure/1nvw.pdb",
                                                                     chain_query = "S",chain_target = "Q")
  raf_distance_bp[,binding_partner:="RAF"]
  K27_distance_bp[,binding_partner:="K27"]
  PI3_distance_bp[,binding_partner:="PI3"]
  RAL_distance_bp[,binding_partner:="RAL"]
  SOS1_distance_bp[,binding_partner:="SOS1"]
  SOS2_distance_bp[,binding_partner:="SOS2"]
  K55_distance_bp[,binding_partner:="K55"]
  all_long_distance<-rbind(raf_distance_bp,
                           K27_distance_bp,
                           PI3_distance_bp,
                           RAL_distance_bp,
                           K55_distance_bp,
                           SOS1_distance_bp,
                           SOS2_distance_bp)
  all_long_distance<-rbind(all_long_distance,data.table(Pos=188,HAmin_ligand=NA,scHAmin_ligand=NA,binding_partner="RAF"))
  all_long_distance<-within(all_long_distance, binding_partner <- factor(binding_partner, levels = c("RAF","PI3","RAL","SOS1","SOS2","K55","K27")))
  all_long_distance
}