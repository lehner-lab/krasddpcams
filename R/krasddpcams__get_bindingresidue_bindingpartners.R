
#' Get binding residues of binding partners
#' 
#' This function allows you to get binding residues of binding partners.
#' @param input_list PDB file list
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__get_bindingresidue_bindingpartners<-function(
  input_list
  ){
  RAF_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["RAF"]][["pdb"]], input_list[["RAF"]][["chain_query"]], input_list[["RAF"]][["chain_target"]])
  K27_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["K27"]][["pdb"]], input_list[["K27"]][["chain_query"]], input_list[["K27"]][["chain_target"]])
  PI3_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["PI3"]][["pdb"]], input_list[["PI3"]][["chain_query"]], input_list[["PI3"]][["chain_target"]])
  RAL_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["RAL"]][["pdb"]], input_list[["RAL"]][["chain_query"]], input_list[["RAL"]][["chain_target"]])
  K55_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["K55"]][["pdb"]], input_list[["K55"]][["chain_query"]], input_list[["K55"]][["chain_target"]])
  SOS1_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["SOS1"]][["pdb"]], input_list[["SOS1"]][["chain_query"]], input_list[["SOS1"]][["chain_target"]])
  SOS2_distance_bp<-krasddpcams__minimum_interchain_distances_from_PDB(input_list[["SOS2"]][["pdb"]], input_list[["SOS2"]][["chain_query"]], input_list[["SOS2"]][["chain_target"]])

  RAF_distance_bp[,binding_partner:="RAF"]
  K27_distance_bp[,binding_partner:="K27"]
  PI3_distance_bp[,binding_partner:="PI3"]
  RAL_distance_bp[,binding_partner:="RAL"]
  SOS1_distance_bp[,binding_partner:="SOS1"]
  SOS2_distance_bp[,binding_partner:="SOS2"]
  K55_distance_bp[,binding_partner:="K55"]
  all_long_distance<-rbind(RAF_distance_bp,
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