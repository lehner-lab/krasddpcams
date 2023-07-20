
#' To get sotorasib pocket residues
#' 
#' This function allows you to get sotorasib pocket residues.
#' @param intput_file drug pdb file 6oim
#' @return sortorasib distance
#' @export
krasddpcams__get_sotorasib_pocket<-function(
  input_file=input_file
  ){
  sortorasib_dis<-krasddpcams__minimum_drug_distances_from_PDB(input_file = input_file)
  sortorasib_dis
}