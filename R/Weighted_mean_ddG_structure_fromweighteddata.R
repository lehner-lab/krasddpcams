#' weighted mean free energy to structure
#' 
#' This function allows you to set b factor of sturcture file(PDB) as weighted mean free energy changes.
#' @param weighted_mean_ddG free energy data
#' @param input_PDB PDB file
#' @param chain_KRAS PDB file chain of KRAS
#' @param output_PDB_file export PDB file
#' @return Nothing
#' @export 
Weighted_mean_ddG_structure_fromweighteddata<-function(weighted_mean_ddG=weighted_mean_ddG,
                                      input_PDB=input_PDB,
                                      chain_KRAS=chain_KRAS,
                                      Pos_correction=0,
                                      output_PDB_file=output_PDB_file){
  weighted_mean_ddG_sel<-weighted_mean_ddG
  input_structure<-read.pdb(input_PDB)
  output_PDB<-input_structure
  for(i in unique(output_PDB$atom$resno[output_PDB$atom$chain==chain_KRAS])){
    output_PDB$atom$b[output_PDB$atom$resno==i&output_PDB$atom$chain==chain_KRAS]<-0
  }
  for(i in 2:166){
    if(!is.na(weighted_mean_ddG_sel[Pos_real==i,mean])){
      output_PDB$atom$b[output_PDB$atom$resno==(i+Pos_correction)&output_PDB$atom$chain==chain_KRAS]<-weighted_mean_ddG_sel[Pos_real==i,mean]
    }
  }
  write.pdb(output_PDB,file=output_PDB_file)
}
