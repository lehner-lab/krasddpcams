#' Maximum free energy to structure
#' 
#' This function allows you to set b factor of sturcture file(PDB) as Maximum free energy changes.
#' @param ddG free energy data
#' @param input_PDB PDB file
#' @param chain_KRAS PDB file chain of KRAS
#' @param output_PDB_file export PDB file
#' @return Nothing
#' @export
Maximum_ddG_structure<-function(ddG=ddG,
                                input_PDB = input_PDB,
                                chain_KRAS = chain_KRAS,
                                output_PDB_file=output_PDB_file){
  ddG<-fread(ddG)
  ddG[,Pos_real:=Pos_ref+1]
  ddG[id!="WT",wt_codon:=substr(id,1,1)]
  ddG[id!="WT",mt_codon:=substr(id,nchar(id),nchar(id))]
  ddG[,mt:=paste0(wt_codon,Pos_real,mt_codon)]
  input_structure<-read.pdb(input_PDB)
  output_PDB<-input_structure
  custom_max <- function(x) {if (length(x)>0) max(x) else 0}
  for(i in unique(output_PDB$atom$resno[output_PDB$atom$chain==chain_KRAS])){
    output_PDB$atom$b[output_PDB$atom$resno==i&output_PDB$atom$chain==chain_KRAS]<-0
  }
  for(i in 2:188){
    output_PDB$atom$b[output_PDB$atom$resno==i&output_PDB$atom$chain==chain_KRAS]<-custom_max(ddG[Pos==i,`mean_kcal/mol`])
  }
  write.pdb(output_PDB,file=output_PDB_file)
}