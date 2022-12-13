#' import ddG file
#' 
#' This function allows you to get ddG file.
#' @param ddG free energy data
#' @param assay_sele assay_sele
#' 
#' @return ddG data
#' @export
Read_ddG<-function(ddG=ddG,
                   assay_sele=assay_sele){
  ddG<-fread(ddG)
  ddG[,Pos_real:=Pos_ref+1]
  ddG[id!="WT",wt_codon:=substr(id,1,1)]
  ddG[id!="WT",mt_codon:=substr(id,nchar(id),nchar(id))]
  ddG[,mt:=paste0(wt_codon,Pos_real,mt_codon)]
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:188,each=20),
                           mt_codon = unlist(aa_list))
  ddG<-merge(ddG,heatmap_tool,by=c("Pos_real","wt_codon","mt_codon"),all=T)
  ddG[,assay:=assay_sele]
  ddG
}