
#' write table s5
#' 
#' This function allows you to write table s5 as csv.
#' @param ddG1 free energy data1
#' @param assay1 data1's name
#' @param ddG2 free energy data2
#' @param assay2 data2's name
#' @param ddG3 free energy data3
#' @param assay3 data3's name
#' @param ddG4 free energy data4
#' @param assay4 data4's name
#' @param ddG5 free energy data5
#' @param assay5 data5's name
#' @param ddG6 free energy data6
#' @param assay6 data6's name
#' @param ddG7 free energy data7
#' @param assay7 data7's name
#' @param ddG8 free energy data8
#' @param assay8 data8's name
#' @param output output file name ".csv"
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__write_tableS5<-function(
  ddG1=ddG1,
  assay1=assay1,
  ddG2=ddG2,
  assay2=assay2,
  ddG3=ddG3,
  assay3=assay3,
  ddG4=ddG4,
  assay4=assay4,
  ddG5=ddG5,
  assay5=assay5,
  ddG6=ddG6,
  assay6=assay6,
  ddG7=ddG7,
  assay7=assay7,
  ddG8=ddG8,
  assay8=assay8,
  output = output
  ){
  Get_fdr<-function(ddG=ddG){
    ddG[!is.na(Pos_real),ddG_fdr:=p.adjust(krasddpcams__pvalue(`mean_kcal/mol`,`std_kcal/mol`),method="BH")]
    ddG[is.na(Pos_real),ddG_conf:=NA]
    ddG[!is.na(Pos_real),ddG_conf:=FALSE]
    ddG[ddG_fdr<0.05,ddG_conf:=TRUE]
    ddG
  }
  ddG1<-krasddpcams__read_ddG(ddG1,assay1)
  ddG1<-Get_fdr(ddG=ddG1)
  ddG2<-krasddpcams__read_ddG(ddG2,assay2)
  ddG2<-Get_fdr(ddG=ddG2)
  ddG3<-krasddpcams__read_ddG(ddG3,assay3)
  ddG3<-Get_fdr(ddG=ddG3)
  ddG4<-krasddpcams__read_ddG(ddG4,assay4)
  ddG4<-Get_fdr(ddG=ddG4)
  ddG5<-krasddpcams__read_ddG(ddG5,assay5)
  ddG5<-Get_fdr(ddG=ddG5)
  ddG6<-krasddpcams__read_ddG(ddG6,assay6)
  ddG6<-Get_fdr(ddG=ddG6)
  ddG7<-krasddpcams__read_ddG(ddG7,assay7)
  ddG7<-Get_fdr(ddG=ddG7)
  ddG8<-krasddpcams__read_ddG(ddG8,assay8)
  ddG8<-Get_fdr(ddG=ddG8)
  
  all_ddG<-rbind(ddG1,ddG2,ddG3,ddG4,ddG5,ddG6,ddG7,ddG8)
  all_ddG[,id:=mt]
  all_ddG[mt=="NANANA",id:="WT"]
  openxlsx::write.xlsx(all_ddG[,.(Pos_real,wt_codon,mt_codon,`mean_kcal/mol`,`std_kcal/mol`,`ci95_kcal/mol`,id,assay,ddG_conf)],
             output,keepNA=TRUE,na.string="NA")
}
