#' write MAVEdb from Mochi output ddG
#' 
#' This function allows you to write MAVEdb data.
#' @param input MoCHI output
#' @param assay_name assay_name
#' @param output output csv
#' @return Nothing
#' @export
Write_MAVEdb_scoretable_ddG<-function(input=input,assay_name=assay_name,output=output){
  aa_mapping <- c("A" = "Ala",
                  "R" = "Arg",
                  "N" = "Asn",
                  "D" = "Asp",
                  "C" = "Cys",
                  "Q" = "Gln",
                  "E" = "Glu",
                  "G" = "Gly",
                  "H" = "His",
                  "I" = "Ile",
                  "L" = "Leu",
                  "K" = "Lys",
                  "M" = "Met",
                  "F" = "Phe",
                  "P" = "Pro",
                  "S" = "Ser",
                  "T" = "Thr",
                  "W" = "Trp",
                  "Y" = "Tyr",
                  "V" = "Val")
  convert_aa <- function(sequence) {
    full_names <- aa_mapping[sequence]
    return(full_names)
  }
  data<-Read_ddG(ddG = input,assay_sele = assay_name)
  data[,hgvs_nt:=NA]
  data[,hgvs_splice:=NA]
  data[,wt_codon_3cha:=convert_aa(wt_codon)]
  data[,mt_codon_3cha:=convert_aa(mt_codon)]
  data[,hgvs_pro:=paste0("p.",wt_codon_3cha,Pos_real,mt_codon_3cha)]
  data[,score:=`mean_kcal/mol`]
  data[,std:=`std_kcal/mol`]
  output<-data[!is.na(Pos_real)&wt_codon!=mt_codon,.(hgvs_nt,hgvs_splice,hgvs_pro,score,std)]
  return(output)
}