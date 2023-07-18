#' write MAVEdb from DiMSum output fitness
#' 
#' This function allows you to write MAVEdb data.
#' @param input DiMSum output
#' @param assay_name assay_name
#' @param output output csv
#' @return Nothing
#' @export
Write_MAVEdb_scoretable<-function(input=input,
                                  assay_sele=assay_sele,
                                  output=output){
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
                  "V" = "Val",
                  "*" = "Ter")
  convert_aa <- function(sequence) {
    full_names <- aa_mapping[sequence]
    return(full_names)
  }
  output<-input[assay==assay_sele,]
  output[,wtcodon1_3cha:=convert_aa(wtcodon1)]
  output[,wtcodon2_3cha:=convert_aa(wtcodon2)]
  output[,codon1_3cha:=convert_aa(codon1)]
  output[,codon2_3cha:=convert_aa(codon2)]
  output<-output[Nham_aa>0,]
  output[Nham_aa==1,hgvs_pro:=paste0("p.",wtcodon1_3cha,AA_Pos1,codon1_3cha)]
  output[Nham_aa==2,hgvs_pro:=paste0("p.[",wtcodon1_3cha,AA_Pos1,codon1_3cha,";",wtcodon2_3cha,AA_Pos2,codon2_3cha,"]")]
  output[,hgvs_nt:=NA]
  output[,hgvs_splice:=NA]
  output[,score:=nor_fitness]
  output[,SE:=nor_fitness_sigma]
  output<-output[,.(hgvs_nt,hgvs_splice,hgvs_pro,score,SE,block)]
  colnames(output)<-c("hgvs_nt","hgvs_splice","hgvs_pro","score","SE","block")
  return(output)
}