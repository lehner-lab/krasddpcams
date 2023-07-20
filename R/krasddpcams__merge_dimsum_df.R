
#' A Function to merge dimsum output
#' 
#' This function allows you to merge dimsum output.
#' @param merge_1 1st dimsum output
#' @param merge_2 2nd dimsum output
#' @param merge_3 3rd dimsum output
#' @param merge_4 4th dimsum output
#' @param merge_5 5th dimsum output
#' @param merge_6 6th dimsum output
#' @param merge_7 7th dimsum output
#' 
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__merge_dimsum_df <- function(
  merge_1 = merge_1, 
  merge_2 = merge_2, 
  merge_3 = merge_3,
  merge_4 = merge_4, 
  merge_5 = merge_5,
  merge_6 = merge_6,
  merge_7 = merge_7
  ){
  a1 <- as.character(substitute(merge_1))
  a2 <- as.character(substitute(merge_2))
  a3 <- as.character(substitute(merge_3))
  a4 <- as.character(substitute(merge_4))
  a5 <- as.character(substitute(merge_5))
  a6 <- as.character(substitute(merge_6))
  a7 <- as.character(substitute(merge_7))
  merge_1[,assay:=a1]
  merge_2[,assay:=a2]
  merge_3[,assay:=a3]
  merge_4[,assay:=a4]
  merge_5[,assay:=a5]
  merge_6[,assay:=a6]
  merge_7[,assay:=a7]
  output <- rbind(merge_1,merge_2,merge_3,merge_4,merge_5,merge_6,merge_7)
  return(output)
}