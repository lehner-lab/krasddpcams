
#' Plot fitness correlation across triplicates for block1 only
#' 
#' This function allows you to plot fitness correlation across triplicates for block1 only。
#' @param fitness dimsum output
#' @param assay_name assay_name
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_fitness_correlation_block1<-function(
  fitness=fitness
  ){
  # filter NA
  fitness_all<-fitness[!is.na(fitness1_uncorr)&!is.na(fitness2_uncorr)&!is.na(fitness3_uncorr),]
  #plot
  d <- GGally::ggpairs(fitness_all,
          columns = c(17,19,21),
          columnLabels = c("replicate 1", "replicate 2", "replicate 3"),
          lower = list(continuous=function(data,mapping,...){
            p<-ggplot2::ggplot(data=data,mapping=mapping)+
              ggplot2::geom_bin2d(bins=20)+
              ggplot2::scale_fill_gradient(low="white",high="black",trans = "log")
            p
          }),
          diag = list(continuous = "blankDiag"))+
    ggpubr::theme_classic2()+
    ggplot2::theme(text = ggplot2::element_text(size = 5),
          strip.text.x = ggplot2::element_text(size = 5),
          strip.text.y = ggplot2::element_text(size = 5),
          legend.text = ggplot2::element_text(size = 5))
    return(d)
}

