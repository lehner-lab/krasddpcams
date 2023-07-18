#' Plot fitness correlation across triplicates for block1 only
#' 
#' This function allows you to plot fitness correlation across triplicates for block1 only。
#' @param fitness dimsum output
#' @param assay_name assay_name
#' @return Nothing
#' @export
Plot_fitness_correlation_block1<-function(fitness=fitness){
  # filter NA
  fitness_all<-fitness[!is.na(fitness1_uncorr)&!is.na(fitness2_uncorr)&!is.na(fitness3_uncorr),]
  #plot
  ggpairs(fitness_all,
          columns = c(17,19,21),
          columnLabels = c("replicate 1", "replicate 2", "replicate 3"),
          lower = list(continuous=function(data,mapping,...){
            p<-ggplot(data=data,mapping=mapping)+
              geom_bin2d(bins=20)+
              scale_fill_gradient(low="white",high="black",trans = "log")
            p
          }),
          diag = list(continuous = "blankDiag"))+
    theme_classic2()+
    theme(text = element_text(size = 5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5))
}

