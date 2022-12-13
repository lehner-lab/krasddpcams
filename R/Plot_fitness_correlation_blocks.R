#' A Function to plot fitness correlation between blocks
#' 
#' This function allows you to lot fitness correlation between blocks.
#' @param input dimsum data
#' @param assay_name binding or abundance fitness
#' 
#' @return Nothing
#' @export
Plot_fitness_correlation_blocks<-function(input=input,assay_name=assay_name){
  ggpairs(input,
          columns = c(18,20,22),
          columnLabels = c("replicate 1", "replicate 2", "replicate 3"),
          mapping=aes(color=block),
          lower = list(continuous=function(data,mapping,...){
            p<-ggplot(data=data,mapping=mapping)+
              geom_bin2d(aes(color=block),bins = 100,alpha=0.1) +
              scale_fill_gradient(low="white",high="black")
            p
          }),
          upper = list(continuous = wrap("cor",mapping=aes(color=block), size = 5*0.35)
                       ),
          diag = list(continuous = "blankDiag"))+
    scale_color_manual(values = c(colour_scheme[["red"]],colour_scheme[["green"]],colour_scheme[["blue"]]))+
    ggtitle(assay_name)+
    theme_classic2()+
    theme(text = element_text(size = 5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5))
}
