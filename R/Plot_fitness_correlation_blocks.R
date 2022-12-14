#' A Function to plot fitness correlation between blocks
#' 
#' This function allows you to lot fitness correlation between blocks.
#' @param input dimsum data
#' @param assay_name binding or abundance fitness
#' 
#' @return Nothing
#' @export
Plot_fitness_correlation_blocks<-function(input=input,assay_name=assay_name){
  ggplot2::ggpairs(input,
          columns = c(18,20,22),
          columnLabels = c("replicate 1", "replicate 2", "replicate 3"),
          mapping=ggplot2::aes(color=block),
          lower = list(continuous=function(data,mapping,...){
            p<-ggplot2::ggplot(data=data,mapping=mapping)+
              ggplot2::geom_bin2d(ggplot2::aes(color=block),bins = 100,alpha=0.1) +
              ggplot2::scale_fill_gradient(low="white",high="black")
            p
          }),
          upper = list(continuous = ggplot2::wrap("cor",mapping=ggplot2::aes(color=block), size = 5*0.35)
                       ),
          diag = list(continuous = "blankDiag"))+
    ggplot2::scale_color_manual(values = c(colour_scheme[["red"]],colour_scheme[["green"]],colour_scheme[["blue"]]))+
    ggplot2::ggtitle(assay_name)+
    ggplot2::theme_classic2()+
    ggplot2::theme(text = ggplot2::element_text(size = 5),
          strip.text.x = ggplot2::element_text(size = 5),
          strip.text.y = ggplot2::element_text(size = 5),
          legend.text = ggplot2::element_text(size = 5))
}
