#' A Function to plot fitness correlation between single and double mutations
#' 
#' This function allows you to lot fitness correlation between single and double mutations
#' @param input dimsum data
#' @param mut mut
#' @param highest_fitness highest_fitness
#' 
#' @return Nothing
#' @export
Plot_fitness_correlation_doublevssingle<-function(input=input,mut=mut,highest_fitness=highest_fitness){
  input_WT<-input[Nham_aa<=1,]
  input_WT[,mt_sec:=mt1]
  input_WT[Nham_aa==0,mt_sec:="base"]
  input_MT<-input[Nham_aa>=1&(mt1==mut|mt2==mut),]
  input_MT[mt2==mut,mt_sec:=mt1]
  input_MT[mt1==mut,mt_sec:=mt2]
  input_MT[mt1==mut&Nham_aa==1,mt_sec:="base"]
  input_MT[mt2==mut&Nham_aa==1,mt_sec:="base"]
  input_melt<-input_WT[input_MT,on=.(mt_sec),nomatch = NULL]
  ggplot(input_melt[AA_Pos1<65,],aes(x=nor_fitness,y=i.nor_fitness))+
    geom_point()+
    xlab("WT")+
    ylab(mut)+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = input_WT[mt1==mut,fitness])+
    theme_classic2()+
    theme(text = element_text(size = 5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5))+
    scale_x_continuous(limits = c(-1.5, highest_fitness)) +
    scale_y_continuous(limits = c(-1.5, highest_fitness)) +
    coord_fixed()
}