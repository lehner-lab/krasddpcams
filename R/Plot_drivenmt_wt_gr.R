#' Plot growth rate correlation between wt background and mt background
#' 
#' This function allows you to plot fitness correlation between wt background and mt background.
#' @param input dimsum output
#' @param driven_mt driven mutation
#' @return Nothing
#' @export
Plot_drivenmt_wt_gr<-function(input = input,
                              driven_mt = driven_mt,
                              colour_r = colour_r){
  input_WT<-input[Nham_aa<=1,]
  input_WT[,mt_sec:=mt1]
  input_WT[Nham_aa==0,mt_sec:="base"]
  input_MT<-input[Nham_aa>=1&(mt1==driven_mt|mt2==driven_mt),]
  input_MT[mt2==driven_mt,mt_sec:=mt1]
  input_MT[mt1==driven_mt,mt_sec:=mt2]
  input_MT[mt1==driven_mt&Nham_aa==1,mt_sec:="base"]
  input_MT[mt2==driven_mt&Nham_aa==1,mt_sec:="base"]
  input_melt<-input_WT[input_MT,on=.(mt_sec),nomatch = NULL]
  r_model<-lm(nor_gr~i.nor_gr,input_melt)
  ggplot(input_melt[AA_Pos1<65,],aes(x=nor_gr,y=i.nor_gr))+
    geom_point(size=0.1,color=colour_r)+
    annotate("text",x=0.05,y=0.2,label = paste0("r = ",round(sqrt(summary(r_model)$r.squared),2)) ,size=7*0.35,color=colour_r)+
    geom_smooth(method=lm, se=T,size=0.1,color=colour_r)+
    xlab("WT")+
    ylab(driven_mt)+
    geom_vline(xintercept = input_WT[WT==TRUE,nor_gr])+
    geom_hline(yintercept = input_WT[mt1==driven_mt,nor_gr])+
    theme_classic2()+
    theme(text = element_text(size = 5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5))+
    #scale_x_continuous(limits = c(-1.5, highest_gr)) +
    #scale_y_continuous(limits = c(-1.5, highest_gr)) +
    geom_abline(linetype="dashed")+
    coord_fixed()
}