
#' Plot common binding interface residues
#' 
#' This function allows you to plot common binding interface residues.
#' @param anno_input anno_input final_anno3
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__get_common_interface<-function(
  anno_input=anno_input
  ){
  matrix_common<-anno_input[binding_RAF%in%c("RAF","both")|
                              binding_RAL%in%c("RAL","both")|
                              binding_PI3%in%c("PI3","both")|
                              binding_SOS%in%c("SOS","both")|
                              binding_K27%in%c("K27","both")|
                              binding_K55%in%c("K55","both")
                            ,.(Pos,WT_AA,binding_RAF,binding_RAL,binding_PI3,binding_SOS,binding_K27,binding_K55)]
  matrix_common[,binding_common:=F]
  matrix_common[binding_RAF%in%c("RAF","both")&
                  binding_RAL%in%c("RAL","both")&
                  binding_PI3%in%c("PI3","both")&
                  binding_SOS%in%c("SOS","both")&
                  binding_K27%in%c("K27","both")&
                  binding_K55%in%c("K55","both"),binding_common:=T]
  #RAF
  matrix_common[,binding_RAF_common:="no"]
  matrix_common[binding_RAF%in%c("RAF","both"),binding_RAF_common:="contact"]
  matrix_common[binding_common==T,binding_RAF_common:="common"]
  #RAL
  matrix_common[,binding_RAL_common:="no"]
  matrix_common[binding_RAL%in%c("RAL","both"),binding_RAL_common:="contact"]
  matrix_common[binding_common==T,binding_RAL_common:="common"]
  #PI3
  matrix_common[,binding_PI3_common:="no"]
  matrix_common[binding_PI3%in%c("PI3","both"),binding_PI3_common:="contact"]
  matrix_common[binding_common==T,binding_PI3_common:="common"]
  #SOS
  matrix_common[,binding_SOS_common:="no"]
  matrix_common[binding_SOS%in%c("SOS","both"),binding_SOS_common:="contact"]
  matrix_common[binding_common==T,binding_SOS_common:="common"]
  #K27
  matrix_common[,binding_K27_common:="no"]
  matrix_common[binding_K27%in%c("K27","both"),binding_K27_common:="contact"]
  matrix_common[binding_common==T,binding_K27_common:="common"]
  #K55
  matrix_common[,binding_K55_common:="no"]
  matrix_common[binding_K55%in%c("K55","both"),binding_K55_common:="contact"]
  matrix_common[binding_common==T,binding_K55_common:="common"]
  
  
  matrix_common_long_RAF<-data.table(wt_pos=matrix_common[,paste0(WT_AA,Pos)],
                                     binding=matrix_common[,binding_RAF_common],
                                     Pos=matrix_common[,Pos],
                                     assay="RAF1")
  matrix_common_long_RAL<-data.table(wt_pos=matrix_common[,paste0(WT_AA,Pos)],
                                     binding=matrix_common[,binding_RAL_common],
                                     Pos=matrix_common[,Pos],
                                     assay="RALGDS")
  matrix_common_long_PI3<-data.table(wt_pos=matrix_common[,paste0(WT_AA,Pos)],
                                     binding=matrix_common[,binding_PI3_common],
                                     Pos=matrix_common[,Pos],
                                     assay="PIK3CG")
  matrix_common_long_SOS<-data.table(wt_pos=matrix_common[,paste0(WT_AA,Pos)],
                                     binding=matrix_common[,binding_SOS_common],
                                     Pos=matrix_common[,Pos],
                                     assay="SOS1")
  matrix_common_long_K27<-data.table(wt_pos=matrix_common[,paste0(WT_AA,Pos)],
                                     binding=matrix_common[,binding_K27_common],
                                     Pos=matrix_common[,Pos],
                                     assay="DARPin K27")
  matrix_common_long_K55<-data.table(wt_pos=matrix_common[,paste0(WT_AA,Pos)],
                                     binding=matrix_common[,binding_K55_common],
                                     Pos=matrix_common[,Pos],
                                     assay="DARPin K55")
  matrix_common_long_all<-rbind(matrix_common_long_RAF,
                                matrix_common_long_RAL,
                                matrix_common_long_PI3,
                                matrix_common_long_SOS,
                                matrix_common_long_K27,
                                matrix_common_long_K55)
  matrix_common_long_all<-within(matrix_common_long_all,
                                 wt_pos<-factor(wt_pos,
                                                levels=unique(matrix_common_long_all[order(Pos),wt_pos])))
  # matrix_common_long_all<-within(matrix_common_long_all,
  #                                binding<-factor(binding,
  #                                               levels=c("common","contact","no")))
  output<-list()
  output[["p"]]<-ggplot2::ggplot(matrix_common_long_all,ggplot2::aes(x=wt_pos,y=assay,fill=binding))+
    ggplot2::geom_tile(color="black")+
    ggplot2::scale_fill_manual(values=c("orange","yellow","white"))+
    ggpubr::theme_classic2()+
    ggplot2::coord_fixed()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7,angle = 90, vjust = .5,hjust = 0.5),
          axis.text.y = ggplot2::element_text(size = 7,vjust = .5,hjust = .5),
          text = ggplot2::element_text(size=7),
          legend.text = ggplot2::element_text(size=7),
          axis.ticks.x=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position="right",
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width= ggplot2::unit(3.1, 'mm'))+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)
  output[["matrix"]]<-matrix_common_long_all
  output
}
