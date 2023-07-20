
#' A Function to plot fitness correlation with kuriyan dataset
#' 
#' This function allows you to plot fitness correlation with kuriyan dataset.
#' @param kuriyan_list kuriyan_list"kuriyan_list.RData"
#' @param colour_scheme colour scheme list
#' @param all_data_pos all data positions
#' 
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_kuriyan_cor<-function(
  kuriyan_list=kuriyan_list,
  colour_scheme,
  all_data_pos
  ){
  kuriyan_list<-readRDS(kuriyan_list)
  kuriyan_long_df_dc<-krasddpcams__get_dc_long_kuriyan_data(kuriyan_list_data=kuriyan_list)
  kuriyan_all<-krasddpcams__plot_fitness_kuriyan_data(fitness=all_data_pos,
                                                 kuriyan_data=kuriyan_long_df_dc)
  kuriyan_lm<-lm(nor_fitness~KRas173_RBD,kuriyan_all)
  ggplot2::ggplot(kuriyan_all,ggplot2::aes(x=nor_fitness,y=KRas173_RBD))+
    ggplot2::stat_binhex(bins = 50,size=0,color="black") +
    ggplot2::scale_fill_gradient(low="white",high="black",trans="log10") +
    ggplot2::geom_hline(yintercept=0)+
    ggplot2::geom_vline(xintercept=0)+
    ggplot2::geom_abline(intercept = 0,slope=1,linetype="dashed")+
    ggplot2::geom_smooth(method=lm,se=T,size=0.3,color=colour_scheme[["red"]]) +
    ggplot2::annotate("text",x=-0.5,y=0.5,,label = paste0("r = ",round(sqrt(summary(kuriyan_lm)$r.squared),3)) ,size=7*0.35)+
    ggplot2::ylab("Mean E score(Bandaru et al. 2017)")+
    ggplot2::xlab("Binding fitness")+ggpubr::theme_classic2()+
    ggplot2::theme(text = ggplot2::element_text(size=7),
          legend.position="right",
          legend.text = ggplot2::element_text(size=7),
          axis.text.x = ggplot2::element_text(size =7,vjust=.5, hjust=.5),
          axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5,margin=ggplot2::margin(0,0,0,0,"mm")),
          legend.key.height= ggplot2::unit(3.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"),
          plot.margin=ggplot2::margin(0,0,0,0))+
    ggplot2::coord_fixed()
  
}