#' A Function to plot fitness correlation with kuriyan dataset
#' 
#' This function allows you to plot fitness correlation with kuriyan dataset.
#' @param kuriyan_list kuriyan_list"kuriyan_list.RData"
#' 
#' @return Nothing
#' @export
Plot_kuriyan_cor<-function(kuriyan_list=kuriyan_list){
  kuriyan_list<-readRDS(kuriyan_list)
  kuriyan_long_df_dc<-Get_dc_long_kuriyan_data(kuriyan_list_data=kuriyan_list)
  kuriyan_all<-Plot_fitness_kuriyan_data(fitness=all_data_pos,
                                                 kuriyan_data=kuriyan_long_df_dc)
  kuriyan_lm<-lm(nor_fitness~KRas173_RBD,kuriyan_all)
  ggplot(kuriyan_all,aes(x=nor_fitness,y=KRas173_RBD))+
    stat_binhex(bins = 50,size=0,color="black") +
    scale_fill_gradient(low="white",high="black",trans="log10") +
    geom_smooth(method=lm,se=T,size=0.3,color=colour_scheme[["red"]]) +
    annotate("text",x=-0.5,y=0.5,,label = paste0("r = ",round(sqrt(summary(kuriyan_lm)$r.squared),3)) ,size=7*0.35)+
    ylab("Mean E score(Bandaru et al. 2017)")+
    xlab("Binding fitness")+theme_classic2()+
    theme(text = element_text(size=7),
          legend.position="right",
          legend.text = element_text(size=7),
          axis.text.x = element_text(size =7,vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          plot.margin=margin(0,0,0,0))+
    coord_fixed()
  
}