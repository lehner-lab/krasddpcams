#' A Function to plot binding fitness against kuriyan data
#' 
#' This function allows you to plot binding fitness against kuriyan data
#' @param fitness fitness data
#' @param kuriyan_data kuriyan_data
#' 
#' @return merged data binding fitness and kuriyan score
#' @export
Plot_fitness_kuriyan_data<-function(fitness=fitness,
                                    kuriyan_data=kuriyan_data){
  # kuriyan_assay<-substitute(kuriyan_assay)
  fitness_single<-fitness[Nham_aa==1&assay=="RAF",]
  escore<-kuriyan_data
  merge_dt<-merge(fitness_single,escore,by.x=c("mt1"),by.y=c("mt"),all=T)
  # cor_fit<-cor(merge_dt[,nor_fitness],
  #              merge_dt[,eval(kuriyan_assay)],method="pearson",use="pairwise.complete.obs")
  merge_dt
  # ggplot(merge_dt,aes(y=nor_fitness,x=KRas173_RBD))+
  #   geom_bin2d(bins = 200)+
  #   scale_fill_gradient(low="white",high="black",trans="log")+
  #   geom_smooth(method="lm") +
  #   stat_cor(method = "pearson",aes(label=..rr.label..), label.x=-1)+
  #   xlab("Mean E score")+
  #   ylab("Binding fitness")+theme_bw()
}