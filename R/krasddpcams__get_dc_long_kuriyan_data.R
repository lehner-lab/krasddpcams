
#' kuriyan dcast data table
#' 
#' This function allows you to get kuriyan dcast data table.
#' @param kuriyan_list_data kuriyan_list_data
#' @return kuriyan_dcast_data
#' @export
#' @import data.table
krasddpcams__get_dc_long_kuriyan_data<-function(
  kuriyan_list_data=kuriyan_list_data
  ){
  kuriyan_long_df<-data.table()
  for (i in names(kuriyan_list_data)){
    kuriyan_long_list<-data.table(kuriyan_list_data[[i]])
    kuriyan_long_list[,assay:=i]
    kuriyan_long_df<-rbind(kuriyan_long_df,kuriyan_long_list)
  }
  kuriyan_long_df[,mt:=paste0(wt_aa,Pos,mt_aa)]
  kuriyan_long_df_dc<-dcast(kuriyan_long_df,mt~assay,value.var = c("mean_E_score"))
  kuriyan_long_df_dc
}