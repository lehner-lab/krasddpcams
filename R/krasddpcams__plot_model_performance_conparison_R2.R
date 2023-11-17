#' A Function to plot model performance
#' 
#' This function allows you to plot model performance
#' @param phenotype a list of phenotype RAF1, folding,,,,
#' @param model_type model type old/new in text file name
#' 
#' @return nothing
#' @export
#' @import data.table
krasddpcams__plot_model_performance_comparison_R2<-function(phenotype = phenotype,
                                               model_type = model_type){
  read_and_merge_data <- function(phenotype, model_type) {
    # data.table
    combined_data <- data.table()
    for (j in phenotype) {
      # read and merge
      for (i in model_type) {
        file_name <- sprintf("./RESULTS/model_%s_%s.txt", j, i)
        
        # test file address
        if (file.exists(file_name)) {
          # read and set model and phenotypes
          current_data <- fread(file_name, header = TRUE)
          current_data[,model:=i]
          current_data[,phenotypes:=j]
          # merge data.table
          combined_data <- rbindlist(list(combined_data, current_data), fill = TRUE)
        } else {
          warning(sprintf("File %s not found.", file_name))
        }
      }
    }
    return(combined_data)
  }
  merged_data <- read_and_merge_data(phenotype, model_type)
  old_model1_name<-c("training","^test$","orig","^randomtest$","RAF1","PIK3CG","RALGDS","SOS1","DARPinK27","DARPinK55","folding")
  new_model1_name<-c("Training set","Test set","Random doubles","Random singles & doubles","KRAS-RAF1 BindingPCA","KRAS-PIK3CG BindingPCA","KRAS-RALGDS BindingPCA","KRAS-SOS1 BindingPCA","KRAS-DARPin K27 BindingPCA","KRAS-DARPin K55 BindingPCA","KRAS AbundancePCA")
  for (i in seq_along(old_model1_name)) {
    merged_data[,model:=gsub(old_model1_name[i],new_model1_name[i],model)]
    merged_data[,dataset:=gsub(old_model1_name[i],new_model1_name[i],dataset)]
    merged_data[,phenotypes:=gsub(old_model1_name[i],new_model1_name[i],phenotypes)]
  }
  merged_data<-within(merged_data,
                      dataset<-factor(dataset,
                                      levels=c("Training set",
                                               "Test set")))
  merged_data<-within(merged_data,
                      model<-factor(model,
                                    levels=c("Random doubles",
                                             "Random singles & doubles"
                                    )
                      ))
  merged_data<-within(merged_data,
                      phenotypes<-factor(phenotypes,
                                         levels=c("KRAS-RAF1 BindingPCA","KRAS-PIK3CG BindingPCA","KRAS-RALGDS BindingPCA","KRAS-SOS1 BindingPCA","KRAS-DARPin K27 BindingPCA","KRAS-DARPin K55 BindingPCA","KRAS AbundancePCA")
                      ))
  plot_fill <- c(colour_scheme[["red"]],colour_scheme[["blue"]])
  names(plot_fill) <- c("Random singles & doubles", "Random doubles")
  plot_alpha <- c(1, 0.8)
  names(plot_alpha) <- c("Training set", "Test set")
  ggplot2::ggplot(merged_data,ggplot2::aes(x=model,fill=model,y=R_squared, alpha=dataset))+
    ggplot2::geom_col(width=0.5, position = ggplot2::position_dodge(0.8))+
    ggplot2::facet_wrap(~ phenotypes)+
    ggplot2::coord_fixed()+
    ggplot2::coord_cartesian(ylim = c(0.6,1))+
    ggplot2::theme_classic()+
    ggplot2::ylab("Model performance (R-squared)")+
    ggplot2::scale_fill_manual(values=plot_fill)+
    ggplot2::scale_alpha_manual(values=plot_alpha)+
    ggplot2::theme(text = ggplot2::element_text(size=7),
          axis.text = ggplot2::element_text(size=7),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust=0.5),
          legend.text = ggplot2::element_text(size=7),
          plot.title = ggplot2::element_text(size=7))
}