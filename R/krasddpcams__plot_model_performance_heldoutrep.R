
#' krasddpcams__plot_model_performance_heldoutrep
#'
#' Plot model performance on held out replicate.
#'
#' @param input_dt data.table with model free energy estimates (required)
#'
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot_model_performance_heldoutrep <- function(
  input_dt
){
  
  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(observed_fitness = fitness, predicted_fitness = fitness_pred)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness (Rep. 3)") +
    ggplot2::ylab("Predicted fitness (Reps. 1&2)") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete")^2, 2), sep="")),], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()+
    ggplot2::geom_abline(linetype = "dashed")
  d
}