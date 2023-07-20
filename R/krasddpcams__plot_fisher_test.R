
#' fisher's test for decrease folding ddG in the core
#' 
#' This function allows you to plot fisher's test
#' @param target_group1 target_group1
#' @param other_group2 other_group2
#' @param condition1 condition1
#' @param condition2 condition2
#' @param g1c1 group 1 condition 1 sample size
#' @param g1c2 group 1 condition 2 sample size
#' @param g2c1 group 2 condition 1 sample size
#' @param g2c2 group 2 condition 2 sample size
#' @param x_y x axis name and y axis name
#' @param alternative "greater"...
#' @return Nothing
#' @export
#' @import data.table 
krasddpcams__plot_fisher_test<-function(
  target_group1=target_group1,
  other_group2=other_group2,
  condition1=condition1,
  condition2=condition2,
  g1c1 = g1c1,
  g1c2 = g1c2,
  g2c1 = g2c1,
  g2c2 = g2c2,
  x_y = x_y,
  alternative="greater"
  ){
  #target_group1<-substitute(target_group1)
  #target_group2<-substitute(target_group2)
  dat<-data.frame(target_group1=c(g1c1,g1c2),target_group2=c(g2c1,g2c2),row.names = c(condition1, condition2),stringsAsFactors = F)
  colnames(dat)<-c(target_group1,other_group2)
  x <- c()
  for (row in rownames(dat)) {
    for (col in colnames(dat)) {
      x <- rbind(x, matrix(rep(c(row, col), dat[row, col]), ncol = 2, byrow = TRUE))
    }
  }
  df <- as.data.frame(x)
  colnames(df) <- x_y
  df
  test <- fisher.test(table(df),alternative=alternative)
  ggstatsplot::ggbarstats(
    df, structure, ddG,
    results.subtitle = FALSE,
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      test$p.value,
      "Odd ratios = ",
      test$estimate
    )
  )
}