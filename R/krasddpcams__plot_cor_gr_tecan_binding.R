
#' validation tecan vs selection for bindingPCA
#' 
#' This function allows you to plot correlation between tecan culture growthrate and selection culture fitness.
#' @param tecandata tecandata output
#' @param plasmidid plasmid id
#' @param single single mutation fitness from DiMSum output
#' @param colour_scheme colour scheme list
#' @return Nothing
#' @export
#' @import data.table 
#' @import dplyr
krasddpcams__plot_cor_gr_tecan_binding<-function(
  tecandata=tecandata, 
  plasmidid = plasmidid,
  single=single,
  colour_scheme
  ){
  # read OD data
  a = read.delim(tecandata)
  od = as.matrix(a[,-1])
  
  # get the real OD by subracting the blank
  real_od = do.call("rbind",apply(od, 1, function(x){ x - a[a$well == "O23", -1]}))
  pls2id = read.csv(plasmidid, sep=';')
  names(pls2id)[1] <- "well"
  
  # creat file with real ODs and plasmid name
  df_od = data.frame(real_od)
  df_od$well = a$well
  df_od = df_od[!(df_od$well == "O23"),]
  
  # merge with plasmid info
  final_od = merge(df_od, pls2id, by="well")
  final_od_melt = melt(final_od, id.vars = c("well","plasmid","mutation","background","assay"), variable.name = "time")
  final_od_melt$time = sapply(as.character(final_od_melt$time), function(x){as.numeric(substr(x,2,nchar(x)))/3600})
  
  final_od_melt$value[final_od_melt$value < 0] <- 0
  final_od_melt$value[final_od_melt$value <= 0] <- 0.001
  final_od_melt<-as.data.table(final_od_melt)
  final_od_melt<-final_od_melt[plasmid!="blank",]
  #
  gr_all<-do.call("rbind", lapply(unique(final_od_melt$well), function(x){
    ds = final_od_melt[final_od_melt$well == x & final_od_melt$time <= 70,]
    max_od = max(ds$value)
    mid_od = max_od - (max_od-min(ds$value))/2
    time_odmax = ds$time[ds$value == max_od]
    time_od50 = max(ds$time[ds$value <= mid_od & ds$time <= time_odmax[1]])
    lmds = data.frame(t = ds$time[ds$time <= time_od50+4 & ds$time >= time_od50-4])
    lmds$ods = log(ds$value[ds$time %in% lmds$t])
    lmfit = lm(formula = ods ~ t, data = lmds)
    data.table(plasmid = unique(ds$plasmid), well = x, mutation = unique(ds$mutation),background=unique(ds$background),assay= unique(ds$assay),
               growth_rate_slope = lmfit$coefficients[2]) 
  }))
  gr_summary <- gr_all %>%
    group_by(plasmid) %>%
    summarize(
      mean = mean(growth_rate_slope),
      variance = var(growth_rate_slope),
      sd = sd(growth_rate_slope),
      mutation = unique(mutation),
      background = unique(background),
      assay= unique(assay)
    )
  gr_summary<-as.data.table(gr_summary)
  #
  single[,mutation:=mt1]
  single[Nham_aa==0,mutation:="wt"]
  data_gr<-merge(as.data.table(gr_summary[background=="WT",]),single,by = "mutation")
  data_gr[mutation=="wt",mutation:="wt2"]
  data_gr[Nham_nt==0&block=="block1",mutation:="wt"]
  bind_lm_gr<-lm(nor_fitness~mean,data_gr[assay.x=="binding"&assay.y=="RAF"&mutation!="wt2",])
  ggplot2::ggplot(data_gr[assay.x=="binding"&assay.y=="RAF"&mutation!="wt2",],ggplot2::aes(x=mean,y=nor_fitness,label=mutation))+
    ggplot2::geom_smooth(method=lm, se=T,size=0.1,color=colour_scheme[["red"]])+
    ggplot2::xlab("Growth rate of individual measurements")+
    ggplot2::ylab("fitness of BindingPCA")+
    ggplot2::geom_pointrange(ggplot2::aes(ymin=nor_fitness-nor_fitness_sigma, ymax=nor_fitness+nor_fitness_sigma),size=0.3,color="gray")+
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=mean-sd, xmax=mean+sd,height = 0),size=0.3,color="gray")+
    ggplot2::geom_point(size=2)+
    ggrepel::geom_text_repel(size=7*0.35)+
    ggplot2::annotate("text",x=0.08,y=0.23,label = paste0("r = ",round(sqrt(summary(bind_lm_gr)$r.squared),2)) ,size=7*0.35)+
    ggpubr::theme_classic2()+
    ggplot2::theme(text = ggplot2::element_text(size=7),
          legend.position="right",
          legend.text = ggplot2::element_text(size=7),
          axis.text.x = ggplot2::element_text(size =7,vjust=.5, hjust=1),
          axis.text.y = ggplot2::element_text(size=7, vjust = .5,hjust = .5,margin=ggplot2::margin(0,0,0,0,"mm")),
          legend.key.height= ggplot2::unit(1.1, 'mm'),
          legend.key.width = ggplot2::unit(3.1, 'mm'),
          legend.key.size = ggplot2::unit(1,"mm"),
          plot.margin=ggplot2::margin(0,0,0,0))
}
  
