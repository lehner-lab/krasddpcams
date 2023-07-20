
#' A Function to plot3D binding fitness against folding ddG and binding ddG
#' 
#' This function allows you to plot binding fitness against folding ddG and binding ddG.
#' @param binding_input binding_input
#' @param folding_assay free energy of binding(Abundance1)
#' @param binding_assay binding assay's name(Binding1_RAF)
#' @param mochi_parameters mochi_parameters
#' @param RT RT
#' @param colour_scheme colour scheme list
#' 
#' @return Nothing
#' @export
#' @import data.table
krasddpcams__plot3d_ddGf_ddGb_ob_pre_fitness<-function(
  binding_input=binding_input,
  binding_assay=binding_assay,
  folding_assay=folding_assay,
  mochi_parameters=mochi_parameters,
  RT=RT,
  colour_scheme
  ){
  binding_assay<-substitute(binding_assay)
  folding_assay<-substitute(folding_assay)
  bind_ddg1<-binding_input[eval(binding_assay)==1,.(nt_seq,Nham_aa,aa_seq,fold_1_additive_trait1,fitness)]
  fold_ddg1<-binding_input[eval(folding_assay)==1,.(nt_seq,Nham_aa,aa_seq,fold_1_additive_trait0)]
  ddg_merge1<-merge(fold_ddg1,bind_ddg1,by=c("nt_seq","aa_seq","Nham_aa"))

  binding_range1<-bind_ddg1[Nham_aa>0,range(fold_1_additive_trait1)]
  binding_range1
  binding_energy_grid1<-seq(binding_range1[1],
                           binding_range1[2],
                           (binding_range1[2]-binding_range1[1])/15)
  folding_range1<-fold_ddg1[Nham_aa>0,range(fold_1_additive_trait0)]
  folding_energy_grid1<-seq(folding_range1[1],
                           folding_range1[2],
                           (folding_range1[2]-folding_range1[1])/15)
  folding_energy_grid1
  energy_grid_dt<-as.data.table(expand.grid(folding_energy_grid = folding_energy_grid1, binding_energy_grid = binding_energy_grid1))
  energy_grid_dt
  fitness_binding<-krasddpcams__bind_predict_fitness(mochi_parameters=mochi_parameters,
                                                      folding_energy = energy_grid_dt[,folding_energy_grid],
                                                      binding_energy = energy_grid_dt[,binding_energy_grid],RT=RT)
  # fitness_binding
  pred_fitness_dt_binding <- data.table(
    f_dg_pred = rep(energy_grid_dt[,folding_energy_grid], 2),
    b_dg_pred = rep(energy_grid_dt[,binding_energy_grid], 2),
    observed_fitness = fitness_binding,
    mut_order = rep(c(1, 2), each = energy_grid_dt[,.N]))
  pred_fitness_dt_binding
  plot3D::persp3D(x=binding_energy_grid1*RT,
          y=folding_energy_grid1*RT,
          z=matrix(data=pred_fitness_dt_binding[,observed_fitness],nrow=length(binding_energy_grid1),ncol=length(folding_energy_grid1)),
          r=2, shade=0.4, axes=TRUE,scale=TRUE, box=TRUE, nticks=5, ticktype="detailed", colvar=F, col="white", alpha = 0, border=colour_scheme[["red"]], lwd=0.2,
          cex.lab=1,cex.main=1,cex.axis=1,
          xlab="ddG Folding",
          ylab="ddG Binding",
          zlab="Binding(observed)"
  )
  plot3D::scatter3D(x=ddg_merge1[,fold_1_additive_trait1*RT],
            y=ddg_merge1[,fold_1_additive_trait0*RT],
            z=ddg_merge1[,fitness],
            add = T, col = "black", alpha = 0.2, cex = 0.2
  )
}