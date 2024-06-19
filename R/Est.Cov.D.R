library(kernlab)
##n1                    sample size of stage 1
#num.protein            number of exposures
##cor.Z.null.Z1         the list of correlation for the common null snps used and IVs used for each exposure
##cor.Z.sample1         the list of correlation for the IVs used for each exposure
##cor.Z.null1           the correlation of null snps
##gamma.hat             stage 1 of weights into list
##protein.clump.null1   the GWAS result for each null snps with "beta.exposure","se.exposure","z.exposure","pval"
Est.Cov.D.fun=function(n1,num.protein, cor.Z.null.Z1,
                   cor.Z.sample1, cor.Z.null1,
                   protein.clump.null1,gamma.hat){
 sigma_jk.matrix=function(j,k){
      return(sigma_jk_fun(n1=n1,cor.Z.null.Z1= cor.Z.null.Z1[[j]], cor.Z.null.Z2= cor.Z.null.Z1[[k]],
                      cor.Z1=cor.Z.sample1[[j]], cor.Z2=cor.Z.sample1[[k]],
                      cor.Z.null= cor.Z.null1,protein.clump1.null=protein.clump.null1[[j]],
                      protein.clump2.null=protein.clump.null1[[k]],gamma_hat1= gamma.hat[[j]],gamma_hat2= gamma.hat[[k]]))
}
 class(sigma_jk.matrix)="kernel"

Estimated.Cov.D =kernelMatrix(sigma_jk.matrix,1:num.protein)
return (Estimated.Cov.D)
}
