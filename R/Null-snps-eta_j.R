#function to calculate z-statistics for Zgamma_hat
Zgamma=function(cor.Z.null.Z,
                cor.Z.null,
                gamma_hat, denominator.sigma){
cov.x.gamma_hat.a=c()
cov.x.gamma_hat=function(x){
return(sum(gamma_hat * cor.Z.null.Z[x,]))
}

cov.x.gamma_hat.a=apply(as.matrix(1:nrow(cor.Z.null.Z)),1, cov.x.gamma_hat)


denominator=function(x){
return(sqrt(x* denominator.sigma))
}

denominator.result=apply(as.matrix(diag(cor.Z.null)), 1, denominator)
return(cov.x.gamma_hat.a/denominator.result)
}

#main funciton to calculate the eta_j
eta_j_fun=function(j,cor.Z.null.Z,
                      cor.Z, cor.Z.snps.set, cor.Z.null.Z.snps,
                      cor.Z.null, num.protein,
                      protein.clump.null, protein.clump.null.Y,
                      gamma_hat,beta_hat, alpha_hat){

library(kernlab)
 
    covfn=function(x,y){
      return(x*y)
    }
    class(covfn)="kernel"
denominator.sigma=c()

Sigma=kernelMatrix(covfn,gamma_hat[[j]])
denominator.sigma=sqrt(sum(cor.Z[[j]]* Sigma))



A_j=Zgamma(cor.Z.null.Z =cor.Z.null.Z[[j]],
          cor.Z.null= cor.Z.null,gamma_hat=gamma_hat[[j]],denominator.sigma= denominator.sigma)

D=list()
for (k in 1: num.protein){
D[[k]]= protein.clump.null[[k]]$z.exposure
}
Y=protein.clump.null.Y$z.exposure



b=c()
c=c()
for (k in 1:num.protein){
b[k]=beta_hat[k]*cor(D[[j]],D[[k]])
c[k]=beta_hat[k]*cor(A_j,D[[k]])*denominator.sigma
}
if(is.null(cor.Z.null.Z.snps)){
eta_j=cor(D[[j]],Y)-sum(b)+sum(c)-denominator.sigma*cor(A_j,Y)
}else{
Sigma_alpha=kernelMatrix(covfn,alpha_hat)
denominator.sigma_alpha=sqrt(sum(cor.Z.snps.set* Sigma_alpha))
Z_alpha=Zgamma(cor.Z.null.Z =cor.Z.null.Z.snps,
          cor.Z.null= cor.Z.null,gamma_hat=alpha_hat,denominator.sigma= denominator.sigma_alpha)

eta_j=cor(D[[j]],Y)-sum(b)+sum(c)-cor(D[[j]],Z_alpha)*denominator.sigma_alpha-denominator.sigma*cor(A_j,Y)+denominator.sigma*cor(A_j,Z_alpha)*denominator.sigma_alpha
}

return(eta_j)
}


