sigma_jk_fun=function(n1, cor.Z.null.Z1, cor.Z.null.Z2,
                      cor.Z1, cor.Z2,
                      cor.Z.null, 
                      protein.clump1.null,
                      protein.clump2.null, gamma_hat1, gamma_hat2){

#cor.Z.null.Z1=cor(Z.null[,j],Z1)

####for Z_gamma1

#calculate the sd of Zgamma
library(kernlab)
 
    covfn=function(x,y){
      return(x*y)
    }
    class(covfn)="kernel"

Sigma1=kernelMatrix(covfn,gamma_hat1)
denominator.sigma1=sqrt(sum(cor.Z1* Sigma1))

Sigma2=kernelMatrix(covfn,gamma_hat2)
denominator.sigma2=sqrt(sum(cor.Z2* Sigma2))



############################
###for Zgamma_1 Z statistics
############################


Zgamma=function(cor.Z.null.Z,
                cor.Z,cor.Z.null,
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


A1=Zgamma(cor.Z.null.Z =cor.Z.null.Z1,
          cor.Z.null= cor.Z.null,gamma_hat=gamma_hat1,denominator.sigma= denominator.sigma1 )
A2=Zgamma(cor.Z.null.Z =cor.Z.null.Z2,
          cor.Z.null= cor.Z.null,gamma_hat=gamma_hat2,denominator.sigma= denominator.sigma2)





D1=protein.clump1.null$z.exposure
D2=protein.clump2.null$z.exposure
sigma_jk=cor(D1,D2)-cor(D1,A2)* denominator.sigma2-cor(A1,D2)*denominator.sigma1 +cor(A1,A2)* denominator.sigma1* denominator.sigma2

return(sigma_jk)
}
