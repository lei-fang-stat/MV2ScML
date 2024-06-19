library(kernlab)
#input
#n1 sample size of stage 1
#n3 sample size of stage 2
#num.protein number of exposures
#num.snp number of unique IVs for all exposures
#protein.clump GWAS summary statistics
#gamma.hat stage 1 of weights into list
#cor.Z.sample3 correlation matrix of in-sample LD for stage 2 or from reference panel with colnames of SNPs id.
#exposure.of.interest a set of interested exposured included for multivariate analysis
#SNP list of IV names
#SNP.comb unique of the combined IV names
#outcome.GWAS.i is the stage2 GWAS statistics for each exposure
#outcome.GWAS.comb is the combined stage2 GWAS statistics for all IVs (must be the same order with SNP.comb)
#K.vec the tunning parameter for BIC selection, default is 0:ceiling (num.snp/2)
#res.D residual for each exposure
#Est.Cov.D estimated covariance for each residuals in stage 1
#snp.set.provide the invalid IVs, default is NULL unless user specified

##output
##the estimated coefficient
##the standard error
##the naive standard error with considering stage 1 variability
##the exposures included in the analysis
##the invalid IVs identified
MV.2ScML.nooverlap=function(n1,n3,num.protein, num.snp,gamma.hat,
                            cor.Z.sample3,SNP,SNP.comb,outcome.GWAS.i,outcome.GWAS.comb,
                            K.vec,Est.Cov.D,snp.set.provide){
  n2=0 # the number of overlapping samples
  Est.Cov.DY = rep(0,length(c(1:num.protein)))
  gamma.hat.ind=list()
  gamma.hat.trans=list()
  for ( j in 1:num.protein){
    gamma.hat.ind[[j]]=match(SNP[[j]][which(gamma.hat[[j]]!=0)],SNP.comb )
    gamma.hat.trans[[j]]=rep(0,length(SNP.comb))
    gamma.hat.trans[[j]][gamma.hat.ind[[j]] ]=gamma.hat[[j]]
  }
  gamma.hat.mat=do.call(cbind, gamma.hat.trans) #make gamma.hat into a matrix form


	##combined snps

      	cor.DjDk.matrix=function(x,y){
         	return(gamma.hat[[x]]%*%cor.Z.sample3[SNP[[x]], SNP[[y]]]%*%gamma.hat[[y]])
      	}
      	class(cor.DjDk.matrix)="kernel"

      	DjhatDkhat=kernelMatrix(cor.DjDk.matrix,c(1:num.protein))


      	vec.DjhatZ.list=list()
      	for (j in 1:num.protein){
       	vec.DjhatZ.list[[j]]=as.numeric(gamma.hat[[j]]%*% cor.Z.sample3[SNP[[j]], SNP.comb])
      	}
     	#get design matrix
      	A11= DjhatDkhat
      	A12= do.call(rbind,vec.DjhatZ.list)
      	A21=t(A12)
      	A22=cor.Z.sample3[SNP.comb,SNP.comb]
      	Cov.Design = rbind(cbind(A11,A12),cbind(A21,A22))

      	Cov.Design= Cov.Design + diag( dim(Cov.Design)[2] )*posdef # num.gene + num.snp
      	Eigen.Cov = eigen(Cov.Design)

      	Eigen.Cov$values[(num.snp+1): (num.protein+num.snp)] = 1e-5 ### add small number to avoid singularity
      	Half.Cov =
        	Eigen.Cov$vectors%*%
        	diag(sqrt(Eigen.Cov$values))%*%
        	t(Eigen.Cov$vectors)
      	inv.Half.Cov = solve(Half.Cov,tol=0)
      	cov_D_Y=c()
      	for (j in 1: num.protein){
       	cov_D_Y[j]=sum(outcome.GWAS.i[[j]]$Beta/sqrt( (outcome.GWAS.i[[j]]$Beta)^2+(n3-2)*(outcome.GWAS.i[[j]]$SE)^2)* gamma.hat[[j]])
      	}

      	cov_Z_Y=c(outcome.GWAS.comb$Beta/sqrt( (outcome.GWAS.comb$Beta)^2+(n3-2)*(outcome.GWAS.comb$SE)^2))


      	cov_DesignX_sim.Y=as.numeric(c(cov_D_Y, cov_Z_Y))

      	inputY = as.numeric(inv.Half.Cov%*% as.matrix(cov_DesignX_sim.Y))
      	inputX = Half.Cov
#####TLC

      	TLC.res1.mat = NULL
      	tlc.weight = rep(1,num.snp+num.protein)
      	tlc.weight[1:num.protein] = 0
      	Est.error.vec1=NULL
	for(K in K.vec)
      	{
        	TLC.res1 = TLC(Y = inputY, X = inputX, K = K,tlc_weight = tlc.weight)
        	TLC.res1.mat = rbind(TLC.res1.mat,TLC.res1)
        	Est.error= as.numeric(1- 2*sum(cov_DesignX_sim.Y*TLC.res1) +
                       TLC.res1%*%Cov.Design%*% TLC.res1)
      	if (Est.error<0){
        	Est.error=1}
        	Est.error.vec1=c(Est.error.vec1,Est.error)
      	}

	BIC.vec1 = (n2+n3)*log(Est.error.vec1) + log((n2+n3))*K.vec
	##use minimum BIC to identify invalid IVs
	nonzero.set = which(TLC.res1.mat[which.min(BIC.vec1),]!=0)
      	#we do not select proteins
	protein.set =c(1:num.protein)
      	snp.set = nonzero.set[nonzero.set> num.protein]-num.protein
        if (length(snp.set.provide)==0){
         snp.set= snp.set
        }else{
        snp.set= snp.set.provide
        }
      	##check whether snp.set=0
      	if (length(snp.set)==0){
      		Cov.Design.BIC=A11
      		cov_DesignX_sim.Y.BIC=cov_D_Y
      	}else{

      		################Obtain Inference after variable selection
      		### estimate covariance of Z
      		#use selected invalid IVs to refit the model
      		SNP.BIC=SNP.comb[snp.set]
      		vec.DjhatZ.BIC.list=list()
        	#recalculate the design matrix
      		for (j in 1:num.protein){
       		vec.DjhatZ.BIC.list[[j]]=as.numeric(gamma.hat[[j]]%*% cor.Z.sample3[SNP[[j]], SNP.BIC])
      		}

      		A12.BIC= do.call(rbind, vec.DjhatZ.BIC.list)
      		A21.BIC=t(A12.BIC)
      		if (length(snp.set)==1){
      		A22.BIC=1
      		}else{
      			A22.BIC=cor.Z.sample3[SNP.BIC,SNP.BIC]
      		}
      		Cov.Design.BIC = rbind(cbind(A11,A12.BIC),cbind(A21.BIC,A22.BIC))

      		outcome.GWAS.comb.BIC=outcome.GWAS.comb[snp.set,]
      		cov_Z_Y.BIC= c(outcome.GWAS.comb.BIC$Beta/sqrt( (outcome.GWAS.comb.BIC$Beta)^2+(n3-2)*(outcome.GWAS.comb.BIC$SE)^2))
      		cov_DesignX_sim.Y.BIC=as.numeric(c(cov_D_Y, cov_Z_Y.BIC))
  	}

  	BetaAlpha.hat = solve(Cov.Design.BIC)%*%cov_DesignX_sim.Y.BIC
  	BetaAlpha.hat

  	BetaHat = BetaAlpha.hat[1:length(protein.set),]

  	AlphaHat = BetaAlpha.hat[-(1:length(protein.set)),]

	### estimates of variance component
	gamma.hat.mat= gamma.hat.mat
	Est.Cov.Z = cor.Z.sample3[SNP.comb, SNP.comb]
	Est.Cov.D= Est.Cov.D
  	Est.Cov = rbind(cbind(Est.Cov.D,Est.Cov.DY),
                  c(Est.Cov.DY,0))
  	##stage2 error variance
  	covX.Y= sum(BetaAlpha.hat*cov_DesignX_sim.Y.BIC)
  	covX= t(BetaAlpha.hat)%*%Cov.Design.BIC%*% BetaAlpha.hat
  	covX= covX+diag( dim(covX)[2] )*posdef
  	Est.Var.total = as.numeric(1-2*covX.Y+covX)

  	### calculate matrix W
  	set.list = list()
  	for(ind in protein.set)
  	{
    	set.list = c(set.list,list(set = which(gamma.hat.mat[,ind]!=0)))
  	}
  	set.list = c(set.list,list(set = snp.set))

        #--
  	Mat.V = NULL
 	if (length(snp.set)==0){
 	r1=length(protein.set)
 	s1=length(protein.set)
	}else{
  	r1=length(protein.set)+1
  	s1=length(protein.set)+1
	}

  	for(r in 1:r1){
  	  Mat.V1 = NULL
  	  for(s in 1:s1)
  	  {
  	    #-
  	    Mat.V.rs = matrix(0,length(set.list[r]$set),length(set.list[s]$set))
  	    for(j in 1:length(protein.set)){
  	      Mat.V.rs =
  	        Mat.V.rs +
  	        (-2*n2*(n2+n3)/(n1+n2)*BetaHat[j]*sum(Est.Cov[j,-j]*c(BetaHat[-j],1))+
  	           (n3^2-n2^2)/(n1+n2)*BetaHat[j]^2*Est.Cov[j,j])*
  	        Est.Cov.Z[set.list[r]$set,set.list[j]$set]%*%
  	        solve(Est.Cov.Z[set.list[j]$set,set.list[j]$set])%*%
  	        Est.Cov.Z[set.list[j]$set,set.list[s]$set]
  	    }
  	    for(j in 1:length(protein.set)){
  	      for(k in (1:length(protein.set))[-j]){
  	        Mat.V.rs =
  	          Mat.V.rs +
  	          BetaHat[j]*BetaHat[k]*Est.Cov[j,k]*(n2+n3)^2/(n1+n2)*
  	          Est.Cov.Z[set.list[r]$set,set.list[j]$set]%*%
  	          solve(Est.Cov.Z[set.list[j]$set,set.list[j]$set])%*%
  	          Est.Cov.Z[set.list[j]$set,set.list[k]$set]%*%
  	          solve(Est.Cov.Z[set.list[k]$set,set.list[k]$set])%*%
  	          Est.Cov.Z[set.list[k]$set,set.list[s]$set]
  	      }
  	    }
  	    Mat.V.rs =
  	      Mat.V.rs + Est.Var.total*Est.Cov.Z[set.list[r]$set,set.list[s]$set]*(n2+n3)
  	    #-
  	    Mat.V1 = cbind(Mat.V1,Mat.V.rs)
  	  }
  	  Mat.V = rbind(Mat.V,Mat.V1)
  	}
  	#--
  	Mat.W = Mat.V/(n2+n3)

  	### calculate matrix Gamma
  	Mat.Gamma = matrix(0,0,0)
  	for(ind in protein.set)
  	{
  	  nonzero.ind = which(gamma.hat.mat[,ind]!=0)
  	  Mat.Gamma =
  	    Matrix::bdiag(Mat.Gamma,
  	                  gamma.hat.mat[nonzero.ind,ind])
  	}
  	Mat.Gamma = Matrix::bdiag(Mat.Gamma,diag(length(snp.set)))
  	Mat.Gamma = as.matrix(Mat.Gamma)

  	### calculate matrix S
  	Mat.S = NULL


  	for(ind1 in 1:r1){
  	  Mat.S1 = NULL
  	  for(ind2 in 1:s1){
  	    Mat.S1 =
  	      cbind(Mat.S1,Est.Cov.Z[set.list[ind1]$set,set.list[ind2]$set,drop=F])
  	  }
  	  Mat.S = rbind(Mat.S,Mat.S1)
  	}

 	 ### calculate covariance matrix
  	Final.Est.Cov =
    	solve(t(Mat.Gamma)%*%Mat.S%*%Mat.Gamma)%*%
    	t(Mat.Gamma)%*%Mat.W%*%Mat.Gamma%*%
   	 solve(t(Mat.Gamma)%*%Mat.S%*%Mat.Gamma)

   	delta=as.numeric(1- 2*sum(cov_DesignX_sim.Y.BIC* BetaAlpha.hat) +
                       t(BetaAlpha.hat)%*% Cov.Design.BIC%*% BetaAlpha.hat)

    	A=list(Est = BetaAlpha.hat,
              Sd = sqrt(diag(Final.Est.Cov)/(n2+n3)),
              Sd.naive1 =sqrt(delta*diag(solve(Cov.Design.BIC))/(n2+n3)),
              protein.set = protein.set,snp.set = snp.set
              )

	return(A)

}




