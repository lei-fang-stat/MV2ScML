# MV2ScML

MV2ScML is designed to perform the *Multivariate Two-stage Constrained Maximum Likelihood (MV-2ScML)* method to conduct the causal inference for multiple exposures. Here we provide a tutorial on how to run this R package.

## Installation
You can install it using the following command:
```r
# install.packages("devtools")
devtools::install_github("lei-fang-stat/MV2ScML")
```
## Load R packages
We also need to install the following R package:

- **MASS**: This package is available on CRAN. We use `MASS::mvrnorm()`
  function to generate random variables.
- **kernlab**: This package is to get kernel matrix calculation  
  
- **lasso2**: This package is not available on CRAN for R 4.2.2, it is available on GitHub at https://github.com/cran/lasso2.

## Prepare the summary input 
We first need the harmonized GWAS summary statistics for multiple exposures and the outcome. We also need the LD matrix for all the SNPs (IVs) selected among all exposures.

The GWAS data for each exposure and 
```
exposure.GWAS.list[[1]]
##GWAS summary statistics for exposure 1 with selected SNPs as IVs
      SNP   Beta          SE           Z_STAT         pval
   snp_1    0.05156978  0.01412614   3.650662 2.642384e-04
   snp_2    0.09877343  0.01407580   7.017254 2.564895e-12
   snp_3    0.15395708  0.01397632  11.015565 6.716782e-28
   snp_4    0.17606737  0.01392399  12.644890 4.242292e-36
   snp_5    0.07100900  0.01410926   5.032795 5.003345e-07
   snp_6    0.22396436  0.01378565  16.246201 7.065795e-58
   snp_7    0.20514517  0.01384412  14.818213 1.188866e-48
   snp_8    0.15570615  0.01397244  11.143802 1.659854e-28
   snp_9    0.14119920  0.01400325  10.083317 1.101571e-23
   snp_10   0.02936373  0.01413887   2.076810 3.787010e-02
```

```
outcome.GWAS.comb
##GWAS summary statistics for the outcome
   SNP     Beta          SE         Z_STAT        pval
   snp_1  0.1197827052 0.006279271 19.07589497  1.497420e-80
   snp_2  0.1361674495 0.006265898 21.73151385  9.439951e-104
   snp_3  0.1248120054 0.006275351 19.88924721  2.391318e-87
   snp_4  0.0647646659 0.006311530 10.26132618  1.178026e-24
        ⋮
   snp_45 -0.0026410088 0.006324786 -0.41756491  6.762688e-01
   snp_46 -0.0040127625 0.006324757 -0.63445319  5.257909e-01
```
```
exposure.GWAS.NULL.list[[1]]
##GWAS summary statistics for exposure 1 with shared NULL SNPs (no association with all exposures)
        SNP       Beta          SE           Z_STAT         pval
     snp_NULL1 -6.370629e-03  0.01414468 -4.503905e-01 0.6524484093
     snp_NULL2  4.952169e-03  0.01414479  3.501055e-01 0.7262742980
     snp_NULL3  1.458810e-03  0.01414495  1.031329e-01 0.9178616062
     snp_NULL4 -1.475113e-04  0.01414496 -1.042854e-02 0.9916797957
                                  ⋮
    snp_NULL499 -5.802417e-03  0.01414473 -4.102177e-01 0.6816638605
    snp_NULL500  2.874170e-03  0.01414491  2.031947e-01 0.8389911042
```

```
#Est.Cov.D estimated covariance for each residuals in stage 1
```
parameter input for run the two sample version of MV2ScML
```
#n1 sample: size of stage 1
#n3 sample: size of stage 2
#num.protein: number of exposures
#num.snp: number of unique SNPs(IVs) for all exposures
#exposure.GWAS.list: list of GWAS result for exposures with preselected IVs for each exposure.
#exposure.GWAS.NULL.list: list of GWAS result for exposures with NULL SNPs for each exposure.
#cor.Z: correlation matrix of in-sample LD for stage 2 or from reference panel with colnames of SNPs id.
#SNP: list of IV names for each exposure
#SNP.comb: unique of the combined IV names
#outcome.GWAS.i is the stage2 GWAS statistics for each exposure
#outcome.GWAS.comb is the stage2 GWAS statistics for all IVs 
#K.vec the tunning parameter for BIC selection, default is 0:(num.snp-2)
#res.D residual for each exposure

#snp.set.provide the invalid IVs, default is NULL unless user specified

MV.2ScML.nooverlap(n1,n3,num.snp,exposure.GWAS.list,
                            cor.Z.sample3,
                            outcome.GWAS.comb,
                            K.vec,Est.Cov.D,snp.set.provide)
```

The Stage 1 fitted models for UKB-PPP are available on https://www.synapse.org/Synapse:syn64041817 
