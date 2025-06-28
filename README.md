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

- **kernlab**: This package is to get kernel matrix type of calculation  
  
- **lasso2**: This package is not available on CRAN for R 4.2.2, it is available on GitHub at https://github.com/cran/lasso2.

## Prepare the input 
We first need the harmonized GWAS summary statistics for multiple exposures and the outcome. We also need the LD matrix for all the SNPs selected among all exposures and selected NULL SNPs for each exposure.

We first need the list of GWAS summary statistics for each exposure for some preselected SNPs as IVs, as shown below. In terms of format, we require the input colnames to match "SNP", "Beta", and "SE". 
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
We then need the GWAS summary statistics for the outcome, combining all the selected SNPs from each exposure.
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
We also need the list of GWAS summary statistics for each exposure for NULL SNPs(selected based on large p-values for all exposures). We need this to compute the variance component between each exposure.
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

We need the LD matrix with rownames and colnames for all the SNPs, including NULL SNPs.
```
              snp_1         snp_2        snp_3        snp_4        snp_5
snp_1   1.000000000  0.1519173896  0.014933198 -0.002617276 -0.021856415
snp_2   0.151917390  1.0000000000  0.140473615  0.015120372 -0.020503246
snp_3   0.014933198  0.1404736149  1.000000000  0.187000017 -0.041684178
snp_4  -0.002617276  0.0151203717  0.187000017  1.000000000 -0.008562873
snp_5  -0.021856415 -0.0205032463 -0.041684178 -0.008562873  1.000000000
snp_6   0.001080836 -0.0157352531 -0.018348674 -0.001580189  0.153201075
snp_7  -0.004526601 -0.0003082802  0.012465715 -0.015170380  0.061496605
snp_8   0.010369612  0.0061495440 -0.005988474 -0.003906237 -0.020700674
snp_9  -0.012574050 -0.0121604418  0.011151235 -0.004270432 -0.018394141
snp_10  0.011202148  0.0135503486 -0.003665583 -0.024114023  0.019604671
```

## Run MV2ScML
```
#n1 sample: size of stage 1
#n3 sample: size of stage 2
#exposure.GWAS.list: list of GWAS result for exposures with preselected IVs for each exposure.
#exposure.GWAS.NULL.list: list of GWAS result for exposures with NULL SNPs for each exposure.
#cor.Z: correlation matrix of in-sample LD for stage 2 or from reference panel with colnames of SNPs id.
#outcome.GWAS.comb is the stage2 GWAS statistics for all IVs 
#K.vec range of number of invalid IVs for BIC selection, default is 0:(num.snp-2)
#snp.set.provide the invalid IVs, default is NULL unless user specified

MV.2ScML(n1,n3,exposure.GWAS.list,
               exposure.GWAS.NULL.list,
               outcome.GWAS.comb,
               cor.Z,
               K.vec)
```
## Output 
The output gives the estimated causal effects, standard errors and p-values for all exposures, along with the identified invalid IVs. 
```
$Estimate
                   [,1]
exposure1   0.030156256
exposure2   0.255923311
exposure3   0.341756025
exposure4   0.231838494
exposure5   0.215322935
exposure6  -0.050364047
exposure7  -0.035213049
exposure8   0.006757947
exposure9  -0.014931405
exposure10  0.008498412
snp_1       0.066412569
snp_7       0.058610966
snp_11      0.066428206
snp_15      0.052966577
snp_18     -0.025768041
snp_19      0.037000769
snp_23      0.052286420
snp_27      0.060834744
snp_31      0.056964122
snp_35      0.068057394
snp_39      0.054895088
snp_40      0.030009996
snp_43      0.056762508

$Standard.error
 [1] 0.032834655 0.030873660 0.049409701 0.026713126 0.025515408 0.029996989
 [7] 0.030191396 0.021486927 0.026838121 0.034150693 0.012773010 0.008163508
[13] 0.006946620 0.007624294 0.009777631 0.007692203 0.007233707 0.008172872
[19] 0.006544262 0.006122865 0.006262617 0.006486777 0.006109152

$pval.exposure
 [1] 3.583950e-01 2.220446e-16 4.620304e-12 0.000000e+00 0.000000e+00
 [6] 9.315788e-02 2.434821e-01 7.531304e-01 5.779712e-01 8.034765e-01
```
## Reference
Lei Fang, Haoran Xue, Zhaotong Lin, Wei Pan,
Multivariate proteome-wide association study to identify causal proteins for Alzheimer disease,
The American Journal of Human Genetics,
https://doi.org/10.1016/j.ajhg.2024.12.010.
The Stage 1 fitted models for UKB-PPP in the paper are available on https://www.synapse.org/Synapse:syn64041817.
