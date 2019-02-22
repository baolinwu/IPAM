# IPAM
 - An R package for statistical tests in various ANOVA models
 - Reference
    - X,X., X,X., and Wu,B. (2019) A generalized Tukey method for multiple pairwise comaprisons in ANOVA model. *tech rep*.
    -  X,X. and Wu,B. (2019) Generalized confidence interval calculation for the total variance in a single-factor random-effects ANOVA model  with application to medical device comparison problems. *tech rep*.
 - Sample R codes
```R
 ## install the package
 devtools::install_github('baolinwu/IPAM')
 library(IPAM)
 ## multiple testing
 GTukey(c(15,10,5), W=3:1)
 ## CI for total var
 A = rep(1:10, times=5:14)
 s2a = 0.05; s2 = 0.95
 Y = rnorm(10)[A]*sqrt(s2a) + rnorm(length(A))*sqrt(s2)
 BGPsv(Y,A)
 aa = GCIsv(Y,A); aa$GCI
```
