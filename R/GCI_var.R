#' CI calculation for the total variance parameter in a single-factor random-effects ANOVA model
#'
#' Consider a random-effects ANOVA model, \eqn{Y_{ij}=\mu+u_i+\epsilon_{ij}}, where the random-effects \eqn{u_i\sim N(0,\sigma_a^2)} and
#'  the random error components \eqn{\epsilon_{ij}\sim N(0,\sigma^2)}. We want to compute CI for the total variance parameter \eqn{\sigma_a^2+\sigma^2}.
#' An asymptotic approach and a generalized pivotal test approach are implemented. See the two references. 
#' 
#' @param Y observed outcomes in an ANOVA model
#' @param A factor level
#' @param alpha desired significance level for the confidence interval
#' @param Nmax total number of generalized pivotal statistics to be generated
#' @return
#' \describe{
#'   \item{GCI}{ generalized CI from Park and Burdick (2004) }
#'   \item{BCI}{ asymptotic CI from Burdick and Graybill (1984) }
#' }
#'
#' @export
#' @references
#' Burdick, R.K., Graybill, F.A., 1984. Confidence intervals on linear combinations of variance components in the unbalanced one-way classification. Technometrics 26, 131–136.
#'
#' Park, D.J., Burdick, R.K., 2004. Confidence intervals on total variance in a regression model with an unbalanced one-fold nested error structure. Comm. Statist. Theory Methods 33, 2735–2743.
BGPsv <- function(Y,A, alpha=0.1, Nmax=1e5){
  obj = QUADinf0(Y,A)
  r = obj$r;  s = obj$s; 
  SSE = obj$SSE; MSM = obj$MSM
  Ql = obj$Ql; dl = obj$dl; h = s/sum(1/dl)
  ## BG
  MSE = SSE/r
  gh = MSM + (1-1/h)*MSE
  SU2 = MSM*h
  alphaL = alpha/2; alphaU=alpha/2
  G1 = 1-s/qchisq(alphaL,s, lower=FALSE)
  G2 = 1-r/qchisq(alphaL,r, lower=FALSE)
  H1 = s/qchisq(alphaU,s)-1
  H2 = r/qchisq(alphaU,r)-1
  MSE2 = MSE^2
  bL = gh - sqrt(G1^2*SU2^2+G2^2*(h-1)^2*MSE2)/h
  bU = gh + sqrt(H1^2*SU2^2+H2^2*(h-1)^2*MSE2)/h
  bci = c(bL,bU)
  ## GPQ/GCI
  RE = rchisq(Nmax, r)
  U = rchisq(Nmax, s)
  tx = rep(0, Nmax)  ## recommended by PB 2004. Conservative
  id = which( sum(Ql)/(SSE/RE)>U)
  if(length(id)>0){
    for(i in id){
      tx[i] = T2sol(Ql,SSE,dl, RE[i],U[i])
    }
  }
  gci = quantile(tx, c(alpha/2,1-alpha/2))
  return( list(GCI=gci, BCI=bci) )
}
## Essentially: SSE/RE is the simulated \sigma^2, hence
##  by definition, we need to search for S2S \in (\sigma^2, \infty)
T2sol <- function(Ql,SSE,dl, RE,U){
  fc = function(s2s)  sum( Ql/(SSE/RE+dl*(s2s-SSE/RE)) ) - U
  x2 = SSE/RE+1
  while(fc(x2)>0){
    x2 = x2*2
  }
  uniroot(fc, c(SSE/RE, x2))$root
}
QUADinf0 <- function(Y,A){
  B = model.matrix(~factor(A)-1)
  g = dim(B)[2]; s = g-1
  n = length(Y)
  p = 1; r = n - p - s
  E = t( t(B) - colMeans(B) )
  seb = svd(E%*%t(B), nv=0, nu=s); Ub = seb$u
  Ue = cbind(1/sqrt(n),Ub)
  H2E = Ue%*%(t(Ue)%*%E)
  s2e = svd(H2E, nv=0, nu=s); U2 = s2e$u
  ##
  dl = s2e$d[1:s]^2
  h = s/sum(1/dl)
  H2Y = colSums(t(Ue)*(colSums(Ue*Y)))
  Y2 = Y-H2Y
  SSE = sum(Y2^2)
  H1Y = mean(Y)
  Y1 = Y-H1Y
  FY = H2Y-H1Y
  Ql = colSums(FY*U2)^2
  MSM = sum(Ql/dl)/s
  return( list(SSE=SSE,MSM=MSM, dl=dl,Ql=Ql, r=r,s=s) )
}

BGPxv <- function(Y,A, X, alpha=0.1, Nmax=1e5){
  obj = QUADinf(Y,X,A)
  r = obj$r;  s = obj$s; 
  SSE = obj$SSE; MSM = obj$MSM
  Ql = obj$Ql; dl = obj$dl; h = s/sum(1/dl)
  ## BG
  MSE = SSE/r
  gh = MSM + (1-1/h)*MSE
  SU2 = MSM*h
  alphaL = alpha/2; alphaU=alpha/2
  G1 = 1-s/qchisq(alphaL,s, lower=FALSE)
  G2 = 1-r/qchisq(alphaL,r, lower=FALSE)
  H1 = s/qchisq(alphaU,s)-1
  H2 = r/qchisq(alphaU,r)-1
  MSE2 = MSE^2
  bL = gh - sqrt(G1^2*SU2^2+G2^2*(h-1)^2*MSE2)/h
  bU = gh + sqrt(H1^2*SU2^2+H2^2*(h-1)^2*MSE2)/h
  bci = c(bL,bU)
  ## GPQ/GCI
  ## tmax = 1e5
  RE = rchisq(Nmax, r)
  U = rchisq(Nmax, s)
  ## tx = SSE/RE  ## better. same coverage with shorter CI.
  tx = rep(0, Nmax)  ## recommended by PB 2004. too conservative
  id = which( sum(Ql)/(SSE/RE)>U)
  if(length(id)>0){
    for(i in id){
      tx[i] = T2sol(Ql,SSE,dl, RE[i],U[i])
      ## cat(i, tx[i], '\n')
    }
  }
  gci = quantile(tx, c(alpha/2,1-alpha/2))
  return( list(GCI=gci, BCI=bci) )
}


#' Generalized confidence interval for the total variance in a single-factor random-effects ANOVA model
#'
#' Consider a random-effects ANOVA model, \eqn{Y_{ij}=\mu+u_i+\epsilon_{ij}}, where the random-effects \eqn{u_i\sim N(0,\sigma_a^2)} and
#'  the random error components \eqn{\epsilon_{ij}\sim N(0,\sigma^2)}. We want to compute CI for the total variance parameter \eqn{\sigma_a^2+\sigma^2}.
#' We adopt the generalized pivotal test approach to computing the CI, based on
#' sequentially ``calculating" the within-group variance parameter \eqn{\sigma^2} and the between-group variance parameter \eqn{\sigma_a^2}.
#'
#' @param Y observed outcomes in an ANOVA model
#' @param A factor level
#' @param alpha desired significance level for the confidence interval
#' @param Nmax total number of generalized pivotal statistics to be generated
#' @return
#' \describe{
#'   \item{GCI}{ generalized CI }
#'   \item{Qt}{ generated generalized pivotal statistics for the total variance }
#' }
#'
#' @export
#' @references
#' X,X. and Wu,B. (2019) Generalized confidence interval calculation for the total variance in a single-factor random-effects ANOVA model
#'             with application to medical device comparison problems. tech report.
#' @examples
#' A = rep(1:10, times=5:14)
#' s2a = 0.05; s2 = 0.95
#' Y = rnorm(10)[A]*sqrt(s2a) + rnorm(length(A))*sqrt(s2)
#' BGPsv(Y,A)
#' aa = GCIsv(Y,A); aa$GCI
GCIsv <- function(Y,A, alpha=0.1, Nmax=1e5){
  ## sufficient stats
  ng = as.vector(table(A)); m = length(ng); N = sum(ng)
  sse = sum(tapply(Y, A, var)*(ng-1), na.rm=TRUE)
  mus = tapply(Y, A, mean)
  ## inf
  Qe = rchisq(Nmax,N-m)
  Qr = rchisq(Nmax,m-1)
  s2a = rep(0, Nmax)
  for(i in 1:Nmax){
    s2a[i] = S2Asol(sse/Qe[i],ng,mus,Qr[i])
  }
  ## T2
  Qt = sse/Qe + s2a
  ci = quantile(Qt, c(alpha/2,1-alpha/2))
  return( list(GCI=ci,Qt=Qt) )
}
S2Asol <- function(s2,ng,mus,Qr){
  ## s2 = sse/Qe
  fc = function(s2a){
    W = 1/(s2a+s2/ng)
    yw = sum(mus*W)/sum(W)
    sum(W*(mus-yw)^2) - Qr
  }
  v0 = fc(0)
  if(v0<=0) return(0)
  x2 = 1
  while(fc(x2)>0){
    x2 = x2*2
  }
  ans = uniroot(fc, c(0,x2))$root
  return(ans)
}

GCIxv <- function(Y,A, X, alpha=0.1, Nmax=1e5){
  obj = QUADinf(Y,X,A)
  r = obj$r;  s = obj$s; 
  SSE = obj$SSE; MSM = obj$MSM
  Ql = obj$Ql; dl = obj$dl; h = s/sum(1/dl)
  ## BG
  MSE = SSE/r
  gh = MSM + (1-1/h)*MSE
  SU2 = MSM*h
  alphaL = alpha/2; alphaU=alpha/2
  G1 = 1-s/qchisq(alphaL,s, lower=FALSE)
  G2 = 1-r/qchisq(alphaL,r, lower=FALSE)
  H1 = s/qchisq(alphaU,s)-1
  H2 = r/qchisq(alphaU,r)-1
  MSE2 = MSE^2
  bL = gh - sqrt(G1^2*SU2^2+G2^2*(h-1)^2*MSE2)/h
  bU = gh + sqrt(H1^2*SU2^2+H2^2*(h-1)^2*MSE2)/h
  bci = c(bL,bU)
  ## GPQ/GCI
  ## tmax = 1e5
  RE = rchisq(Nmax, r)
  U = rchisq(Nmax, s)
  tx = rep(0, Nmax)
  id = which( sum(Ql)/(SSE/RE)>U)
  if(length(id)>0){
    for(i in id){
      tx[i] = S2A.sol(Ql,SSE,dl, RE[i],U[i])
    }
  }
  tx = tx + SSE/RE
  gci = quantile(tx, c(alpha/2,1-alpha/2))
  return( list(GCI=gci, BCI=bci) )
}
##
S2A.sol <- function(Ql,SSE,dl, RE,U){
  fc = function(s2a)  sum( Ql/(SSE/RE+dl*s2a) ) - U
  x2 = 1
  while(fc(x2)>0){
    x2 = x2*2
  }
  uniroot(fc, c(0, x2))$root
}




## ' Construct quadratic forms that can be used to estimate variance parameters in a random-effects model
## '
## ' For a linear mixed-effects model with one-fold nested error structure, we can construct quadratic
## ' forms to estimate variance parameters and build confidence intervals in the downstream analysis.
## '
## ' Some technical details:
## ' Under reduced model, we construct a projection matrix based on X with \eqn{H_1=X(X'X)^{-1}X }.
## ' Under the augmented model, we construct another porjection matrix based on \eqn{X_2=(X,BB')} with
## ' \eqn{H_2=X_2(X_2'X_2)^{+}X_2'}. Here the superscript '+' denotes the Moore-Penrose inverse, and 
## ' B is contructed as the design matrix of group factors A, s.t., \eqn{Var(Y)=\sigma_a^2BB'+\sigma^2I_n}.
## ' In a matrix notation, we have \eqn{Y = X\beta + BU+\epsilon}, where \eqn{U\sim N(0,\sigma_a^2I_r)}.
## '
## ' Let \eqn{F=H_2-H_1, W=FBB'F}. Define \eqn{SSE=Y'(I_n-H_2)Y} and \eqn{SSR = Y'F'W^+FY}.
## ' Denote eigen values of W as \eqn{(\lambda_1,\cdots,\lambda_s)}. Here \eqn{s = rank(X_2)-rank(X), r=n-rank(X_2)}.
## ' We can show that \eqn{SSE\sim \sigma^2\chi_r^2/r}, and
## ' \eqn{SSR\sim (\sigma_a^2+\sigma^2/h)\sum_{j=1}^s (1-\rho+\lambda_j)/(\lambda_j(1-\rho)/h+\lambda_j\rho)\chi_1^2 }.
## ' Here \eqn{\rho=\sigma_a^2/(\sigma_a^2+\sigma^2)} and \eqn{ h = s/(\sum_j1/\lambda_j) }.
## '
## ' Computationally we can proceed as follows.
## ' Denote regression decomposition \eqn{ B = X\alpha_b + E }, where \eqn{E=(I_n-H_1)X, E'X=0}. We show
## ' (1) \eqn{ B'F = E'F = E'H_2, W=H_2EE'H_2 }
## ' (2) the eigenvalues of E'E are the same as W.
## ' (3) construction of H_2 can be equivalently based on (X,EB'). Specifically denote the SVD, \eqn{EB'=U_bD_bV_b'},
## '     \eqn{X=UDV'}, and \eqn{H_2E = U_eD_eV_e'}.
## ' We have
## ' (1) \eqn{H_1 = UU', H_2 = (U,U_b)(U,U_b)'}.
## ' (2) \eqn{ H_2Y = (U,U_b)(U,U_b)'Y }
## ' (3) \eqn{ W = H_2EE'H_2 = U_eD_e^2U_e'}
## '
## ' Alternatively, we can do the following.
## ' (1) First we estimate \eqn{\sigma^2} using projection matrix based on \eqn{X_2=(X,BB')} with \eqn{H_2=X_2(X_2'X_2)^{+}X_2'}.
## ' (2) Consider \eqn{(H_2-H_1)Y\sim N(0, (B-H_1B)'(B-H_1B)'\sigma_a^2+(H_2-H_1)\sigma^2) },
## '     we can readily construct a \eqn{\chi_s^2} statistic to estimate \eqn{\sigma_a^2}.
## ' (3) Consider \eqn{H_1Y\sim N(X\beta,H_1BB'H_1\sigma_a^2+H_1\sigma^2)},
## '     we can construct a p-dim multivariate normal vector to estimate \eqn{\beta}.
## '
## ' 
## ' @param Y outcome vector
## ' @param X fixed-effects covariates
## ' @param A group indicators
## '
## ' @return
## ' \describe{
## '   \item{SSE,SSR}{the error and regression sum of squares}
## '   \item{r, s}{the error and regression DFs}
## ' }
## ' @references
## ' Eubank, L., Seely, J., Lee, Y. (2001). Unweighted mean squares for the general two variance component mixed model. In Proceedings of Graybill Conference, Ft. Collins, CO., June, pp. 281-290.
## '
## ' El-Bassiouni, M. Y. (1994). Short confidence intervals for variance components. Communications in Statistics—Theory and Methods 23(7):1915–1933.
## '
## ' Olsen, A., Seely, J., Birkes, D. (1976). Invariant quadratic unbiased estimation for two variance components. Annals of Statistics 4:878–890.
QUADinf <- function(Y,X,A){
  B = model.matrix(~factor(A)-1)
  g = dim(B)[2]; s = g-1
  n = length(Y)
  X1 = cbind(rep(1,n), X); p = dim(X1)[2]; r = n - p - s
  sx = svd(X1,nu=p,nv=0); U = sx$u
  E = B - U%*%(t(U)%*%B)
  seb = svd(E%*%t(B), nv=0, nu=s); Ub = seb$u
  Ue = cbind(U,Ub)
  H2E = Ue%*%(t(Ue)%*%E)
  s2e = svd(H2E, nv=0, nu=s); U2 = s2e$u
  ##
  dl = s2e$d[1:s]^2
  h = s/sum(1/dl)
  H2Y = colSums(t(Ue)*(colSums(Ue*Y)))
  Y2 = Y-H2Y
  SSE = sum(Y2^2)
  H1Y = colSums(t(U)*(colSums(U*Y)))
  Y1 = Y-H1Y
  FY = H2Y-H1Y
  Ql = colSums(FY*U2)^2
  MSM = sum(Ql/dl)/s
  return( list(SSE=SSE,MSM=MSM, dl=dl,Ql=Ql, r=r,s=s) )
}



