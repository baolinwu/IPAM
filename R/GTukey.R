#' A generalized Tukey method for multiple pairwise comparisons in ANOVA model
#'
#' Given an ANOVA model with r levels, we are interested in testing the pairwise mean differences (1 vs 2,...; 2 vs 3,... etc). 
#' We consider a weighted multiple hypothesis testing problem, where we allow different tests
#' to get unequal type I errors, and in the meantime, we want to control the overall type I errors at a pre-specified level.
#'
#' @param nr the sample sizes for each factor level.
#' @param alpha overall significance level
#' @param W the relative type I errors assigned to each comparison.
#' @return
#' \describe{
#'   \item{Bonf}{ critial values computed from Bonferroni correction }
#'   \item{Tukey}{ Tukey HSD approximate value (common to all) }
#'   \item{WT}{ weighted Tukey HSD critical values }
#'   \item{FP}{ overall familywise type I errors }
#'   \item{FPi}{ individual test type I errors }
#' }
#' @export
#' @references
#' X,X., and Wu,B. (2019) A generalized Tukey method for multiple pairwise comparisons in ANOVA model.
#' @examples
#' GTukey(c(15,10,5), W=3:1)
GTukey <- function(nr, alpha=0.05, W=NULL){
  r = length(nr); N = sum(nr); K = r*(r-1)/2
  if(is.null(W)){
    W = rep(1/K, K)
  } else{
    W = W/sum(W)
  }
  id = which(W>0)
  W1 = W[id]
  ## Joint prob: best/exact/smallest critical value
  require(mvtnorm)
  V = diag(1/nr); pcm = matrix(0, K,r)
  it = 0
  for(i in 1:(r-1)){
    ij = it+1:(r-i)
    pcm[ij,i] = 1; pcm[cbind(ij,(i+1):r)] = -1
    it = it+r-i
  }
  ##
  CV = pcm%*%(t(pcm)/nr); a1 = sqrt(diag(CV))
  sig = t(CV/a1)/a1; sig1 = sig[id,id]
  ##
  func = function(ai){
    x3 = qt(ai*W1/2,N-r,lower=FALSE)
    1-pmvt(-x3, x3, df=N-r,sigma=sig1,algorithm=GenzBretz(maxpts=round(1e4/alpha),abseps=alpha*1e-5))-alpha
  }
  mc3 = uniroot(func, c(alpha,alpha*K))$root
  x3 = qt(mc3*W1/2,N-r,lower=FALSE)
  ## Bonferroni
  Bf.alpha = alpha*W1
  brK = qt(Bf.alpha/2,N-r,lower=FALSE)
  Bf.fp = c( 1-pmvt(-brK, brK, df=N-r,sigma=sig1,algorithm=GenzBretz(maxpts=round(1e4/alpha),abseps=alpha*1e-5)) )
  ## Tukey HSD (approx) for all tests
  tk = qtukey(alpha,r,N-r, lower=FALSE)/sqrt(2)
  tk.alpha = 2*pt(-tk,N-r)
  tk.fp = c( 1-pmvt(-rep(tk,length(W1)), rep(tk,length(W1)), df=N-r,sigma=sig1,algorithm=GenzBretz(maxpts=round(1e4/alpha),abseps=alpha*1e-5)) )
  FP = list(Bonf=Bf.fp,Tukey=tk.fp, WT=alpha)
  FPi = list(Bonf=Bf.alpha, Tukey=tk.alpha, WT=mc3*W1)
  ##
  ## return( list(Bonf=brK, Tukey=tk, WT=x3, WTalpha=mc3, FP=FP,FPi=FPi) )
  return( list(Bonf=brK, Tukey=tk, WT=x3, FP=FP,FPi=FPi) )
}

