#' Calculate individual PRS for each trait based on refined effect size estimates from PolyPred
#'
#' Calculate individual PRS for each trait based on refined effect size estimates from PolyPred
#' @param K number of traits
#' @param base list of disease GWAS summary statistics of K traits from base cohort
#' @param pcut p-value cutoff
#' @param G.target genotype in target cohort
#' @details calculate_PRS needs disease GWAS summary statistics from K traits, assuming the effect size estimates have been refined by PolyPred after fine-mapping
#' @return a matrix of K individual PRSs
#' @references Zhai, S., Guo, B., and Shen, J., 2025. Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning.
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "CS", sparseness = "no",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.5, rho_C = 0.2,
#' K=4, m=2000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' re <- calculate_PRS(K = 4, base = dat$base, pcut = 1, G.target = dat$target[,-c(1,2)])
#' head(re)
#' }
#'
calculate_PRS <- function(K, base, pcut, G.target){
  X <- double()
  for(k in 1:K){
    ss <- base[[k]]
    index <- which(ss$pval > pcut)

    beta <- ss$beta; beta[index] <- 0

    prs <- G.target %*% beta %>% as.vector()

    X <- cbind(X, prs)
  }
  for (k in 1:K) {
    if(sum(X[,k] != 0) > 0){X[,k] <- scale(X[,k])}
  }
  colnames(X) <- names(base)
  return(X)
}
