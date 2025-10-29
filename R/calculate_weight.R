#' Calculate weights by by performing principal component analysis (PCA) on genetic correlation matrix among traits
#'
#' Calculate weights by by performing principal component analysis (PCA) on genetic correlation matrix among traits
#' @param K number of traits
#' @param corr genetic correlation matrix
#' @param varcut variance cutoff
#' @details Calculate weights by by performing principal component analysis (PCA) on genetic correlation matrix among traits
#' @return a vector of weights
#' @references Zhai, S., Guo, B., and Shen, J., 2025. Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning.
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "CS", sparseness = "no",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.5, rho_C = 0.2,
#' K=4, m=2000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' re <- calculate_weight(K = 4, corr = dat$corr, varcut = 0.8)
#' hist(re$mtPRS)
#' }
#'
calculate_weight <- function(K, corr, varcut){
  D <- nearPD(corr[1:K,1:K], corr = TRUE); mat1 <- D$mat %>% as.matrix()
  x = princomp(covmat=mat1, cor=TRUE)
  pc <- x$sdev^2/sum(x$sdev^2)
  pcc <- 0; ct <- 0
  while(pcc < varcut){
    ct <- ct + 1
    pcc <- pcc + pc[ct]
  }
  e1 <- eigen(mat1)
  if(ct > 1){w <- apply(e1$vectors[,1:ct], 1, sum)}
  if(ct == 1){w <- e1$vectors[,1]}
  w[which(is.na(w))] <- 0
  w
}
