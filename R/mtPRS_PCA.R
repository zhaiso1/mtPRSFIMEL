#' Construct mtPRS by performing PCA on genetic correlation matrix among traits
#'
#' Construct mtPRS by performing PCA on genetic correlation matrix among traits
#' @param base list of disease GWAS summary statistics of K traits
#' @param target data frame in target cohort, including Y, T (if the target cohort is from a PGx study), and G
#' @param validation data frame from an independent validation cohort, including Y, T (if the validation cohort is from a PGx study), and G
#' @param corr genetic correlation matrix among K traits
#' @param pcut p-value cutoff for C+T method to construct individual stPRS
#' @param varcut variance cutoff to choose top principal components (PCs)
#' @param K number of traits
#' @param phenotype indicator of phenotype in target cohort, either "dis" or "pgx"
#' @details Construct mtPRS by performing PCA on genetic correlation matrix among traits
#' @return a list of weights, individual stPRSs, and mtPRS-PCA
#' @references Zhai, S., Guo, B., and Shen, J., 2025. Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning.
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "CS", sparseness = "no",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.5, rho_C = 0.2,
#' K=4, m=2000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' re <- mtPRS_PCA(dat$base, dat$target, dat$validation, dat$corr, pcut = 1, varcut = 0.8, K = 4, phenotype = "pgx")
#' hist(re$mtPRS_validation)
#' }
#'
mtPRS_PCA <- function(base, target, validation, corr, pcut=0.05, varcut=0.8, K=4, phenotype="pgx"){
  if(phenotype == "pgx"){G.target <- target[,-c(1:2)] %>% as.matrix(); Y.target <- target[,1] %>% as.vector(); T.target <- target[,2] %>% as.vector()}
  if(phenotype == "dis"){G.target <- target[,-1] %>% as.matrix(); Y.target <- target[,1] %>% as.vector()}

  if(phenotype == "pgx"){G.validation <- validation[,-c(1:2)] %>% as.matrix(); Y.validation <- validation[,1] %>% as.vector(); T.validation <- validation[,2] %>% as.vector()}
  if(phenotype == "dis"){G.validation <- validation[,-1] %>% as.matrix(); Y.validation <- validation[,1] %>% as.vector()}

  X.target <- calculate_PRS(K, base, pcut, G.target)
  X.validation <- calculate_PRS(K, base, pcut, G.validation)

  X.target <- X.target %>% scale(); X.validation <- X.validation %>% scale()

  w <- calculate_weight(K, corr, varcut)

  PRS.target <- X.target %*% w %>% as.vector()
  PRS.validation <- X.validation %*% w %>% as.vector()

  X.target <- X.target %>% as.data.frame(); X.validation <- X.validation %>% as.data.frame()

  re <- list(w=w, stPRS_target=X.target, stPRS_validation=X.validation, mtPRS_target=PRS.target, mtPRS_validation=PRS.validation)
  re
}
