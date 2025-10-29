#' Construct mtPRS by performing ensemble learning
#'
#' Construct mtPRS by performing ensemble learning
#' @param base list of disease GWAS summary statistics of K traits
#' @param target data frame in target cohort, including Y, T (if the target cohort is from a PGx study), and G
#' @param validation data frame from an independent validation cohort, including Y, T (if the validation cohort is from a PGx study), and G
#' @param corr genetic correlation matrix among K traits
#' @param pcut p-value cutoff for C+T method to construct individual stPRS
#' @param varcut variance cutoff to choose top principal components (PCs)
#' @param K number of traits
#' @param phenotype indicator of phenotype in target cohort, either "dis" or "pgx"
#' @details Construct mtPRS by performing ensemble learning
#' @return a list of fitted model, stPRSs and mtPRSs, and mtPRS-FIMEL
#' @references Zhai, S., Guo, B., and Shen, J., 2025. Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning.
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "CS", sparseness = "no",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.5, rho_C = 0.2,
#' K=4, m=2000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' re <- mtPRS_FIMEL(dat$base, dat$target, dat$validation, dat$corr, pcut = 1, varcut = 0.8, K = 4, phenotype = "pgx")
#' hist(re$slPRS_validation)
#' }
#'
mtPRS_FIMEL <- function(base, target, validation, corr, pcut=0.05, varcut=0.8, K=4, phenotype="pgx"){
  if(phenotype == "pgx"){G.target <- target[,-c(1:2)] %>% as.matrix(); Y.target <- target[,1] %>% as.vector(); T.target <- target[,2] %>% as.vector()}
  if(phenotype == "dis"){G.target <- target[,-1] %>% as.matrix(); Y.target <- target[,1] %>% as.vector()}

  print("mtPRS-PCA...")
  re1 <- mtPRS_PCA(base, target, validation, corr, pcut, varcut, K, phenotype)
  print("mtPRS-ML...")
  re2 <- mtPRS_ML(base, target, validation, pcut, K, phenotype)

  X.target <- cbind(re1$mtPRS_target, re2$mtPRS_target, re1$stPRS_target)
  X.validation <- cbind(re1$mtPRS_validation, re2$mtPRS_validation, re1$stPRS_validation)

  X.target <- X.target %>% scale(); X.validation <- X.validation %>% scale()

  X.target <- X.target %>% as.data.frame() %>% as.matrix(); X.validation <- X.validation %>% as.data.frame() %>% as.matrix()

  colnames(X.target) <- colnames(X.validation) <- c("mtPRS_PCA","mtPRS_ML",names(base))

  if(phenotype == "pgx"){G.validation <- validation[,-c(1:2)] %>% as.matrix(); Y.validation <- validation[,1] %>% as.vector(); T.validation <- validation[,2] %>% as.vector()}
  if(phenotype == "dis"){G.validation <- validation[,-1] %>% as.matrix(); Y.validation <- validation[,1] %>% as.vector()}

  print("Ensemble learning...")
  mod <- cv.glmnet(x = X.validation, y = Y.validation, alpha = 0.5)

  PRS.target <- predict(mod, newx = X.target, s = "lambda.min") %>% as.vector()
  PRS.validation <- predict(mod, newx = X.validation, s = "lambda.min") %>% as.vector()

  re <- list(mod=mod, indPRS_target=X.target, indPRS_validation=X.validation, slPRS_target=PRS.target, slPRS_validation=PRS.validation)
  re
}
