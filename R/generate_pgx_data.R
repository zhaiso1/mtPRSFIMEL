#' Simulate PGx GWAS data for mtPRS analysis
#'
#' Simulate PGx GWAS data for mtPRS analysis
#' @param structure genetic correlation structure, either uniformly correlated (i.e., "CS" or "AR1") or clustered correlated ("clustered")
#' @param sparseness effect sparsity, either "no" or "half" or "more"
#' @param rho_DT effect correlations between traits in base cohort and the phenotype in target cohort
#' @param rho_T genetic correlation among traits in base cohort
#' @param rho_E correlation between prognostic and predictive effects in target cohort
#' @param rho_C between-cluster correlation when the genetic correlation structure is clustered correlated
#' @param K number of traits
#' @param m number of SNPs
#' @param pcausal proportion of causal SNPs
#' @param blocksize LD block size (i.e., the number of SNPs in each LD block)
#' @param gamma prognostic-to-predictive effect size ratio
#' @param samplesize sample size of target cohort
#' @param h2_base heritability in base cohort
#' @param h2_target heritability in target cohort
#' @details Simulate PGx GWAS data for mtPRS analysis
#' @return a list of disease GWAS summary statistics in base cohort, individual-level data in validation and target cohorts, genetic correlation matrix among traits, and true effect sizes in base and target cohorts
#' @references Zhai, S., Guo, B., and Shen, J., 2025. Improving multi-trait polygenic risk score prediction using fine-mapping and ensemble learning.
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "CS", sparseness = "no",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.5, rho_C = 0.2,
#' K=4, m=2000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' plot(dat$truesize[,1])
#' print(dat$corr)
#' }
#'
generate_pgx_data <- function(structure="CS", sparseness="no", rho_DT, rho_T, rho_E, rho_C, K=4, m=1000, pcausal=0.1, blocksize=100, gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3){
  map <- fread("https://github.com/zhaiso1/mtPRSFIMEL/blob/main/data/map.csv?raw=true") %>% as.data.frame()
  snpnames <- paste0(map$chromosome,":",map$position,"_",map$snp.name)

  print("Step 1. Prepare correlation matrix")
  Omega <- simu_pgx_cor(structure, rho_E, rho_C, rho_T, rho_DT, K)
  Sigma <- h2_base/(m*pcausal)*Omega
  Sigma <- nearPD(Sigma)$mat %>% as.matrix()

  print("Step 2. Simulate true effect size")
  TRUESIZE <- simu_pgx_truesize(m, blocksize, pcausal, Sigma, K, sparseness)

  print("Step 3. Simulate base cohort")
  se_pool <- fread("https://github.com/zhaiso1/mtPRSFIMEL/blob/main/data/se.csv?raw=true") %>% as.matrix()
  if(K > 4){
    add <- sample(as.vector(se_pool), m*(K-4), replace = TRUE)
    add <- matrix(add, nrow = m)
    se_pool <- cbind(se_pool, add)
    colnames(se_pool) <- paste0("Trait",1:K)
  }
  if(K < 4){
    se_pool <- se_pool[,1:K]
  }
  base <- simu_pgx_base(TRUESIZE, se_pool, K, m, map)

  print("Step 4. Simulate validation and target cohort")
  G <- fread("https://github.com/zhaiso1/mtPRSFIMEL/blob/main/data/geno.csv?raw=true") %>% as.matrix()
  colnames(G) <- snpnames
  target <- simu_pgx_target(TRUESIZE, G, gamma, m, samplesize, h2_target, K)

  validation <- target[1:400,]
  target <- target[401:nrow(target),]

  print("Step 5. Output results")
  xx <- t(TRUESIZE)
  rownames(xx) <- snpnames; colnames(xx) <- c(paste0("mu",1:K),"beta","alpha")
  re <- list(base=base, validation=validation, target=target, corr=Omega, truesize=xx)
  re
}

simu_pgx_cor <- function(structure, rho_E, rho_C, rho_T, rho_DT, K){
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
    rho^exponent
  }

  cs_cor <- function(n, rho){
    re <- matrix(rho, ncol = n, nrow = n); re <- re + diag(1-rho, nrow = n, ncol = n)
    re
  }

  if(structure == "CS"){
    corr11 <- cs_cor(K, rho_T)
  }
  if(structure == "AR1"){
    corr11 <- ar1_cor(K, rho_T)
  }
  if(structure == "clustered"){
    p11 <- p22 <- cs_cor(K/2, rho_T)
    p12 <- p21 <- matrix(rep(rho_T*rho_C,(K/2)^2),ncol=K/2)
    corr11 <- rbind(cbind(p11, p12), cbind(p21, p22))
  }

  corr12 <- cbind(rho_DT, rho_DT)
  corr21 <- rbind(rho_DT, rho_DT)
  corr22 <- matrix(rho_E, ncol = 2, nrow = 2); corr22 <- corr22 + diag(1-rho_E, nrow = 2, ncol = 2)
  corr <- rbind(cbind(corr11, corr12), cbind(corr21, corr22))

  colnames(corr) <- rownames(corr) <- NULL
  corr
}

simu_pgx_truesize <- function(m, blocksize, pcausal, Sigma, K, sparseness){
  nblock <- ceiling(m/blocksize)
  pi <- rbeta(nblock, pcausal, 1-pcausal)

  TRUESIZE <- double()
  for (k in 1:nblock) {
    for (j in 1:blocksize) {
      label <- runif(1)
      if(label <= pi[k]){truesize <- rmvnorm(1, mean = rep(0,K+2), sigma = Sigma) %>% as.vector()}
      if(label > pi[k]){truesize <- rep(0,K+2)}

      if(sparseness == "half"){phi <- c(rbinom(K, 1, 0.5), 1, 1); truesize <- phi*truesize}
      if(sparseness == "more"){phi <- c(rbinom(K, 1, 0.25), 1, 1); truesize <- phi*truesize}

      TRUESIZE <- cbind(TRUESIZE, truesize)
    }
  }

  return(TRUESIZE)
}

simu_pgx_base <- function(TRUESIZE, se_pool, K, m, map){
  SumStat <- list()
  for (k in 1:K) {
    mat <- matrix(NA, ncol=3, nrow=m)
    for (j in 1:m) {
      mat[j,1] <- rnorm(1, mean = TRUESIZE[k,j], sd = se_pool[j,k])
    }
    mat[,2] <- se_pool[,k]
    mat[,3] <- 2*(1-pnorm(abs(mat[,1]/mat[,2])))
    mat <- as.data.frame(mat)

    mat <- cbind.data.frame(map[,c(1,2,4,5,6)], mat)

    colnames(mat) <- c("chr","rsid","bp","A1","A2","beta","se","p"); rownames(mat) <- NULL

    SumStat[[k]] <- mat
  }
  names(SumStat) <- colnames(se_pool)
  return(SumStat)
}

simu_pgx_target <- function(TRUESIZE, G, gamma, m, samplesize, h2_target, K){
  index <- sample(1:nrow(G), samplesize) %>% sort()
  G_sub <- G[index,]
  nsubj <- nrow(G_sub)
  beta <- TRUESIZE[K+1,]
  alpha <- gamma*TRUESIZE[K+2,]
  Tr <- rbinom(nsubj, 1, 0.5)
  main_part <- (G_sub%*%beta + (Tr*G_sub)%*%alpha) %>% as.matrix() %>% as.vector()
  sigma <- ((1-h2_target)/h2_target*var(main_part) - var(Tr)) %>% sqrt()
  Y <- Tr + main_part + rnorm(nsubj, mean = 0, sd=sigma)

  ind <- cbind.data.frame(Y=Y, Tr=Tr, G_sub)
  return(ind)
}


