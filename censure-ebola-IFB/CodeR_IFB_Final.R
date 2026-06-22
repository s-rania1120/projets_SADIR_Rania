setwd("/Users/rania/Desktop/Travail de groupe-20251117")
data <- read.csv("Données Ebola.txt", sep = "")

### QUESTION 1 

y    <- data$ct
time <- data$time

mod0 <- lm(y ~ time)
summary(mod0)

beta0_mod0 <- coef(mod0)[1]
beta1_mod0 <- coef(mod0)[2]
sigma_mod0 <- summary(mod0)$sigma

beta0_mod0  
beta1_mod0 
sigma_mod0 


confint(mod0)


### QUESTION 2

C <- 40

# Log-vraisemblance censurée
f_loglik_inv <- function(theta, y, time, C)
{
  beta0 <- theta[1]
  beta1 <- theta[2]
  sigma <- theta[3]
  mu <- beta0 + beta1 * time
  
  obs  <- (y < C)
  lVr_obs <- sum( dnorm(y[obs], mean=mu[obs], sd=sigma, log=TRUE) )
  
  cens <- (y == C)
  zc   <- (C - mu[cens]) / sigma
  surv <- 1 - pnorm(zc)
  
  lVr_cens <- sum(log(surv))
  
  return(-(lVr_obs + lVr_cens))
}

# Valeurs initiales : celles du modèle sans censure
theta_init <- c(beta0_mod0, beta1_mod0, sigma_mod0)

# Maximisation
nlm <- nlm(f_loglik_inv,
           p = theta_init,
           y = y, time = time, C = C
)

beta0_cens  <- nlm$estimate[1]
beta1_cens  <- nlm$estimate[2]
sigma_cens  <- nlm$estimate[3]

beta0_cens   
beta1_cens  
sigma_cens  

### QUESTION 3

n <- length(time)

beta0_true  <- beta0_cens
beta1_true  <- beta1_cens
sigma_true  <- sigma_cens


# Estimation sans censure 
estimer_noncensuree <- function(y_obs, time) {
  mod <- lm(y_obs ~ time)
  beta0_hat <- coef(mod)[1]
  beta1_hat <- coef(mod)[2]
  sigma_hat <- summary(mod)$sigma
  return(c(beta0_hat, beta1_hat, sigma_hat))
}

# Estimation avec censure 
estimer_censure <- function(y_obs, time, C, theta_init) {
  res <- nlm(f_loglik_inv,
             p    = theta_init,
             y    = y_obs,
             time = time,
             C    = C,
             hessian = FALSE)
  return(res$estimate)  # (beta0, beta1, sigma)
}


## Simulation pour un seuil de censure donné  

design_ebola <- function(C, theta_init) {
  # 1) Génération des valeurs latentes de Ct :
  mu_true <- beta0_true + beta1_true * time
  y_true  <- rnorm(n, mean = mu_true, sd = sigma_true)
  
  # 2) Application de la censure à droite au seuil C :
  y_obs <- pmin(y_true, C)
  
  # 3) Estimation avec le modèle non censuré :
  est_noncens <- estimer_noncensuree(y_obs, time)
  
  # 4) Estimation avec le modèle censuré :
  est_cens <- estimer_censure(y_obs, time, C, theta_init)
  
  # 5) Réunir les 6 estimateurs dans un seul vecteur :
  out <- c(est_noncens, est_cens)
  names(out) <- c("beta0_noncens", "beta1_noncens", "sigma_noncens",
                  "beta0_cens",    "beta1_cens",    "sigma_cens")
  return(out)
}


## Boucle de simulation 

set.seed(123)               
K      <- 1000              
C_vals <- c(40, 38, 35)     

theta_init_sim <- c(beta0_true, beta1_true, sigma_true)

simuls_list <- list()

for (C_sim in C_vals) {
  simuls_C <- replicate(K, design_ebola(C = C_sim,
                                        theta_init = theta_init_sim))
  simuls_list[[as.character(C_sim)]] <- simuls_C
}


rowMeans(simuls_list[["40"]])
apply(simuls_list[["40"]], 1, sd)

rowMeans(simuls_list[["38"]])
apply(simuls_list[["38"]], 1, sd)

rowMeans(simuls_list[["35"]])
apply(simuls_list[["35"]], 1, sd)

