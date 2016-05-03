##--------------------------------------------------------------------------------------------------------
## SCRIPT : Estimate bycatch from strandings with a Discrete Weibull Likelihood
##
## Authors : Matthieu Authier
## Last update : 2016-03-24
## R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
##--------------------------------------------------------------------------------------------------------

rm(list=ls())

### load relevant libraries
lapply(c("coda", "mvtnorm", "rstan", "dplyr", "ggplot2", "reshape", "colorspace", "DiscreteWeibull"), library, character.only=TRUE)

### change the working directory
WorkDir <- "D:/Dossier mathieu/Desktop/bycatch"
setwd(WorkDir)

### load the data from a text file
data <- read.table ("BoBWestChannel_bycatch.txt", header = TRUE,
                    colClasses = c ("character", rep ("numeric", 5)))

### have a look at the data
str (data)
summary (data)

### compute Nb bycaught and stranding proba per month/year
# color palette
cbPalette <- c (rev (heat_hcl (20,h = c (200, 260), 
                               c. = c(100, 30), 
                               l = c(30, 90), 
                               power = c(1/5, 1))),
                "black")
dodge <- position_dodge (width = 0.1)

### strandings time serie
stranded <-
   data %.% 
   group_by (Year,Month) %.% 
   summarise (Stranded = sum (deldel))

ggplot(data = stranded, 
       aes (x = Month, y = Stranded, group = as.character (Year), 
            colour = as.character (Year)
            ), 
       position = dodge) +
  geom_line (size = 1, position = dodge) +
  geom_point (size = 2.5, shape = 21, fill = "white", position = dodge) +
  scale_colour_manual (values = cbPalette[-21]) +
  ylab (quote (y[jt])) + xlab ("Month") +
  scale_x_continuous (breaks = 1:12,
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  # log10 scaling of the y axis (with visually-equal spacing)
  scale_y_log10 () +
  theme_set (theme_gray (base_size = 18))
#ggsave ("raw_data_obs.pdf", dpi = 600)

### probabilities from MOTHY computer experiments
prob <- 
  data %.% 
  group_by (Year, Month) %.% 
  summarise (pjt = mean ((mothy + 2) / (89 + 4)))

ggplot (data = prob,
        aes (x = Month, y = pjt, group = as.character (Year), 
             colour = as.character (Year)
             ),
        position = dodge) +
  geom_line (size = 1, position = dodge) +
  geom_point (size = 2.5, shape = 21, fill = "white", position = dodge) +
  scale_colour_manual (values = cbPalette[-21]) +
  ylab (quote (p[jt])) + xlab ("Month") +
  scale_x_continuous (breaks = 1:12,
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  scale_y_continuous (breaks = seq (0, 1, 0.1)) +
  theme_set (theme_gray (base_size = 18))
#ggsave("raw_data_proba.pdf", dpi = 600)

### Fit model for direct drift approach with Huang & Wand (2013) prior for Omega
# use great software stan

bycatch <- '
functions {
 real discrete_weibull_log (int y, real lambda, real omega){
  real q;
  real pmf;
  q <- 1 - inv_cloglog(lambda);
  pmf <- pow(q, pow(y, omega)) - pow(q, pow(y+1, omega)); 
  if (pmf > 0){
   return log(pmf);
  }
  else{
   return negative_infinity();
  }
 }
}

data {
 int<lower=1> n_year;
 int<lower=1> n_month;
 int<lower=1> n_obs;
 int<lower=1,upper=n_month> MONTH[n_obs];
 int<lower=1,upper=n_year> YEAR[n_obs];
 int<lower=0> MOTHY[n_obs];
 int<lower=0> STRANDINGS[n_obs];
 real<lower=1> DAYS[n_obs];
 vector[n_month] PRIOR_SIGMA2; // set to 1 for standard half-Cauchy
 int<lower=1> nu; // set to 1 for half-Cauchy
}

parameters {
 vector[2] intercept;
 vector[n_month] alpha;
 vector[n_month] epsilon;
 vector[n_year] beta;
 vector[n_month] scale;
 real<lower=0> sigma_month;
 real<lower=0> tau_month;
 real<lower=0> sigma_year;
 real logit_omega;
 vector[n_month] delta[n_year]; // 20 x 12 matrix
 cov_matrix[n_month] Omega;
}

transformed parameters{
 real omega;
 omega <- 2 * inv_logit(logit_omega);
}

model {
 logit_omega ~ student_t(7, 0.0, 5.0);
 intercept ~ student_t(7, 0.0, 5.0);
 // standard deviations
 sigma_month ~ cauchy(0.0, 1.0);
 tau_month ~ cauchy(0.0, 1.0);
 sigma_year ~ cauchy(0.0, 1.0);
 //month effects
 alpha ~ normal(0.0, sigma_month);
 epsilon ~ normal(0.0, tau_month);
 // year effects
 beta ~ normal(0.0, sigma_year);
 // half-t priors on the standard deviations
 // see Huang & Wand 2013 Bayesian Analysis
 for (j in 1:n_month){
  scale[j] ~ inv_gamma(0.5, 1/PRIOR_SIGMA2[j]);
 };
 Omega ~ inv_wishart(n_month+nu-1, diag_matrix(2*nu*scale));
 // month by year interactions
 for(i in 1:n_year){
  delta[i] ~ multi_normal(epsilon, Omega);
 };
 // likelihood
 for ( k in 1:n_obs ){
  MOTHY[k]+2 ~ binomial_logit(89+4, intercept[1] + alpha[MONTH[k]] + beta[YEAR[k]]);
  STRANDINGS[k] ~ discrete_weibull(intercept[2] + delta[YEAR[k], MONTH[k]] + log(DAYS[k]) - log10(), omega);
 };
}
'

hinit <- function (n_year) {
  Q <- runif (12, 0.1, 0.5) * diag (12)
  return (list (intercept = rnorm (2),
                sigma_year = runif (1, 0.1, 0.5),
                beta = rnorm (n_year, 0, 0.2),
                sigma_month = runif (1, 0.1, 0.5),
                tau_month = runif (1, 0.1, 0.5),
                alpha = rnorm (12, 0, 0.2),
                epsilon = rnorm (12, 0, 0.1),
                logit_omega = rnorm (1),
                delta = rmvnorm (n_year, rnorm (12), Q),
                Omega = Q,
                scale = diag (Q)
  ))
}

### generate initial values
my_inits <- vector(4, mode = 'list')
set.seed(19811021)
my_inits[[1]] <- hinit(20)
set.seed(314)
my_inits[[2]] <- hinit(20)
set.seed(81)
my_inits[[3]] <- hinit(20)
set.seed(19790428)
my_inits[[4]] <- hinit(20)

### fit model
system.time (
  fit <- stan (model_code = bycatch, 
               model_name = "Discrete Weibull Likelihood with Huang & Wand prior",
               data = with (data,
                           list (n_year = length (unique (Year)),
                                 n_month = length (unique (Month)),
                                 n_obs = nrow (data),
                                 STRANDINGS = deldel,
                                 MOTHY = mothy,
                                 MONTH = Month,
                                 YEAR = Year - min (Year) + 1,
                                 DAYS = NbDays,
                                 PRIOR_SIGMA2 = rep (1, 12),
                                 nu = 1# df for half-t prior
                                 )
              ),
              iter = 12000, 
              warmup = 2000, 
              chains = 4, 
              thin = 10,
              init = my_inits
  )
)

plot (fit)

sink (paste ("bycatch_dw_summary.txt", sep = ""))
print (fit, digit = 3)
sink ()

### plots with estimated total number of bycaught dolphins
get_ci <- function (x, alpha = 0.95) {
  return (c (HPDinterval (as.mcmc (x), prob = alpha)[1],
             mean (x, na.rm = TRUE),
             HPDinterval (as.mcmc (x), prob = alpha)[2])
          )
}

### compute mean of the discrete Weibull
dw_missed <- function (iter) {
  obs <- data.frame (STRANDINGS = data$deldel,
                     MOTHY = data$mothy,
                     MONTH = data$Month,
                     YEAR = data$Year-min(data$Year)+1,
                     DAYS = data$NbDays
  )
  
  p_float <- rbeta (1, 10.622, 47.674)
  
  intercept <- as.numeric(extract(fit, "intercept")$intercept[iter, ])
  delta <- as.matrix(extract(fit, "delta")$delta[iter, , ])
  alpha <- as.numeric(extract(fit, "alpha")$alpha[iter, ])
  beta <- as.numeric(extract(fit, "beta")$beta[iter, ])
  omega <- as.numeric(extract(fit, "omega")$omega[iter])
  
  obs$MISSED <- numeric(nrow(obs))
  obs$BYCATCH_TOT <- numeric(nrow(obs))
  
  for ( k in 1:nrow(obs) ){
    mu <- intercept[2] + delta[obs$YEAR[k], obs$MONTH[k]] + log(obs$DAYS[k]) - log(10)
    obs$MISSED[k] <- Edweibull (q = exp (-exp (mu)), beta = omega, zero = TRUE) * exp(-intercept[1] - alpha[obs$MONTH[k]] - beta[obs$YEAR[k]])
    obs$BYCATCH_TOT[k] <- (obs$STRANDINGS[k] + obs$MISSED[k])/p_float
  }
  
  tot <- obs %>% group_by (YEAR) %>% summarise (Stranded = sum(BYCATCH_TOT))
  return (as.numeric (tot$Stranded))
  
}

tot_bycatch <- t (sapply (1:4000, dw_missed))
### convert to a dataframe
get_ci <- function(x, alpha = 0.95){
  library(coda)
  return(c(HPDinterval(as.mcmc(x),prob=alpha)[1], mean(x, na.rm = TRUE), HPDinterval(as.mcmc(x),prob=alpha)[2]))
}
tot_bycatch <- as.data.frame(t(apply(tot_bycatch, 2, get_ci)))
names(tot_bycatch) <- c("Dead_lower", "Dead_mean", "Dead_upper")
tot_bycatch$Year <- 1990:2009

ggplot (data = tot_bycatch, 
        aes (x = as.numeric (Year), y = mean, col = "green")) +
  geom_line (size = 0.8) +
  geom_point (size = 3, shape = 21, fill = "white") +
  geom_line (aes (x = as.numeric (Year), y = lower), size = 0.6, alpha = 0.5) + 
  geom_line (aes (x = as.numeric (Year), y = upper), size = 0.6, alpha = 0.5) + 
  scale_fill_hue (c = 45, l = 80) +
  ylab ("Total Nb Bycaught Dolphins") + xlab ("Year") +
  # log10 scaling of the y axis (with visually-equal spacing)
  scale_y_log10 (breaks = c (1000, seq (2500, 15000, 2500))) +
  scale_x_continuous (breaks = c (1990:2009)) +
  guides (colour = FALSE) +
  theme_set (theme_gray (base_size = 18))
#ggsave("discreteweibull_bycatch.pdf", dpi = 600)
