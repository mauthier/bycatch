##--------------------------------------------------------------------------------------------------------
## SCRIPT : Estimate bycatch from strandings with a Negative Binomial Likelihood
##
## Authors : Matthieu Authier
## Last update : 2016-03-24
## R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
##--------------------------------------------------------------------------------------------------------

rm(list=ls())

### load relevant libraries
lapply(c("coda", "mvtnorm", "rstan", "dplyr", "ggplot2", "reshape", "colorspace"), library, character.only=TRUE)

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
 real<lower=0,upper=1> InvOverdispersion;
 vector[n_month] delta[n_year]; // 20 x 12 matrix
 cov_matrix[n_month] Omega;
}

transformed parameters {
 real omega;
 omega <- 1/InvOverdispersion;
}

model {
 real kappa;
 kappa <- 1/(omega-1);
 intercept ~ student_t(7,0,10);
 // standard deviations
 sigma_month ~ cauchy(0,1);
 tau_month ~ cauchy(0,1);
 sigma_year ~ cauchy(0,1);
 //month effects
 alpha ~ normal(0,sigma_month);
 epsilon ~ normal(0,tau_month);
 // year effects
 beta ~ normal(0,sigma_year);
 // half-t priors on the standard deviations
 // see Huang & Wand 2013 Bayesian Analysis
 for (j in 1:n_month){
  scale[j] ~ inv_gamma(0.5,1/PRIOR_SIGMA2[j]);
 };
 Omega ~ inv_wishart(n_month+nu-1,diag_matrix(2*nu*scale));
 // month by year interactions
 for(i in 1:n_year){
  delta[i] ~ multi_normal(epsilon,Omega);
 }
 // likelihood
 for ( k in 1:n_obs ){
  real lambda;
  MOTHY[k]+2 ~ binomial_logit(89+4,intercept[1] + alpha[MONTH[k]] + beta[YEAR[k]]);
  lambda <- DAYS[k]*exp(intercept[2] + delta[YEAR[k],MONTH[k]])/10;
  STRANDINGS[k] ~ neg_binomial(kappa,kappa/lambda);
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
                InvOverdispersion = runif (1, 0.2, 0.8),
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
               model_name = "Negative Binomial likelihood with Huang & Wand prior",
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

sink (paste ("bycatch_negbin_summary.txt", sep = ""))
print (fit, digit = 3)
sink ()

### plots with estimated total number of bycaught dolphins
get_ci <- function (x, alpha = 0.95) {
  return (c (HPDinterval (as.mcmc (x), prob = alpha)[1],
             mean (x),
             HPDinterval (as.mcmc (x), prob = alpha)[2])
          )
}

### predicted number of bycaught cetaceans at sea 
# that did not strand given that they float
Stranded_notStranded <- function (model, data) {
  Stranded <- array (NA, dim = c (4000, nrow (data)))
  notStranded <- array (NA, dim = c (4000, nrow (data)))
  for (i in 1:nrow (data)) {
    notStranded[, i] <- extract (model, "intercept")$intercept[, 2] + 
      extract (model, "delta")$delta[, data$Year[i] - min (data$Year) + 1, data$Month[i]] -
      (extract (model, "intercept")$intercept[, 1] + 
         extract (model, "alpha")$alpha[,data$Month[i]] + 
         extract (model, "beta")$beta[, data$Year[i] - min(data$Year) + 1]
       ) 
    Stranded[, i] <- extract (model, "intercept")$intercept[, 2] + 
      extract(model, "delta")$delta[, data$Year[i] - min (data$Year) + 1, data$Month[i]]
  }
  return (list (Stranded = Stranded, notStranded = notStranded))
}

system.time (pred <- Stranded_notStranded (model = fit, data)) 

system.time (observed <- with (pred,
                               sapply (1:nrow (data),
                                       function (i) {
                                         return (Stranded[, i] + log (data$NbDays[i] / 10))
                                         }
                                       )
                               )
             )

system.time (missed <- with (pred,
                             sapply (1:nrow (data),
                                     function (i) {
                                       return (notStranded[, i] + log (data$NbDays[i] / 10))
                                       }
                                     )
                             )
             )

missed_bycatch <- data.frame (Year = numeric (20),
                              y_obs = numeric (20),
                              y_miss_lower = numeric (20),
                              y_miss_mean = numeric (20),
                              y_miss_upper = numeric (20)
                              )
system.time (
  for(t in 1990:2009) {
    missed_bycatch[t-1989, ] <- c (t,
                                   sum (subset (data, Year == t)$deldel),
                                   round (get_ci (apply (exp (missed)[, with (data, which (Year == t))], 1, sum)), 0))
  }
)
rm (t)

### predicted number of bycaught cetaceans given that they float
lambda <- as.data.frame (matrix (NA, nrow = 20 * 12, ncol = 5))
names (lambda) <- c ("Year", "Month", "lower", "mean", "upper")
mu <- lambda
system.time(
  for(t in 1:20) {
    tt <- with (data, which (Year == 1989 + t))
    for(j in 1:12) {
      jt <- with (data,
                  which (Month == j))[with (data, which (Month == j)) %in% tt]
      lambda[12 * (t - 1) + j, ] <- c (1989 + t, j, round (get_ci (apply (exp (observed)[, jt], 1, mean)), 0))
      mu[12 * (t - 1) + j, ] <- c (1989 + t, j, round (get_ci (apply (exp (missed)[, jt], 1, mean)), 0))
    }
  }
)
#   user  system elapsed 
# 276.79   10.33  287.26
rm (t, tt, j, jt)

### some plots to visualize the results
ggplot (data = lambda, 
        aes (x = Month, y = mean, group = as.character(Year), 
             colour = as.character (Year)), 
        position = dodge) +
  geom_line (size = 0.8, position = dodge) +
  geom_point (size = 2.5, shape = 21, fill = "white", position = dodge) +
  scale_colour_manual (values = cbPalette[-21]) +
  ylab (quote (lambda[jt])) + xlab ("Month") +
  # log10 scaling of the y axis (with visually-equal spacing)
  scale_y_log10 () +
  scale_x_continuous (breaks = c(1:12),
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme_set (theme_gray (base_size = 18))
#ggsave("lambda.pdf", dpi = 600)

ggplot (data = mu, 
        aes (x = Month, y = mean, group = as.character(Year), 
             colour = as.character (Year)), 
        position = dodge) +
  geom_line (size = 0.8, position = dodge) +
  geom_point (size = 2.5, shape = 21, fill = "white", position = dodge) +
  scale_colour_manual (values = cbPalette[-21]) +
  ylab (quote (mu[jt])) + xlab ("Month") +
  # log10 scaling of the y axis (with visually-equal spacing)
  scale_y_log10 () +
  scale_x_continuous (breaks = c(1:12),
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme_set (theme_gray (base_size = 18))
#ggsave("mu.pdf", dpi = 600)

#### include floating proba
p_float <- function (n) { return (rbeta (n, 10.622, 47.674)) }
round (get_ci (p_float (4000)), 2)

# discovery probability: this beta distribution has Pr(X<0.80) = 0.025 et Pr(X<0.975) = 0.975
qbeta (c (0.025, 0.975), 36, 3.71)
p_discovery <- function (n) { return (rbeta (n, 36, 3.71)) }

# plot with p_discovery and p_float
#pdf(file="float_discovery.pdf",width=10,height=8)
par (cex = 2, mar = c(3, 3, 2, 1), mgp = c(2.0, 0.7, 0.0), tck = -0.01)
plot (c (0, 1), c (0, 10), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n',
      xlab = "probability", ylab = "density"
      )
axis (1, at = seq (0, 1, 0.1))
axis (2, at = seq (0, 10, 2), las = 1)
polygon (c (seq (0, 1, 0.001), seq (1, 0, -0.001)),
         c (dbeta (seq (0, 1, 0.001), 10.622, 47.674),
            rep (0, length (seq (0, 1, 0.001)))),
         border = NA, col = "lightblue")
lines (seq (0, 1, 0.001), dbeta (seq (0, 1, 0.001), 10.622, 47.674),
       col = "blue",lwd = 1)
text (0.175, 5, quote (p^float))

polygon (c (seq (0, 1, 0.001), seq (1, 0, -0.001)),
         c (dbeta (seq(0, 1, 0.001), 36, 3.71),
            rep (0, length (seq (0, 1, 0.001)))),
         border = NA, col = "lightgreen")
lines (seq (0, 1, 0.001), dbeta (seq (0, 1, 0.001), 36, 3.71), 
       col = "green", lwd = 1)
text(0.95, 5, quote (p^discovery))
#dev.off()

### predicted number of bycaught cetaceans (including those that did not float)
tot_bycatch <- sapply (1990:2009,
                       function (t) {
  apply ((subset (data, Year == t)$deldel + exp (missed[, which (data$Year==t)])) / p_float (4000), 1, sum)
                         }
  )
tot_bycatch <- as.data.frame (t (apply (tot_bycatch, 2, get_ci)))
names (tot_bycatch) <- c ("lower", "mean", "upper")
tot_bycatch$Year <- as.character (c (1990:2009))

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
#ggsave("negbin_bycatch.pdf", dpi = 600)

