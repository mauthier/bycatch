rm(list=ls())
library(MCMCpack)

n_sim <- 10000
nu <- 12

### inverse Wishart prior
wish_prior <- array(NA,dim=c(n_sim,nu,nu))
for(i in 1:n_sim){
  wish_prior[i,,] <- riwish(nu+1, S=diag(nu))
}

### cholesky prior
chol_mat <- function(n_sim,nu){
  Da <- array(0,dim=c(n_sim,nu,nu))
  La <- Da
  M <- Da
  for ( i in 1:nu ) {
    Da[,i,i] <- abs(rnorm(n_sim,0,1.5))
    La[,i,i] <- 1
  }
  for( i in 1:(nu-1) ){
    for( j in (i+1):nu ) {La[,j,i] <- rnorm(n_sim,0,0.5)}
  }
 for ( k in 1:n_sim ) { M[k,,] <- Da[k,,]%*%La[k,,]%*%t(La[k,,])%*%Da[k,,]}
 return(M)
}
chol_prior <- chol_mat(n_sim=n_sim,nu=nu)

### barnard prior
barnard_mat <- function(n_sim,nu){
  M <- array(NA,dim=c(n_sim,nu,nu))
  for ( k in 1:n_sim ) { 
    Da <- exp(rnorm(nu,-0.5,1))*diag(nu)
    draw <- riwish(nu+1, S=diag(nu))
    Dq <- 1/sqrt(diag(draw))*diag(nu)
    M[k,,] <- Da%*%Dq%*%draw%*%Dq%*%Da
  }
  return(M)
}
barn_prior <- barnard_mat(n_sim=n_sim,nu=nu)

### huang prior
huang_mat <- function(n_sim,nu){
  M <- array(NA,dim=c(n_sim,nu,nu))
  for ( k in 1:n_sim ) { 
    D <- rgamma(nu,shape=0.5,scale=1)*diag(12)
    M[k,,] <- riwish(nu, S=2*D)
  }
  return(M)
}
huang_prior <- huang_mat(n_sim=n_sim,nu=nu)

### plot
f <- function(mat,param="c",i,j=0){
  if(param=="v"){
    x <- sqrt(mat[,i,i])
    return(data.frame(data.frame(x=density(x,kernel="g",to=10)$x,
                                 y=density(x,kernel="g",to=10)$y)
    )
    )
  }
  else{
    x <- mat[,i,j]/sqrt(mat[,i,i]*mat[,j,j])
    return(data.frame(data.frame(x=density(x,kernel="g")$x,
                                 y=density(x,kernel="g")$y)
    )
    )
  }

}

dd <- f(mat=wish_prior,param="v",i=1)
dd$prior <- "Inverse Wishart"
dd$param <- "Std Deviation"
ddw <- f(mat=wish_prior,i=1,j=2)
ddw$prior <- "Inverse Wishart"
ddw$param <- "Correlation"
dd <- rbind(dd,ddw)

ddw <- f(mat=chol_prior,param="v",i=1)
ddw$prior <- "Cholesky"
ddw$param <- "Std Deviation"
dd <- rbind(dd,ddw)
ddw <- f(mat=chol_prior,i=1,j=2)
ddw$prior <- "Cholesky"
ddw$param <- "Correlation"
dd <- rbind(dd,ddw)

ddw <- f(mat=barn_prior,param="v",i=1)
ddw$prior <- "Barnard"
ddw$param <- "Std Deviation"
dd <- rbind(dd,ddw)
ddw <- f(mat=barn_prior,i=1,j=2)
ddw$prior <- "Barnard"
ddw$param <- "Correlation"
dd <- rbind(dd,ddw)

ddw <- f(mat=huang_prior,param="v",i=1)
ddw$prior <- "Huang"
ddw$param <- "Std Deviation"
dd <- rbind(dd,ddw)
ddw <- f(mat=huang_prior,i=1,j=2)
ddw$prior <- "Huang"
ddw$param <- "Correlation"
dd <- rbind(dd,ddw)
rm(ddw)

library(ggplot2)
cbPalette <- c("black","red","green","blue")

ggplot(subset(dd,param=="Correlation"), aes(x=x, y=y, colour=prior)) + 
  geom_line(size=0.8) + 
  ylab("Density") + xlab("Correlation") + 
  scale_colour_manual(values=cbPalette) +
  #facet_wrap(~ param, ncol = 2) +
  theme_set(theme_gray(base_size=18)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave("prior_cor.pdf",width=10,height=8)

ggplot(subset(dd,param=="Std Deviation"), aes(x=x, y=y, colour=prior)) + 
  geom_line(size=0.8) + 
  ylab("Density") + xlab("Standard deviation") + 
  scale_colour_manual(values=cbPalette) +
  #facet_wrap(~ param, ncol = 2) +
  theme_set(theme_gray(base_size=18)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave("prior_sd.pdf",width=10,height=8)