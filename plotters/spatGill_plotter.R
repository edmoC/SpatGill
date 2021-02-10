library(deSolve)
library(ggplot2)

load('outputs/spatGill_h100m300mh10mm1.RData')

# data.frames to store extract data into
shDF <- data.frame(matrix(nrow=101, ncol=50))
shDF['t'] <- seq(0,100)
ihDF <- data.frame(matrix(nrow=101, ncol=50))
ihDF['t'] <- seq(0,100)
smDF <- data.frame(matrix(nrow=101, ncol=50))
smDF['t'] <- seq(0,100)
imDF <- data.frame(matrix(nrow=101, ncol=50))
imDF['t'] <- seq(0,100)
# extract data
for (i in 1:length(out)) {
  df <- out[[i]]
  shDF[,i] <- df$S_H
  ihDF[,i] <- df$I_H
  smDF[,i] <- df$S_M
  imDF[,i] <- df$I_M
}

# average behaviour
shDF['median'] <- apply(shDF[,1:50], 1, median)
ihDF['median'] <- apply(ihDF[,1:50], 1, median)
smDF['median'] <- apply(smDF[,1:50], 1, median)
imDF['median'] <- apply(imDF[,1:50], 1, median)

shDF['bot'] <- apply(shDF[,1:50], 1, quantile, probs=0.1)
ihDF['bot'] <- apply(ihDF[,1:50], 1, quantile, probs=0.1)
smDF['bot'] <- apply(smDF[,1:50], 1, quantile, probs=0.1)
imDF['bot'] <- apply(imDF[,1:50], 1, quantile, probs=0.1)

shDF['top'] <- apply(shDF[,1:50], 1, quantile, probs=0.9)
ihDF['top'] <- apply(ihDF[,1:50], 1, quantile, probs=0.9)
smDF['top'] <- apply(smDF[,1:50], 1, quantile, probs=0.9)
imDF['top'] <- apply(imDF[,1:50], 1, quantile, probs=0.9)


# Comparison with ODE
# parameters
iN_H <- 100 # initial number of humans
iN_M <- 300 # initial number of mosquitoes
u_M <- 0.143 # 1wk mosquito lifespan
gamma_H <- 0.033 # 1mth infectious period 
T_HM <- 0.5 # infection P mosquito to human
T_MH <- 0.8 # infection P human to mosquito 
u_H <- 5.5e-5 # 50yr human lifespan
v_H <- iN_H * u_H # for constant human population
v_M <- iN_M * u_M # for constant mosquito population
b <- 0.5
r <- b/iN_H

# setup
parms <- list(N_H=iN_H, N_M=iN_M, u_H=u_H, u_M=u_M, v_H=v_H, v_M=v_M, 
              gamma_H=gamma_H, T_HM=T_HM, T_MH=T_MH, b=b, r=r)
initial <- c(X_H=iN_H, Y_H=0, X_M=iN_M-1, Y_M=1)
times <- seq(0, 100)
sys <- function(t, y, parms) {
  with(as.list(y, parms), {
    dX_H = v_H - r*T_HM*Y_M*X_H - u_H*X_H
    dY_H = r*T_HM*Y_M*X_H - u_H*Y_H - gamma_H*Y_H
    dX_M = v_M - r*T_MH*Y_H*X_M - u_M*X_M
    dY_M = r*T_MH*Y_H*X_M - u_M*Y_M
    
    list(c(dX_H, dY_H, dX_M, dY_M))
  })
}

# solve ODE
outODE <- ode(initial, times, sys, parms)
outODE <- data.frame(outODE)

shplot <- ggplot() + 
  geom_line(data=outODE, aes(x=time, y=X_H, color='ODE')) + 
  theme_minimal() +
  geom_line(data=shDF, aes(x=t, y=median, color='simulation median')) +
  geom_ribbon(data=shDF, aes(x=t, ymin=bot, ymax=top, fill='simulation 80% CI'), alpha=0.3) +
  labs(x='Time (days)', y='Susceptible Humans', color='', fill='') +
  scale_color_manual(values=c('blue', 'black')) +
  scale_fill_manual(values=c('grey'))

ihplot <- ggplot() + 
  geom_line(data=outODE, aes(x=time, y=Y_H, color='ODE')) + 
  theme_minimal() +
  geom_line(data=ihDF, aes(x=t, y=median, color='simulation median')) +
  geom_ribbon(data=ihDF, aes(x=t, ymin=bot, ymax=top, fill='simulation 80% CI'), alpha=0.3) +
  labs(x='Time (days)', y='Infected Humans', color='', fill='') +
  scale_color_manual(values=c('red', 'black')) +
  scale_fill_manual(values=c('grey'))

smplot <- ggplot() + 
  geom_line(data=outODE, aes(x=time, y=X_M, color='ODE')) + 
  theme_minimal() +
  geom_line(data=smDF, aes(x=t, y=median, color='simulation median')) +
  geom_ribbon(data=smDF, aes(x=t, ymin=bot, ymax=top, fill='simulation 80% CI'), alpha=0.3) +
  labs(x='Time (days)', y='Susceptible Mosquitoes', color='', fill='') +
  scale_color_manual(values=c('purple', 'black')) +
  scale_fill_manual(values=c('grey'))

implot <- ggplot() + 
  geom_line(data=outODE, aes(x=time, y=Y_M, color='ODE')) + 
  theme_minimal() +
  geom_line(data=imDF, aes(x=t, y=median, color='simulation median')) +
  geom_ribbon(data=imDF, aes(x=t, ymin=bot, ymax=top, fill='simulation 80% CI'), alpha=0.3) +
  labs(x='Time (days)', y='Infected Mosquitoes', color='', fill='') +
  scale_color_manual(values=c('green', 'black')) +
  scale_fill_manual(values=c('grey'))

shplot
ihplot
smplot
implot

ggsave('figures/SusceptibleHumans.pdf', shplot, 
       width=10, height=5)
ggsave('figures/InfectedHumans.pdf', ihplot, 
       width=10, height=5)
ggsave('figures/SusceptibleMosquitoes.pdf', smplot, 
       width=10, height=5)
ggsave('figures/InfectedMosquitoes.pdf', implot, 
       width=10, height=5)
