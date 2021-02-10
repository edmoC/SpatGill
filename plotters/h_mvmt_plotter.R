library(ggplot2)
library(deSolve)

################################################################################
# ODE solution
# parameters
iN_H <- 10 # initial number of humans
iN_M <- 30 # initial number of mosquitoes
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

# expected bites
Ebites <- iN_M*b*100

################################################################################
# Deviations from ODE
biteErrorMat <- matrix(nrow=20, ncol=20)
ihErrorMat <- matrix(nrow=20, ncol=20)
for (D in 1:20) {
  for (res in 1:20) {
    # load output for single h and mvmt configuration
    filename <- paste0('outputs/h_mvmt_h10m30/D', D, 'res', res, '.RData')
    load(filename)
    # deviation from ODE
    biteErrorMat[D,res] <- abs(Ebites-medBite)
    ihErrorMat[D,res] <- sum((outODE$Y_H-medIH)^2)
    # plot 
    plt <- ggplot(data=outODE) +
      geom_line(aes(x=time, y=Y_H), color='red') +
      geom_line(aes(x=time, y=medIH))
    outfile <- paste0('figures/', 'D', D, 'res', res, '.png')
    ggsave(outfile, plt, width=5, height=3)
  }
}

################################################################################
# Error plot
# data.frame for plotting
D_vec <- numeric(400)
res_vec <- numeric(400)
biteError_vec <- numeric(400)
ihError_vec <- numeric(400)
i <- 0
for (D in 1:20) {
  for (res in 1:20) {
    i <- i + 1
    D_vec[i] <- D
    res_vec[i] <- res 
    biteError_vec[i] <- biteErrorMat[D,res]
    ihError_vec[i] <- ihErrorMat[D,res]
  }
}
plotDF <- data.frame(D=D_vec, res=res_vec, biteError=biteError_vec, 
                     ihError=ihError_vec)

# plotting
ggplot(plotDF) +
  geom_contour(aes(x=D, y=res, z=biteError, colour=..level..), bins=20)

ggplot(plotDF) +
  geom_contour(aes(x=D, y=res, z=ihError, colour=..level..), bins=50)

library(plotly)
fig <- plot_ly(z=ihErrorMat)
fig <- fig %>% add_surface()
fig <- fig %>% layout(
  scene=list(yaxis=list(title='Human diffusion'),
             xaxis=list(title='Domain resolution'),
             zaxis=list(title='Infected humans error'))
  )
fig