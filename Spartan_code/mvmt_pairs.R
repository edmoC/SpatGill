library(parallel)
set.seed(1)

h <- 1 # distance (m) between lattice points
dimGrid <- round(10/h) # 10m by 10m spatial domain
iN_H <- 10 # initial number of humans
iN_M <- 30 # initial number of mosquitoes
u_M <- 0.143 # 1wk mosquito lifespan
gamma_H <- 0.033 # 1mth infectious period 
T_HM <- 0.5 # infection P mosquito to human
T_MH <- 0.8 # infection P human to mosquito 
b <- 0.5 # biting rate
B <- b/(1-(1-1/dimGrid^2)^iN_H) # bite attempt rate

################################################################################
InitDomain <- function(parms) {
  with(parms, {
    # counts
    domainSH <- matrix(0L, ncol=dimGrid, nrow=dimGrid)
    domainIH <- matrix(0L, ncol=dimGrid, nrow=dimGrid)
    domainRH <- matrix(0L, ncol=dimGrid, nrow=dimGrid)
    domainSM <- matrix(0L, ncol=dimGrid, nrow=dimGrid)
    domainIM <- matrix(0L, ncol=dimGrid, nrow=dimGrid)
    # locations
    S_H_x <- numeric(iN_H)
    S_H_y <- numeric(iN_H)
    I_H_x <- numeric()
    I_H_y <- numeric()
    R_H_x <- numeric()
    R_H_y <- numeric()
    S_M_x <- numeric(iN_M-1)
    S_M_y <- numeric(iN_M-1)
    I_M_x <- numeric(1)
    I_M_y <- numeric(1)
    # humans
    n <- 0
    while (n < iN_H) {
      X <- sample(1:dimGrid, 1)
      Y <- sample(1:dimGrid, 1)
      domainSH[X,Y] <- domainSH[X,Y] + 1
      n <- n + 1
      S_H_x[n] <- X
      S_H_y[n] <- Y
    }
    # mosquitoes
    n <- 0
    while (n < iN_M-1) {
      X <- sample(1:dimGrid, 1)
      Y <- sample(1:dimGrid, 1)
      domainSM[X,Y] <- domainSM[X,Y] + 1
      n <- n + 1
      S_M_x[n] <- X
      S_M_y[n] <- Y
    }
    X <- sample(1:dimGrid, 1)
    Y <- sample(1:dimGrid, 1)
    domainIM[X,Y] <- 1
    I_M_x[1] <- X
    I_M_y[1] <- Y
    
    return(list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                domainSM=domainSM, domainIM=domainIM, S_H_x=S_H_x, S_H_y=S_H_y, 
                I_H_x=I_H_x, I_H_y=I_H_y, R_H_x=R_H_x, R_H_y=R_H_y, 
                S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y))
  })
}

human_recovery <- function(parms, domain, state) {
  with(as.list(c(parms, domain, state)), {
    selected <- sample(1:I_H, 1)
    X <- I_H_x[selected]
    Y <- I_H_y[selected]
    domainIH[X,Y] <- domainIH[X,Y] - 1
    domainRH[X,Y] <- domainRH[X,Y] + 1
    I_H_x <- I_H_x[-selected]
    I_H_y <- I_H_y[-selected]
    I_H <- I_H - 1
    R_H <- R_H + 1
    R_H_x[R_H] <- X
    R_H_y[R_H] <- Y

    domain <- list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                   domainSM=domainSM, domainIM=domainIM, 
                   S_H_x=S_H_x, S_H_y=S_H_y, I_H_x=I_H_x, I_H_y=I_H_y, 
                   R_H_x=R_H_x, R_H_y=R_H_y, 
                   S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y)
    state <- list(t=t, S_H=S_H, I_H=I_H, R_H=R_H, S_M=S_M, I_M=I_M, 
                  N_H=N_H, N_M=N_M)
    return(list(domain=domain, state=state))
  })
}

mosquito_death <- function(parms, domain, state) {
  with(as.list(c(parms, domain, state)), {
    type_M <- runif(1)
    if (type_M<S_M/N_M) { # susceptible mosquito
      selected <- sample(1:S_M, 1)
      X <- S_M_x[selected]
      Y <- S_M_y[selected]
      domainSM[X,Y] <- domainSM[X,Y] - 1  
      S_M_x <- S_M_x[-selected]
      S_M_y <- S_M_y[-selected]
      S_M <- S_M - 1
    } else { # infected mosquito 
      selected <- sample(1:I_M, 1)
      X <- I_M_x[selected]
      Y <- I_M_y[selected]
      domainIM[X,Y] <- domainIM[X,Y] - 1
      I_M_x <- I_M_x[-selected]
      I_M_y <- I_M_y[-selected]
      I_M <- I_M - 1
    }
    
    domain <- list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                   domainSM=domainSM, domainIM=domainIM, 
                   S_H_x=S_H_x, S_H_y=S_H_y, I_H_x=I_H_x, I_H_y=I_H_y, 
                   R_H_x=R_H_x, R_H_y=R_H_y, 
                   S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y)
    state <- list(t=t, S_H=S_H, I_H=I_H, R_H=R_H, S_M=S_M, I_M=I_M, 
                  N_H=N_H, N_M=N_M)
    return(list(domain=domain, state=state))
  })
}

mosquito_birth <- function(parms, domain, state) {
  with(as.list(c(parms, domain, state)), {
    X <- sample(1:dimGrid, 1)
    Y <- sample(1:dimGrid, 1)
    domainSM[X,Y] <- domainSM[X,Y] + 1
    S_M <- S_M + 1
    S_M_x[S_M] <- X
    S_M_y[S_M] <- Y
    
    domain <- list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                   domainSM=domainSM, domainIM=domainIM, 
                   S_H_x=S_H_x, S_H_y=S_H_y, I_H_x=I_H_x, I_H_y=I_H_y, 
                   R_H_x=R_H_x, R_H_y=R_H_y, 
                   S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y)
    state <- list(t=t, S_H=S_H, I_H=I_H, R_H=R_H, S_M=S_M, I_M=I_M, 
                  N_H=N_H, N_M=N_M)
    return(list(domain=domain, state=state))
  })
}

human_movement <- function(parms, domain, state) {
  with(as.list(c(parms, domain, state)), {
    type_H <- runif(1)
    if (type_H < S_H/N_H) { # susceptible human
      selected <- sample(1:S_H, 1)
      X <- S_H_x[selected]
      Y <- S_H_y[selected]
      # movement directions c(right, left, top, bottom)
      moves_X <- c(X%%dimGrid+1, (X-2)%%dimGrid+1, X, X)
      moves_Y <- c(Y, Y, Y%%dimGrid+1, (Y-2)%%dimGrid+1)
      dir <- sample(1:4, 1)
      X_move <- moves_X[dir]
      Y_move <- moves_Y[dir]
      # update domain
      domainSH[X,Y] <- domainSH[X,Y] - 1 
      domainSH[X_move,Y_move] <- domainSH[X_move,Y_move] + 1
      S_H_x[selected] <- X_move
      S_H_y[selected] <- Y_move
    } else if (type_H < (S_H+I_H)/N_H) { # infected human
      selected <- sample(1:I_H, 1)
      X <- I_H_x[selected]
      Y <- I_H_y[selected]
      # movement directions c(right, left, top, bottom)
      moves_X <- c(X%%dimGrid+1, (X-2)%%dimGrid+1, X, X)
      moves_Y <- c(Y, Y, Y%%dimGrid+1, (Y-2)%%dimGrid+1)
      dir <- sample(1:4, 1)
      X_move <- moves_X[dir]
      Y_move <- moves_Y[dir]
      # update domain
      domainIH[X,Y] <- domainIH[X,Y] - 1
      domainIH[X_move,Y_move] <- domainIH[X_move,Y_move] + 1
      I_H_x[selected] <- X_move
      I_H_y[selected] <- Y_move
    } else { # recovered human
      selected <- sample(1:R_H, 1)
      X <- R_H_x[selected]
      Y <- R_H_y[selected]
      # movement directions c(right, left, top, bottom)
      moves_X <- c(X%%dimGrid+1, (X-2)%%dimGrid+1, X, X)
      moves_Y <- c(Y, Y, Y%%dimGrid+1, (Y-2)%%dimGrid+1)
      dir <- sample(1:4, 1)
      X_move <- moves_X[dir]
      Y_move <- moves_Y[dir]
      # update domain
      domainRH[X,Y] <- domainRH[X,Y] - 1
      domainRH[X_move,Y_move] <- domainRH[X_move,Y_move] + 1
      R_H_x[selected] <- X_move
      R_H_y[selected] <- Y_move
    }
    
    domain <- list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                   domainSM=domainSM, domainIM=domainIM, 
                   S_H_x=S_H_x, S_H_y=S_H_y, I_H_x=I_H_x, I_H_y=I_H_y, 
                   R_H_x=R_H_x, R_H_y=R_H_y, 
                   S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y)
    state <- list(t=t, S_H=S_H, I_H=I_H, R_H=R_H, S_M=S_M, I_M=I_M, 
                  N_H=N_H, N_M=N_M)
    return(list(domain=domain, state=state))
  })
}

mosquito_movement <- function(parms, domain, state) {
  with(as.list(c(parms, domain, state)), {
    type_M <- runif(1)
    if (type_M < S_M/N_M) { # susceptible mosquito
      # movement
      selected <- sample(1:S_M, 1)
      X <- S_M_x[selected]
      Y <- S_M_y[selected]
      # movement directions c(right, left, top, bottom)
      moves_X <- c(X%%dimGrid+1, (X-2)%%dimGrid+1, X, X)
      moves_Y <- c(Y, Y, Y%%dimGrid+1, (Y-2)%%dimGrid+1)
      dir <- sample(1:4, 1)
      X_move <- moves_X[dir]
      Y_move <- moves_Y[dir]
      # update domain
      domainSM[X,Y] <- domainSM[X,Y] - 1
      domainSM[X_move,Y_move] <- domainSM[X_move,Y_move] + 1
      S_M_x[selected] <- X_move
      S_M_y[selected] <- Y_move
    } else { # infected mosquito
      # movement
      selected <- sample(1:I_M, 1)
      X <- I_M_x[selected]
      Y <- I_M_y[selected]
      # movement directions c(right, left, top, bottom)
      moves_X <- c(X%%dimGrid+1, (X-2)%%dimGrid+1, X, X)
      moves_Y <- c(Y, Y, Y%%dimGrid+1, (Y-2)%%dimGrid+1)
      dir <- sample(1:4, 1)
      X_move <- moves_X[dir]
      Y_move <- moves_Y[dir]
      # update domain
      domainIM[X,Y] <- domainIM[X,Y] - 1
      domainIM[X_move,Y_move] <- domainIM[X_move,Y_move] + 1
      I_M_x[selected] <- X_move
      I_M_y[selected] <- Y_move
    }
    
    domain <- list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                   domainSM=domainSM, domainIM=domainIM, 
                   S_H_x=S_H_x, S_H_y=S_H_y, I_H_x=I_H_x, I_H_y=I_H_y, 
                   R_H_x=R_H_x, R_H_y=R_H_y, 
                   S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y)
    state <- list(t=t, S_H=S_H, I_H=I_H, R_H=R_H, S_M=S_M, I_M=I_M, 
                  N_H=N_H, N_M=N_M)
    return(list(domain=domain, state=state))
  })
}

bite <- function(parms, domain, state) {
  with(as.list(c(parms, domain, state)), {
    type_M <- runif(1)
    type_H <- runif(1)
    success <- runif(1)
    if (type_M < S_M/N_M) { # susceptible mosquito
      selected <- sample(1:S_M, 1)
      X <- S_M_x[selected]
      Y <- S_M_y[selected]
      nSH <- domainSH[X,Y]
      nIH <- domainIH[X,Y]
      nRH <- domainRH[X,Y]
      if (nSH+nIH+nRH > 0) { # humans biteable
        if (type_H<nIH/(nSH+nIH+nRH) & success<T_MH) { # human to mosquito transmission
          domainSM[X,Y] <- domainSM[X,Y] - 1
          domainIM[X,Y] <- domainIM[X,Y] + 1
          S_M_x <- S_M_x[-selected]
          S_M_y <- S_M_y[-selected]
          S_M <- S_M - 1
          I_M <- I_M + 1
          I_M_x[I_M] <- X
          I_M_y[I_M] <- Y
        }
      }
    } else { # infected mosquito
      selected <- sample(1:I_M, 1)
      X <- I_M_x[selected]
      Y <- I_M_y[selected]
      nSH <- domainSH[X,Y]
      nIH <- domainIH[X,Y]
      nRH <- domainRH[X,Y]
      if (nSH+nIH+nRH > 0) { # humans biteable
        if (type_H<nSH/(nSH+nIH+nRH) & success<T_HM) { # mosquito to human transmission
          domainSH[X,Y] <- domainSH[X,Y] - 1
          domainIH[X,Y] <- domainIH[X,Y] + 1
          target <- which(S_H_x==X & S_H_y==Y)
          if (length(target) > 1) {
            target <- sample(target, 1)
          }
          S_H_x <- S_H_x[-target]
          S_H_y <- S_H_y[-target]
          S_H <- S_H - 1
          I_H <- I_H + 1
          I_H_x[I_H] <- X
          I_H_y[I_H] <- Y
        }
      }
    }
  
    domain <- list(domainSH=domainSH, domainIH=domainIH, domainRH=domainRH, 
                   domainSM=domainSM, domainIM=domainIM, 
                   S_H_x=S_H_x, S_H_y=S_H_y, I_H_x=I_H_x, I_H_y=I_H_y, 
                   R_H_x=R_H_x, R_H_y=R_H_y, 
                   S_M_x=S_M_x, S_M_y=S_M_y, I_M_x=I_M_x, I_M_y=I_M_y)
    state <- list(t=t, S_H=S_H, I_H=I_H, R_H=R_H, S_M=S_M, I_M=I_M, 
                  N_H=N_H, N_M=N_M)
    return(list(domain=domain, state=state))
  })
}

SpatGill <- function(parms) {
  with(parms, {
    domain <- InitDomain(parms)
    # simulation objects
    t_max <- 100
    state <- list(t=0, S_H=iN_H, I_H=0, R_H=0, S_M=iN_M-1, I_M=1, 
                  N_H=iN_H, N_M=iN_M)
    t_vec <- state$t
    S_H_vec <- state$S_H
    I_H_vec <- state$I_H
    R_H_vec <- state$R_H
    S_M_vec <- state$S_M 
    I_M_vec <- state$I_M
    # stepping through time
    while (state$t < t_max) {
      state$N_H <- state$S_H + state$I_H + state$R_H
      state$N_M <- state$S_M + state$I_M 
      R <- state$I_H*gamma_H + # human recovery
        state$N_M*u_M + # mosquito death
        state$N_M*u_M + # mosquito birth
        state$N_H*m_H + # human movement
        state$N_M*m_M + # mosquito movement
        state$N_M*B # bite
      step <- rexp(1, rate=R)
      # events
      event <- runif(1)
      if (event*R < state$I_H*gamma_H) { # human recovery
        out <- human_recovery(parms, domain, state)
        domain <- out$domain
        state <- out$state
      } else if (event*R < state$I_H*gamma_H+state$N_M*u_M) { # mosquito death
        out <- mosquito_death(parms, domain, state)
        domain <- out$domain
        state <- out$state
      } else if (event*R < state$I_H*gamma_H+2*state$N_M*u_M) { # mosquito birth
        out <- mosquito_birth(parms, domain, state)
        domain <- out$domain
        state <- out$state
      } else if (event*R < state$I_H*gamma_H+2*state$N_M*u_M+
                 state$N_H*m_H) { # human movement
        out <- human_movement(parms, domain, state)
        domain <- out$domain
        state <- out$state
      } else if (event*R < state$I_H*gamma_H+2*state$N_M*u_M+
                 state$N_H*m_H+state$N_M*m_M) { # mosquito movement
        out <- mosquito_movement(parms, domain, state)
        domain <- out$domain
        state <- out$state
      } else { # bite
        out <- bite(parms, domain, state)
        domain <- out$domain
        state <- out$state
      }
      # update
      state$t <- state$t + step
      t_vec <- c(t_vec, state$t)
      S_H_vec <- c(S_H_vec, state$S_H)
      I_H_vec <- c(I_H_vec, state$I_H)
      R_H_vec <- c(R_H_vec, state$R_H)
      S_M_vec <- c(S_M_vec, state$S_M)
      I_M_vec <- c(I_M_vec, state$I_M)
      # epidemic check
      if (state$I_M==0 & state$I_H==0) {
        return('dieout')
      }
    }
    return(list(t=t_vec, S_H=S_H_vec, I_H=I_H_vec, S_M=S_M_vec, I_M=I_M_vec))
  })
}

run_rep <- function(i, parms) {
  out <- SpatGill(parms)
  while(out=='dieout') {
    out <- SpatGill(parms)
  }
  # extracting values at fixed times 
  S_H_vec <- numeric(101)
  I_H_vec <- numeric(101)
  S_M_vec <- numeric(101)
  I_M_vec <- numeric(101)
  for (j in 0:100) {
    index <- max(which(out$t<=j))
    S_H_vec[j+1] <- out$S_H[index]
    I_H_vec[j+1] <- out$I_H[index]
    S_M_vec[j+1] <- out$S_M[index]
    I_M_vec[j+1] <- out$I_M[index]
  }
  return(data.frame(t=seq(0,100), S_H=S_H_vec, I_H=I_H_vec, 
                    S_M=S_M_vec, I_M=I_M_vec))
}

################################################################################
# Run in parallel
num_reps <- 50
nCores <- 50

# data.frames for average behaviour
S_H_DF <- data.frame(matrix(nrow=101, ncol=12))
colnames(S_H_DF) <- c('med_mh10mm1', 'top_mh10mm1', 'bot_mh10mm1',
                      'med_mh7mm2',  'top_mh7mm2', 'bot_mh7mm2',
                      'med_mh100mm10', 'top_mh100mm10', 'bot_mh100mm10',
                      'med_mh70mm20', 'top_mh70mm20', 'bot_mh70mm20')
I_H_DF <- data.frame(matrix(nrow=101, ncol=12))
colnames(I_H_DF) <- c('med_mh10mm1', 'top_mh10mm1', 'bot_mh10mm1',
                      'med_mh7mm2',  'top_mh7mm2', 'bot_mh7mm2',
                      'med_mh100mm10', 'top_mh100mm10', 'bot_mh100mm10',
                      'med_mh70mm20', 'top_mh70mm20', 'bot_mh70mm20')
S_M_DF <- data.frame(matrix(nrow=101, ncol=12))
colnames(S_M_DF) <- c('med_mh10mm1', 'top_mh10mm1', 'bot_mh10mm1',
                      'med_mh7mm2',  'top_mh7mm2', 'bot_mh7mm2',
                      'med_mh100mm10', 'top_mh100mm10', 'bot_mh100mm10',
                      'med_mh70mm20', 'top_mh70mm20', 'bot_mh70mm20')
I_M_DF <- data.frame(matrix(nrow=101, ncol=12))
colnames(I_M_DF) <- c('med_mh10mm1', 'top_mh10mm1', 'bot_mh10mm1',
                      'med_mh7mm2',  'top_mh7mm2', 'bot_mh7mm2',
                      'med_mh100mm10', 'top_mh100mm10', 'bot_mh100mm10',
                      'med_mh70mm20', 'top_mh70mm20', 'bot_mh70mm20')
# different combinations of movement parameters
mvmtParmsDF <- data.frame(m_H=c(10,7,100,70), m_M=c(1,2,10,20))
for (i in 1:dim(mvmtParmsDF)[1]) {
  m_H <- mvmtParmsDF$m_H[i]/(h^2)
  m_M <- mvmtParmsDF$m_M[i]/(h^2) 
  parms <- list(dimGrid=dimGrid, iN_H=iN_H, iN_M=iN_M, u_M=u_M,
                gamma_H=gamma_H, T_HM=T_HM, T_MH=T_MH, m_H=m_H, m_M=m_M, B=B)
  out <- mclapply(1:num_reps, run_rep, parms=parms, mc.cores=nCores)
  # data.frames to extract data into
  shDF <- data.frame(matrix(nrow=101, ncol=50))
  ihDF <- data.frame(matrix(nrow=101, ncol=50))
  smDF <- data.frame(matrix(nrow=101, ncol=50))
  imDF <- data.frame(matrix(nrow=101, ncol=50))
  # extract data
  for (j in 1:length(out)) {
    df <- out[[j]]
    shDF[,j] <- df$S_H
    ihDF[,j] <- df$I_H
    smDF[,j] <- df$S_M
    imDF[,j] <- df$I_M
  }
  # average behaviour
  S_H_DF[,((i-1)*3+1)] <- apply(shDF[,1:50], 1, median)
  I_H_DF[,((i-1)*3+1)] <- apply(ihDF[,1:50], 1, median)
  S_M_DF[,((i-1)*3+1)] <- apply(smDF[,1:50], 1, median)
  I_M_DF[,((i-1)*3+1)] <- apply(imDF[,1:50], 1, median)
  S_H_DF[,((i-1)*3+2)] <- apply(shDF[,1:50], 1, quantile, probs=0.9)
  I_H_DF[,((i-1)*3+2)] <- apply(ihDF[,1:50], 1, quantile, probs=0.9)
  S_M_DF[,((i-1)*3+2)] <- apply(smDF[,1:50], 1, quantile, probs=0.9)
  I_M_DF[,((i-1)*3+2)] <- apply(imDF[,1:50], 1, quantile, probs=0.9)
  S_H_DF[,((i-1)*3+3)] <- apply(shDF[,1:50], 1, quantile, probs=0.1)
  I_H_DF[,((i-1)*3+3)] <- apply(ihDF[,1:50], 1, quantile, probs=0.1)
  S_M_DF[,((i-1)*3+3)] <- apply(smDF[,1:50], 1, quantile, probs=0.1)
  I_M_DF[,((i-1)*3+3)] <- apply(imDF[,1:50], 1, quantile, probs=0.1)
}

save(S_H_DF, I_H_DF, S_M_DF, I_M_DF, file='movement_out.RData')





