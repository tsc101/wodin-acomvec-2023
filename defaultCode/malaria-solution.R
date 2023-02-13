# human equations
deriv(S_h) <- - lambda_h * S_h + sigma * I_h
deriv(I_h) <- lambda_h * S_h - sigma * I_h

# mosquito equations
deriv(S_v) <- b_v * N_v - lambda_v * S_v - mu_v * S_v
deriv(I_v) <- lambda_v * S_v - mu_v * I_v

lambda_h = m * b * a * I_v / N_v  # human force of infection
lambda_v = c * a * I_h / N_h      # mosquito force of infection

# initial conditions
initial(S_h) <- N_h0 - I_init_h
initial(I_h) <- I_init_h

initial(S_v) <- N_v0 - I_init_v
initial(I_v) <- I_init_v

N_h <- S_h + I_h
N_v <- S_v + I_v

output(R0) <- (m * a^2 * b * c)/(sigma * mu_v)
output(EIR) <- m * a * (I_v / (S_v+ I_v)) * 365

# parameter values
N_h0 <- user(100)             # total human population size
N_v0 <- user(1000)            # total mosquito population size
I_init_h <- 1                 # number of infectious humans at start of epidemic
I_init_v <- 0                 # number of infectious mosquitoes at start of epidemic

b_v <- -log(p)                # mosquito birth rate
mu_v <- -log(p)               # mosquito death rate

sigma <- user(0.02)           # human recovery rate [0.005-0.05]

m = N_v / N_h                 # density of female mosquitoes per person. 1.93 is R0 = 1
a = user(0.33)                # biting rate per female mosquito
b = 0.2                       # probability of infection in susceptible human given bite from infectious mosquito
c = 0.05                      # probability of infection in susceptible mosquito given bite on an infectious human

p = user(0.9)                 # probability a mosquito survives one day
