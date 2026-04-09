
#These scripts run several functions related to a Virtual Population Analysis of the Scalloped Hammerheads (Sphyrna lewini) in the Eastern Tropical Pacific. 

#This model is desribed in Gomez-Garcia et al (in press) "Virtual Population Analysis of the critically endangered Scalloped Hammerhead (Sphyrna lewini) in the Eastern Tropical Pacific" 

#This code has not been optimized for speed

#Created by Miguel de Jesus Gomez Garcia
#Created: 20-November-2025
#Last edited: 09-Apr-2026
#First uploaded to GitHub on 07-Apr-2026


# Load packages ------------------------------------------------------------

library(doParallel)
library(R2jags)
library(MCMCvis)
library(tidyverse)
library(viridis)   #color pallettes
library(lsmeans)
library(broom)
library(popdemo)    #easy functions for elasticity
library(patchwork)  #multiple plots in single figure

#Set seed. 
set.seed(2026)


# Parallel Computing ------------------------------------------------------

#setup parallel backend to use many processors
cores=detectCores()

cl <- makeCluster(cores/2) #half of your cores not to overload your computer. But you can go higher
#

registerDoParallel(cl)




# Base population Parameters from literature ---------------------------------------------------

# Growth rate k or g 
gs <- mean(c(0.1,0.12,0.161,0.155,0.142,0.165,0.119))              
g_f <- mean(c(0.1,0.161,0.142))       
g_m <-mean(c(0.12,0.155,0.165))  

# Asymptotic growth length from literature
Linf <- mean(c(376,364,289.6,259.8,367.9,317.7,285.51))             

Linf_f <- mean(c(376,289.6,367.9)) 
Linf_m <- mean(c(364,259.8,317.7))

# Length at birth L0 from literature
L0 <- mean(c(53.2,56.8,51,48.5,41.1))             

L0_f <- mean(c(53.2,48.5)) 
L0_m <- mean(c(56.8,48.5))

#Lm_f1 from Alejo plata and Lm_f2 from Estupiñan-montaño et al. 2021  

Lm_f1<-160
Lm_m1<- 155

Lm_f2<-219.4
Lm_m2<- 178.1

At_f1<-log((Linf_f-L0_f)/(Linf_f-Lm_f1))/gs

#Dont have t0. Fishbase t0 -1.16 f, 1.18m. Mexico
At<- -1.16 - (1 / gs) * log(1 - (Lm_f2/ Linf))

At<-L0_f - (1 / gs) * log(1 - (Lm_f1 / Linf_f))

amat <- mean(c(13,8))    # Age at maturity both sexes, literature
amat_f <- 13
amat_m <- 8 

#Fecundity / litter size. Average from literature

f_lit<- mean(c(16,22,25.8,14,41))
sd(c(16,22,25.8,14,41))
#Longevity. 

#Longevity from Formula
amax <- log(
  ((Linf*0.99)*(Linf*L0))/ (L0*(Linf-((Linf*0.99))))
) /gs 

#From vertebrae varies from 11 to 35. Estimations go up to 55. We averaged using these values
amax <- mean(c(12.5,11,35,24.4,21.1)) 

# NaturalMortality calculation Ms

# Define the list of natural mortality estimates
Frisk1 <- exp(0.42 * log(gs) - 0.83)                  # Frisk et al. (2001)
Frisk2 <- 1 / (0.44 * amat + 1.87)                   # Frisk et al. (2001)
Hisano1 <- 1.65 / amat                               # Hisano et al. (2011)
Then1 <- 4.899 * amax^(-0.916)                       # Then et al. (2015)
Then2 <- (4.118 * gs^0.73) * Linf^(-0.33)            # Then et al. (2015)
Dureuil1 <- exp(1.583 - 1.087 * log(amax))           # Dureuil et al. (2021)
Dureuil2 <- -log(0.0178) / amax                      # Dureuil et al. (2021)

#Averaged Natural mortality

M_lit <- mean(c(Frisk1,Frisk2,Hisano1,Then1,Then2,Dureuil1,Dureuil2))

sd(c(Frisk1,Frisk2,Hisano1,Then1,Then2,Dureuil1,Dureuil2))

# Base Model Stage-specific mortalities------------------------------------------------------------------

sink("MatrixModel_Base.txt")
cat("
model {

###PRIORS-----------------------------------
  
  
#Stochastic random mortalities for all stages



#set sd for natural mortality
tau<- 1/(0.05^2) 
  
#Set stage specific mortality priors

 #High neonate mortality Between 50 and 95% Low but existing fishing mortality
 
 M[1] ~ dunif(0.4,0.8) 
 Fs[1] ~ dunif(0.1,0.15)
 
 #Juvenile natural mortality higher than adult mortality. fishing mortality same as adults
 M[2] ~ dnorm(M_lit*2, tau*2) 
 Fs[2] ~ dunif(0.1, 0.4)
 
  #Adult and subadult natural mortalities from literature, diffuse for fishing.
  
 M[3] ~ dnorm(M_lit, tau) 
 Fs[3] ~ dunif(0.1, 0.4)
 
 M[4] ~ dnorm(M_lit, tau) 
 Fs[4] ~ dunif(0.1, 0.4)
 
  #Higher juvenile female natural mortality. same range of potential fishing mortality as adults 
  
 M[5] ~ dnorm(M_lit*2, tau*2) 
 Fs[5] ~ dunif(0.1, 0.4)
 
   #Adult and subadult natural mortalities from literature, diffuse for fishing
   
 M[6] ~ dnorm(M_lit, tau) 
 Fs[6] ~ dunif(0.1, 0.4)
 
 M[7] ~ dnorm(M_lit, tau) 
 Fs[7] ~ dunif(0.1, 0.4)
 
 M[8] ~ dnorm(M_lit, tau) 
 Fs[8] ~ dunif(0.1, 0.4)
 
 #Total mortalities for each stage

for (i in 1:8) {

 Z[i]<- M[i]+Fs[i]

}



###Stochastic fecundity

f ~ dunif(16, 41)

# Survival and transition calculations --------------------------------------

# Male survival for stages 1–4
for (i in 1:4) {
  Sa_m[i] <- exp(-Z[i])
}

# Female survival for stages 1 and 5–8
Sa_f[1] <- exp(-Z[1])  # same shared stage for neaonates
for (i in 2:5) {
  Sa_f[i] <- exp(-Z[i + 3])  # maps to Z[5:8]
}

##Probability of growing to the next step for males y_im[i]

  for (i in 1:4) {


    y_im[i] <- ifelse(i == 1, 
                      Sa_m[i], 
                      ((Sa_m[i]^T_im[i]) - (Sa_m[i]^(T_im[i] - 1))) / 
                      (Sa_m[i]^T_im[i] - 1) )
  }

# Transition probabilities G for males

 for (i in 1:4) {
  G_male[i] <- Sa_m[i] * y_im[i]
  p_male[i] <- ifelse(i == 1, 0, Sa_m[i] * (1 - y_im[i]))  # no self-loop for neonate males
  Dp_m[i]   <- 1 - G_male[i] - p_male[i]
}

##Probability of growing to the next step for females y_if[i]

  for (i in 1:5) {

    y_if[i] <- ifelse(i == 1, 
                      Sa_f[i], 
                      ((Sa_f[i]^T_if[i]) - (Sa_f[i]^(T_if[i] - 1))) / 
                      (Sa_f[i]^T_if[i] - 1) )
  }

# Transition probabilities G for females

 for (i in 1:5) {
  G_female[i] <- Sa_f[i] * y_if[i]
  p_female[i] <- ifelse(i == 1, 0, Sa_f[i] * (1 - y_if[i]))  # no self-loop for neonate females
  Dp_f[i]   <- 1 - G_female[i] - p_female[i]
}

#Calculated fecundity only for reproductive stages:

#Males
fec_a[4] <- Sa_m[4] * f * N_t[4] / (N_t[7] + N_t[4]) 

#Females
fec_a[7] <- Sa_f[4] * f * N_t[7] / (N_t[7] + N_t[4]) 


##BUILD Matrix A_m ----------------------------


# Note: p_m/f entries are shifted one step to the right because Jags does not allow to remove the first entry of p (neonates) with p<- p[-1]
 
A_t[1, 4] <- fec_a[4]
A_t[1, 7] <- G_female[4] * fec_a[7]
A_t[2, 1] <- p * G_male[1]
A_t[2, 2] <- p_male[2]
A_t[3, 2] <- G_male[2]
A_t[3, 3] <- p_male[3]
A_t[4, 3] <- G_male[3]
A_t[4, 4] <- p_male[4]
A_t[5, 1] <- (1 - p) * G_male[1]
A_t[5, 5] <- p_female[2]
A_t[6, 5] <- G_female[1]
A_t[6, 6] <- p_female[3]
A_t[7, 6] <- G_female[2]
A_t[7, 8] <- G_female[4]
A_t[8, 7] <- G_female[5]


}  # end model
")
sink()


## Parameter definition ---------------------------------------------------

#Parameters

# Sex ratio (male/female)
p <- 0.5  

# Simulated input data
N_t <- c(10000, 1000, 500, 400, 1000, 500, 200, 200)


# Time in each life stage (same order as your T_im)

T_i <-c(1, 1,amat_m-1,amax-amat_m, 1, amat_f-1, 1, 1)

T_im<- T_i[c(1:4)]
T_if<- T_i[c(1,5:8)]

#Create prior parameter list

win.data <- list(
  N_t = N_t,
  T_im = T_im,
  T_if = T_if,
  p = p,
  M_lit=M_lit
)

#Define output to save un model object
params <- c("M", "Sa_m", "y_im", "G_male", "p_male", 
            "Sa_f", "y_if", "G_female", "p_female", 
            "A_t", "Fs", "Z", "f" )

#Run model
out_matrix <- jags(data = win.data,
                   inits = NULL,
                   parameters.to.save = params,
                   model.file = "MatrixModel_Base.txt",
                   n.chains = 2,
                   n.iter = 10000,
                   n.burnin = 1000,
                   n.thin = 1,
                   DIC=F)


#MCMC vis to quickly check traceplots and densities
MCMCtrace(out_matrix, 
          params = c("M"), 
          pdf = FALSE)


#Store Sample Distributions of mortalities

out_matrix$BUGSoutput$sims.list$G_female

# Total mortality
Z_all <- out_matrix$BUGSoutput$sims.list$Z  
#Natural mortality
M_all <- out_matrix$BUGSoutput$sims.list$M  
#Fishing mortality
Fs_all <- out_matrix$BUGSoutput$sims.list$Fs  


# Calculate the mean mortality across stages for each draw
mean_Z <- rowMeans(Z_all)
mean_M <- rowMeans(M_all)
mean_Fs <- rowMeans(Fs_all)


# Overall fecundity
f_all <- out_matrix$BUGSoutput$sims.list$f  

# Vectorize fecundity. single value for each draw
mean_f <- rowMeans(f_all)

#Extract matrix entries

# Extract all iterations of monitored parameters
A_raw <- out_matrix$BUGSoutput$sims.matrix

# Get the number of MCMC samples
n_iter <- nrow(A_raw)

# Initialize a list to store matrices per iteration
A_t_list <- vector("list", n_iter)

# Loop over iterations and fill in the matrix
for (i in 1:n_iter) {
  A <- matrix(0, nrow = 8, ncol = 8)
  
  A[1, 4] <- A_raw[i, "A_t[1,4]"]
  A[1, 7] <- A_raw[i, "A_t[1,7]"]
  A[2, 1] <- A_raw[i, "A_t[2,1]"]
  A[2, 2] <- A_raw[i, "A_t[2,2]"]
  A[3, 2] <- A_raw[i, "A_t[3,2]"]
  A[3, 3] <- A_raw[i, "A_t[3,3]"]
  A[4, 3] <- A_raw[i, "A_t[4,3]"]
  A[4, 4] <- A_raw[i, "A_t[4,4]"]
  A[5, 1] <- A_raw[i, "A_t[5,1]"]
  A[5, 5] <- A_raw[i, "A_t[5,5]"]
  A[6, 5] <- A_raw[i, "A_t[6,5]"]
  A[6, 6] <- A_raw[i, "A_t[6,6]"]
  A[7, 6] <- A_raw[i, "A_t[7,6]"]
  A[7, 8] <- A_raw[i, "A_t[7,8]"]
  A[8, 7] <- A_raw[i, "A_t[8,7]"]
  A_t_list[[i]] <- A
}


#Mean Matrix
# Compute mean matrix from list of matrices
A_t_mean <- Reduce("+", A_t_list) / n_iter

#Lambdas

# Calculate population growth rate (lambda)

lambdas <- sapply(A_t_list, function(A) {
  eigvals <- eigen(A)$values
  Re(eigvals[which.max(Re(eigvals))])  # Get the dominant real part
})

x<-sapply(A_t_list, function(A) {
  eigvals <- eigen(A)$values
})


#Mean lambda
# Compute mean matrix from list of matrices
mean_lambda <- mean(lambdas)
sd(lambdas)

#n lambdas, just to chek its one per itteration

length(lambdas)

#Lambdas with a value larger than 1

length(lambdas[lambdas>1])

#Proportion


length(lambdas[lambdas>1])/length(lambdas)

mean(lambdas)
sd(lambdas)

##Append mortalities and lambda to matrix list

A_t_list$M<-mean_M
A_t_list$Fs<-mean_Fs
A_t_list$Z<-mean_Z
A_t_list$fecundity<-mean_f
A_t_list$lambda<-lambdas

# 1. Neonates
A_t_list$M_Neonates         <- M_all[, 1]
A_t_list$Fs_Neonates        <- Fs_all[, 1]
A_t_list$Z_Neonates         <- Z_all[, 1]

# 2. Juvenile Males
A_t_list$M_Juvenile_Males   <- M_all[, 2]
A_t_list$Fs_Juvenile_Males  <- Fs_all[, 2]
A_t_list$Z_Juvenile_Males   <- Z_all[, 2]

# 3. Subadult Males
A_t_list$M_Subadult_Males   <- M_all[, 3]
A_t_list$Fs_Subadult_Males  <- Fs_all[, 3]
A_t_list$Z_Subadult_Males   <- Z_all[, 3]

# 4. Adult Males
A_t_list$M_Adult_Males      <- M_all[, 4]
A_t_list$Fs_Adult_Males     <- Fs_all[, 4]
A_t_list$Z_Adult_Males      <- Z_all[, 4]

# 5. Juvenile Females
A_t_list$M_Juvenile_Females <- M_all[, 5]
A_t_list$Fs_Juvenile_Females<- Fs_all[, 5]
A_t_list$Z_Juvenile_Females <- Z_all[, 5]

# 6. Subadult Females
A_t_list$M_Subadult_Females <- M_all[, 6]
A_t_list$Fs_Subadult_Females<- Fs_all[, 6]
A_t_list$Z_Subadult_Females <- Z_all[, 6]

# 7. Adult Females
A_t_list$M_Adult_Females    <- M_all[, 7]
A_t_list$Fs_Adult_Females   <- Fs_all[, 7]
A_t_list$Z_Adult_Females    <- Z_all[, 7]

# 8. Resting Females
A_t_list$M_Resting_Females  <- M_all[, 8]
A_t_list$Fs_Resting_Females <- Fs_all[, 8]
A_t_list$Z_Resting_Females  <- Z_all[, 8]



## Iterative population projection -----------------------------------------

# PARAMETERS

n_iter <- length(A_t_list$lambda) #Extract number of model itterations
n_years_b <- 36  #Years to project population
n_stages <- length(N_t)  #Number of lifestages
N_t <- c(10000, 1000, 500, 400, 1000, 500, 200, 200) #Initial population


# Preallocate empty list to store results
popb_summary_list <- vector("list", n_iter)


# LOOP: project total population over time for each iteration
for (i in 1:n_iter) {
  N_temp <- N_t  # use a temporary variable
  popb_iter <- numeric(n_years_b + 1)
  popb_iter[1] <- sum(N_temp)
  popb_sub <- numeric(n_years_b + 1)
  popb_sub[1] <- sum(N_temp[c(3,4,6,7,8)])
  
  for (t in 1:n_years_b) {
    N_temp <- A_t_list[[i]] %*% N_temp
    popb_iter[t + 1] <- sum(N_temp)
    popb_sub[t + 1] <- sum(N_temp[c(3,4,6,7,8)])
  }
  
  
  # STORE RESULTS 
  popb_summary_list[[i]] <- data.frame(
    Iteration = i,
    Year = 0:n_years_b,
    
    # Total populations
    Population = popb_iter,
    SubPop = popb_sub,
    
    Natural = A_t_list$M[i],
    Fishing= A_t_list$Fs[i],
    Total = A_t_list$Z[i],
    
    # Natural mortality by stage
    Natural_Neonates         = A_t_list$M_Neonates[i],
    Natural_Juvenile_Males   = A_t_list$M_Juvenile_Males[i],
    Natural_Subadult_Males   = A_t_list$M_Subadult_Males[i],
    Natural_Adult_Males      = A_t_list$M_Adult_Males[i],
    Natural_Juvenile_Females = A_t_list$M_Juvenile_Females[i],
    Natural_Subadult_Females = A_t_list$M_Subadult_Females[i],
    Natural_Adult_Females    = A_t_list$M_Adult_Females[i],
    Natural_Resting_Females  = A_t_list$M_Resting_Females[i],
    
    # Fishing mortality by stage
    Fishing_Neonates         = A_t_list$Fs_Neonates[i],
    Fishing_Juvenile_Males   = A_t_list$Fs_Juvenile_Males[i],
    Fishing_Subadult_Males   = A_t_list$Fs_Subadult_Males[i],
    Fishing_Adult_Males      = A_t_list$Fs_Adult_Males[i],
    Fishing_Juvenile_Females = A_t_list$Fs_Juvenile_Females[i],
    Fishing_Subadult_Females = A_t_list$Fs_Subadult_Females[i],
    Fishing_Adult_Females    = A_t_list$Fs_Adult_Females[i],
    Fishing_Resting_Females  = A_t_list$Fs_Resting_Females[i],
    
    # Total mortality by stage
    Total_Neonates           = A_t_list$Z_Neonates[i],
    Total_Juvenile_Males     = A_t_list$Z_Juvenile_Males[i],
    Total_Subadult_Males     = A_t_list$Z_Subadult_Males[i],
    Total_Adult_Males        = A_t_list$Z_Adult_Males[i],
    Total_Juvenile_Females   = A_t_list$Z_Juvenile_Females[i],
    Total_Subadult_Females   = A_t_list$Z_Subadult_Females[i],
    Total_Adult_Females      = A_t_list$Z_Adult_Females[i],
    Total_Resting_Females    = A_t_list$Z_Resting_Females[i],
    
    # Vital rates
    Fecundity = A_t_list$fecundity[i],
    lambda    = A_t_list$lambda[i]
  )
}



# Combine all iterations into one data frame
popb_proj_df <- do.call(rbind, popb_summary_list)

popb_summary <- popb_proj_df %>%
  group_by(Year) %>%
  summarise(
    mean = mean(Population),
    lower = quantile(Population, 0.025),
    upper = quantile(Population, 0.975)
  )


#Spaghetti plots

#Ease visualization trimming extreme values
popb_proj_dftrim<-popb_proj_df[popb_proj_df$Population<10000 & popb_proj_df$Population> 0,]

#THIS TRIM PICKs UP lAMBDA CRITICAL VALUES ONLY!!
popb_proj_dftrim2<-popb_proj_df[popb_proj_df$lambda<1.025 & popb_proj_df$lambda> 0.975,]

#Untrimmed
pop_proj_full<-ggplot(popb_proj_df, aes(x = Year, y = Population, group = Iteration)) +
  geom_line(color = "grey") +
  stat_summary(fun = mean, geom = "line", color = "blue", aes(group = 1), linetype = "dashed")+ theme_classic()+
  labs(title = "Population Through Years Control Model",
       x = "Years",
       y = "Population") +
  theme_minimal()

print(pop_proj_full)

#Trimmed
ggplot(popb_proj_dftrim, aes(x = Year, y = Population, group = factor(Iteration))) +
  geom_line(color = "grey") +
  stat_summary(fun = mean, geom = "line", color = "blue", linetype = "dashed", aes(group = 1))+ theme_classic()+
  labs(title = "Population Through Years Control Model (trimed)",
       x = "Years",
       y = "Population") +
  theme_minimal()


## Critical values ---------------------------------------------------------


pop_year_critical<-ggplot(popb_proj_dftrim2, aes(x = Year, y = Population, group = factor(Iteration))) +
  geom_line(color = "grey") +
  stat_summary(fun = mean, geom = "line", color = "blue", linetype = "dashed", aes(group = 1))+ theme_classic()+
  labs(title = "Population Through Years ~ Mortality",
       x = "Years",
       y = "Population") +
  theme_minimal()

print(pop_year_critical)

ggsave("Theoretical_Critical_population.jpg", plot = pop_year_critical,
       width = 12, height = 10, dpi = 300)




# Comparing virtual population to data slopes -----------------------------

Hammerheads_Cocos<-read.csv("Lewini_Cocos.csv", header=T)[-1]
Hammerheads_Galapagos<-read.csv("Lewini_Galapagos.csv", header=T)[-1]

Hammerheads_Cocos$Date<-as.Date(Hammerheads_Cocos$Date, format = "%Y-%m-%d")
Hammerheads_Galapagos$Date<-as.Date(Hammerheads_Galapagos$Date, format = "%Y-%m-%d")

Hammerheads_ETP<-rbind(Hammerheads_Cocos,Hammerheads_Galapagos)

#Group by site

Hammerheads_ETP_Site<-Hammerheads_ETP%>%
  group_by(Region, Date, Year, Month, Day, Doy, Site)%>%
  summarise(MeanCount=mean(Count),
            N_divers = n()
  )

#Group by site

Hammerheads_ETP_grouped<-Hammerheads_ETP_Site%>%
  group_by(Region, Date, Year, Month, Day, Doy)%>%
  summarise(MeanCount=mean(MeanCount),
            Mean_divers = mean(N_divers),
            N_sites = n()
  )

#Group by year
Hammerheads_ETP_year <-Hammerheads_ETP_grouped%>%
  group_by(Year,Region)%>%
  summarise(Mean_divers = mean(Mean_divers),
            Mean_sites = mean(N_sites),
            Effort = Mean_divers*Mean_sites,
            MeanCount=mean(MeanCount),
            PonderedCount=MeanCount/Effort
  )



Hammerheads_ETP_year_merged<-Hammerheads_ETP_year%>%group_by(Year)%>%
  summarise(MeanCount=mean(MeanCount),
            PonderedCount=mean(PonderedCount))

#Create a simple lm to have a population trend to compare to


Hammerheads_observed_b<-Hammerheads_ETP_year_merged[c(1,2)]%>%
  mutate(Model="Observed")%>%
  mutate(popb_change = ((MeanCount) - first(MeanCount)) / first(MeanCount)) %>%
  ungroup()

hammer_lm_Obs_Change<-lm(data=Hammerheads_observed_b, popb_change~Year)
summary(hammer_lm_Obs_Change)

#visualize
plot(data = Hammerheads_observed_b, popb_change ~ Year)
abline(hammer_lm_Obs_Change, col = "red", lwd = 2)



#Trend from all simulation iterations

VP_pop_b_all<-popb_proj_df%>%
  dplyr::select(Iteration,Year,SubPop)%>%
  dplyr::rename(MeanCount=SubPop)%>%
  mutate(Model="Virtual")%>%
  mutate(popb_change = ((MeanCount) - first(MeanCount)) / first(MeanCount)) %>%
  ungroup()%>%
  filter(Year >= 5)


## Compare Itterations of slopes -------------------------------------------

##Observed model slope

obs_fit   <- lm(popb_change ~ Year, data = Hammerheads_observed_b)
obs_slope <- coef(obs_fit)["Year"]
obs_se    <- summary(obs_fit)$coefficients["Year", "Std. Error"]

## Fit a slope per iteration 

vp_slopes <- VP_pop_b_all %>%
  group_by(Iteration) %>%
  do(tidy(lm(popb_change ~ Year, data = .))) %>%
  ungroup() %>%
  filter(term == "Year") %>%
  transmute(Iteration,
            vp_slope = estimate,
            vp_se    = std.error)

## Compare each iteration’s slope to the observed slope

comparisons <- vp_slopes %>%
  mutate(obs_slope = obs_slope,
         diff      = vp_slope - obs_slope,
         # simple z using both SEs (independent-slope approximation)
         z         = (vp_slope - obs_slope) / sqrt(vp_se^2 + obs_se^2),
         p_two_sided = 2 * pnorm(-abs(z)))

#Keep models with non signifficant difference
simmilar_slopes<-comparisons[comparisons$p_two_sided>0.05,]


##Evaluate percentage of iterations above and below lambda. ---------------


#proportion of iterations within lambda critical


common_iterations <- intersect(
  unique(simmilar_slopes$Iteration),
  unique(popb_proj_df$Iteration)
)

mean(popb_proj_df$lambda)
mean(popb_proj_dftrim$lambda)
mean(popb_proj_dftrim2$lambda)


posterior_itterations<-length(unique(popb_proj_df$Iteration))

modelled_itterations<-length(unique(simmilar_slopes$Iteration))

lambda1_itterations<-length(unique(popb_proj_dftrim2$Iteration))

lambdaPlus1_itterations<-length(unique(popb_proj_df$Iteration[popb_proj_df$lambda>1]))

lambdaMinus1_itterations<-length(unique(popb_proj_df$Iteration[popb_proj_df$lambda<1]))

#lambda1 / Posterior

lambda1_itterations  / posterior_itterations

#Lambda>1 / posterior
lambdaPlus1_itterations / posterior_itterations

#Lambda<1 / posterior
lambdaMinus1_itterations / posterior_itterations

#simmilar slope/ posterior
modelled_itterations / posterior_itterations


#lambda 1 /simmilar slope
lambda1_itterations / modelled_itterations




#Lambda >1 / lambda1
lambdaPlus1_itterations / lambda1_itterations 


##Overal simmulations with a simmilar slope  -----------------------------



sim_pop<-popb_proj_df%>%
  filter(Iteration %in% simmilar_slopes$Iteration)

#Check lambda distribution of simmilar VPs

sim_pop_lambdas<-sim_pop|>group_by(Iteration)|>
  select(c(-Year, -Population, -SubPop))|>
  unique()


mean(sim_pop_lambdas$lambda)
median(sim_pop_lambdas$lambda)
sd(sim_pop_lambdas$lambda)

pop_year_sim<-ggplot(sim_pop, aes(x = Year, y = Population, group = factor(Iteration))) +
  geom_line(color = "grey") +
  stat_summary(fun = mean, geom = "line", color = "blue", linetype = "dashed", aes(group = 1))+ theme_classic()+
  labs(title = " ",
       x = "Years",
       y = "Population") +
  theme_minimal(base_size = 30)

print(pop_year_sim)

ggsave("Theoretical_Observed_population.jpg", plot = pop_year_sim,
       width = 12, height = 10, dpi = 300)

##Lambda critical simmulations with a simmilar slope  -----------------------------


sim_pop_critical<-popb_proj_dftrim2%>%
  filter(Iteration %in% simmilar_slopes$Iteration)

#Check lambda distribution of simmilar VPs

sim_pop_critical_lambdas<-sim_pop_critical|>group_by(Iteration)|>
  select(c(-Year, -Population, -SubPop))|>
  unique()


mean(sim_pop_critical_lambdas$lambda)
median(sim_pop_critical_lambdas$lambda)
sd(sim_pop_critical_lambdas$lambda)

pop_year_sim_critical<-ggplot(sim_pop_critical, aes(x = Year, y = Population, group = factor(Iteration))) +
  geom_line(color = "grey") +
  stat_summary(fun = mean, geom = "line", color = "blue", linetype = "dashed", aes(group = 1))+ theme_classic()+
  labs(title = " ",
       x = "Years",
       y = "Population") +
  theme_minimal(base_size = 30)

print(pop_year_sim_critical)

ggsave("Theoretical_Observed_Critical_population.jpg", plot = pop_year_sim_critical,
       width = 12, height = 10, dpi = 300)


##Different Projection Visualization --------------------------------------

#Untrimmed
print(pop_proj_full)

#Lambda critical
print(pop_year_critical)

#Adjusted to observed
print(pop_year_sim)

#Adjusted to observed + critical
print(pop_year_sim_critical)


# DeterminePopulationParameters -------------------------------------------



#Long format of filtered estimates
sim_pop_long <- sim_pop%>%
  dplyr::select(-Population, -SubPop,-Year) %>%
  pivot_longer(
    cols = -c(Iteration),  # all columns to pivot
    names_to = "Variable",   # name for the new column holding original column names
    values_to = "Estimate"   # name for the new column holding the values
  )



#Group and get minimum and maximum parameter value per estimate

parameter_summary<-sim_pop_long%>%group_by(Variable)%>%
  summarize(mean=mean(Estimate),
            minimum=min(Estimate),
            maximum=max(Estimate),
            SD=sd(Estimate))%>%
  ungroup()

print(parameter_summary)

write.csv(parameter_summary,"ParameterEstimates.csv")


#With additional columns 

parameter_summary_extra<-parameter_summary%>%mutate(
  # Extract the base parameter 
  Parameter = str_extract(Variable, "^[^_]+"),
  
  # Extract everything after the first underscore (stage + sex)
  Detail = str_remove(Variable, "^[^_]+_"),
  
  # Extract Stage (first word after the parameter)
  Stage = case_when(
    str_detect(Detail, "Neonates") ~ "Neonates",
    str_detect(Detail, "Juvenile") ~ "Juvenile",
    str_detect(Detail, "Subadult") ~ "Subadult",
    str_detect(Detail, "Adult")    ~ "Adult",
    str_detect(Detail, "Resting")  ~ "Resting",
    TRUE ~ NA_character_
  ),
  
  # Extract Sex
  Sex = case_when(
    str_detect(Detail, "Males")   ~ "Male",
    str_detect(Detail, "Females") ~ "Female",
    TRUE ~ NA_character_
  )
) %>%
  select(-Detail)%>%   
  # Remove intermediate column
  mutate(Stage = recode(Stage,
                        "Juvenile" = "Juveniles",
                        "Subadult" = "Subadults",
                        "Adult" = "Adults"))

print(parameter_summary_extra)

write.csv(parameter_summary_extra,"ParameterEstimates_separated.csv")


## HoneycombFecundities ----------------------------------------------------


lambdaFecundity<-mixed_df|>select(Iteration,lambda,Fecundity,Model,Estimate)|>
  unique()


Lamb_F_Honey<-ggplot(lambdaFecundity[lambdaFecundity$Model=="Full",], aes(x = Fecundity, y = lambda)) +
  stat_binhex(bins = 24) +
  scale_fill_viridis_c(option="mako",
                       limits = c(0, 3000)) +
  theme_minimal(base_size = 30) +
  labs(
    x = "Fecundity",
    y = "Population growth rate (λ)",
    fill = "Iteration density",
    title = " "
  )+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.key.width =unit(2, "cm")
        )
print(Lamb_F_Honey)

ggsave("Honeycomb_LAmbdaVsFecundity.jpg", plot = Lamb_F_Honey,
       width = 12, height = 10, dpi = 300)


# Elasticity --------------------------------------------------------------

#Extract A_T Matrices of interest

class(A_t_list)
length(A_t_list)

str(sim_pop)
unique(sim_pop$Iteration)

# Get the iteration indices to keep
iter_keep <- unique(sim_pop$Iteration)

# Extract only those iterations from the list

A_t_extract <- A_t_list[iter_keep]

A_t_extract_mean <-Reduce("+", A_t_extract) / length(A_t_extract)


# Eigen decomposition
eigen_analysis <- eigen(A_t_extract_mean)
lambda <- Re(eigen_analysis$values[1])  # Dominant eigenvalue
w <- Re(eigen_analysis$vectors[, 1])   # Right eigenvector (stable stage distribution)
v <- Re(eigen(t(A_t_extract_mean))$vectors[, 1])    # Left eigenvector (reproductive values)

# Normalize eigenvectors
w <- w / sum(w)  # Normalize stable stage distribution
v <- v / v[1]    # Normalize reproductive values by first element

# Scalar product of w and v
scalar_product <- sum(w * v)

# Elasticity matrix calculation
elasticity_matrix <- matrix(0, nrow = nrow(A_t_extract_mean), ncol = ncol(A_t_extract_mean))  # Initialize

for (i in 1:nrow(A_t_extract_mean)) {
  for (j in 1:ncol(A_t_extract_mean)) {
    a_ij <- A_t_extract_mean[i, j]  # Element of A_t
    elasticity_matrix[i, j] <- (a_ij / lambda) * (v[i] * w[j] / scalar_product)
  }
}

# 5. Sum of elasticities (check if it equals 1)
sum_elasticities <- sum(elasticity_matrix)

#Examine elasticity matrix
print(round(elasticity_matrix, 4))


# Compute elasticity matrices for each iteration
##please note that "elas" is a function from popdemo package. Will not work if package is not loaded. 

elasticity_list <- lapply(A_t_extract, elas)

# Convert list of matrices to data frame
elasticity_df <- map_dfr(
  seq_along(elasticity_list),
  function(i) {
    mat <- elasticity_list[[i]]
    as_tibble(mat) |>
      mutate(row = row_number()) |>
      pivot_longer(
        cols = -row,
        names_to = "col",
        values_to = "elasticity"
      ) |>
      mutate(
        col = as.integer(gsub("V", "", col)),
        iteration = i
      )
  }
)


##Elasticity matrix plot --------------------------------------------------

elasticity_matrix_plot<-elasticity_df |>
  mutate(
    row = factor(row, levels = sort(unique(row), decreasing = FALSE)),
    col = factor(col, levels = sort(unique(col)))
  )|>
  filter(elasticity != 0)|>
  ggplot(aes(x = elasticity)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_grid(row ~ col) +
  theme_minimal() +
  labs(
    title = "Elasticity distributions per matrix element",
    x = "Elasticity",
    y = "Density"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.5, "lines")
  )

plot(elasticity_matrix_plot)

ggsave("elasticity_matrix_plot.jpg", plot = elasticity_matrix_plot,
       width = 12, height = 10, dpi = 300)



##Summirized elasticity ---------------------------------------------------


# Tag important transitions and extract relevant data.
elasticity_df_nonzero <- elasticity_df %>%
  filter(elasticity != 0) %>%
  mutate(cell = paste0("(", row, ",", col, ")"),
         sex = case_when(
           col == 1 ~ "Neonates",
           col %in% 2:4 ~ "Males",
           col %in% 5:8 ~ "Females",
           TRUE ~ NA_character_
         ),
         Stage = case_when(
           col == 1 ~ "Neonates",
           col == 2 ~ "Juveniles",
           col == 3 ~ "Subadults",
           col == 4 ~ "Adults",
           col == 5 ~ "Juveniles",
           col == 6 ~ "Subadults",
           col == 7 ~ "Adults",
           col == 8 ~ "Adults",
           TRUE ~ NA_character_
         )
  )%>%
  group_by(iteration,Stage,sex)%>%
  summarise(elasticity=sum(elasticity))%>%
  mutate(
    Stage = factor(Stage, levels = c("Neonates", "Juveniles", "Subadults", "Adults", "Resting"))
  )

#average by stage
Av_stage_elas<-elasticity_df_nonzero|>group_by(Stage,sex)%>%
  summarise(elasticity=mean(elasticity))
print(Av_stage_elas)

# Elasticities by transition type Check 
elasticity_transition <- elasticity_df %>%
  filter(elasticity != 0) %>%
  mutate(cell = paste0("(", row, ",", col, ")"),
         Transition_type = case_when(
           row == col ~ "Survival",
           
           # Growth transitions
           (col == 1 & row == 2) |
             (col == 2 & row == 3) | 
             (col == 3 & row == 4) |
             (col == 4 & row == 5) |
             (col == 5 & row == 6) | 
             (col == 6 & row == 7) | 
             (col == 7 & row == 8) |
             (col == 8 & row == 7) |
             (col == 1 & row == 5) ~ "Growth",
           
           # Reproduction transitions
           (row == 1 & col == 4) |
             (row == 1 & col == 7) ~ "Reproduction",
           TRUE ~ "Other"),
         sex = case_when(
           col == 1 ~ "Neonates",
           col %in% 2:4 ~ "Males",
           col %in% 5:8 ~ "Females",
           TRUE ~ NA_character_
         )
  )%>%
  group_by(iteration,Transition_type,sex)%>%
  summarise(elasticity=sum(elasticity))



#average by Transition Type
Av_Trans_elas<-elasticity_transition|>group_by(Transition_type,sex)%>%
  summarise(elasticity=mean(elasticity))
print(Av_Trans_elas)


##Elasticity plots --------------------------------------------------------



# Density distributions By Stage


# Plot
Density_stage_elas <- ggplot(elasticity_df_nonzero, aes(x = elasticity)) +
  geom_density(aes(fill = sex), alpha = 0.6) +
  scale_fill_viridis_d(option = "mako",begin=0.5) +   # mako palette for sex categories
  theme_minimal(base_size = 40) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "",
    x = "Elasticity value",
    y = "Density",
    fill = "Sex"
  ) +
  facet_wrap(~ Stage, scales = "free_y", nrow = 1)

print(Density_stage_elas)

# Density distribution by Transition
Density_Trans_elas<-elasticity_transition|>
  filter(sex != "Neonates") |> 
  ggplot(aes(x = elasticity)) +
  geom_density(aes(fill = sex), alpha = 0.7) +
  scale_fill_viridis_d(option = "mako",begin=0.4, end =0.7) +   # mako palette for sex categories
  theme_minimal(base_size = 40) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "",
    x = "Elasticity value",
    y = "Density",
    fill = "Sex"
  ) +
  facet_wrap(~ Transition_type, scales = "free_y", nrow = 1)

print(Density_Trans_elas)


#Barplot by stage

Bar_stage_elas <- Av_stage_elas %>%
  filter(sex != "Neonates") %>%
  ggplot(aes(x = sex, y = elasticity, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill", color = "white") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_viridis_d(option = "mako",end = 0.3) +   # mako palette for stages
  labs(
    title = "",
    x = NULL,
    y = "Proportion of total elasticity",
    fill = "Stage"
  ) +
  theme_minimal(base_size = 40) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 0, r = 200, b = 0, l = 0),
    panel.grid = element_blank()
  )

print(Bar_stage_elas)


#Combined By Stage

Combined_plot_stage <- Density_stage_elas + Bar_stage_elas +
  plot_layout(widths = c(2, 1)) 

print(Combined_plot_stage)


ggsave("ElasticityPlotStages.jpg", plot = Combined_plot_stage,
       width = 25, height = 12, dpi = 300)



#Barplot by Transition

Bar_Trans_elas<-Av_Trans_elas %>%
  filter(sex != "Neonates") %>%
  ggplot(aes(x = sex, y = elasticity, fill = Transition_type)) +
  geom_bar(stat = "identity", position = "fill", color = "white") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_viridis_d(option = "viridis",begin = 0.5, direction = -1) +   # mako palette for stages
  labs(
    title = "",
    x = NULL,
    y = "Proportion of total elasticity",
    fill = "Parameter"
  ) +
  theme_minimal(base_size = 40) +
  theme(
    plot.title = element_text(hjust = 0.6, size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 0, r = 200, b = 0, l = 0),
    panel.grid = element_blank()
  )


print(Bar_Trans_elas)


#Combined By Transition

Combined_plot_Transition <- Density_Trans_elas + Bar_Trans_elas +
  plot_layout(widths = c(2, 1)) 

print(Combined_plot_Transition)

ggsave("ElasticityPlot.jpg", plot = Combined_plot_Transition,
       width = 25, height = 12, dpi = 300)


# Scenario models ---------------------------------------------------------

#Extract ETP parameters + uncertainity

# NEONATES 
sd_Ms_Neo <- parameter_summary$SD[parameter_summary$Variable == "Natural_Neonates"]
sd_Fs_Neo <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Neonates"]

Mean_Ms_Neo <- parameter_summary$mean[parameter_summary$Variable == "Natural_Neonates"]
Mean_Fs_Neo <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Neonates"]

# JUVENILES - Females
sd_Ms_Juv_Females <- parameter_summary$SD[parameter_summary$Variable == "Natural_Juvenile_Females"]

sd_Fs_Juv_Females <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Juvenile_Females"]

Mean_Ms_Juv_Females <- parameter_summary$mean[parameter_summary$Variable == "Natural_Juvenile_Females"]
Mean_Fs_Juv_Females <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Juvenile_Females"]

# JUVENILES - Males
sd_Ms_Juv_Males <- parameter_summary$SD[parameter_summary$Variable == "Natural_Juvenile_Males"]
sd_Fs_Juv_Males <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Juvenile_Males"]

Mean_Ms_Juv_Males <- parameter_summary$mean[parameter_summary$Variable == "Natural_Juvenile_Males"]
Mean_Fs_Juv_Males <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Juvenile_Males"]

# SUBADULTS - Females

sd_Ms_Subadult_Females <- parameter_summary$SD[parameter_summary$Variable == "Natural_Subadult_Females"]
sd_Fs_Subadult_Females <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Subadult_Females"]

Mean_Ms_Subadult_Females <- parameter_summary$mean[parameter_summary$Variable == "Natural_Subadult_Females"]
Mean_Fs_Subadult_Females <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Subadult_Females"]

# SUBADULTS - Males
sd_Ms_Subadult_Males <- parameter_summary$SD[parameter_summary$Variable == "Natural_Subadult_Males"]
sd_Fs_Subadult_Males <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Subadult_Males"]

Mean_Ms_Subadult_Males <- parameter_summary$mean[parameter_summary$Variable == "Natural_Subadult_Males"]
Mean_Fs_Subadult_Males <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Subadult_Males"]

# ADULTS - Females
sd_Ms_Adult_Females <- parameter_summary$SD[parameter_summary$Variable == "Natural_Adult_Females"]
sd_Fs_Adult_Females  <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Adult_Females"]

Mean_Ms_Adult_Females <- parameter_summary$mean[parameter_summary$Variable == "Natural_Adult_Females"]
Mean_Fs_Adult_Females <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Adult_Females"]

# ADULTS - Males
sd_Ms_Adult_Males  <- parameter_summary$SD[parameter_summary$Variable == "Natural_Adult_Males"]
sd_Fs_Adult_Males  <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Adult_Males"]

Mean_Ms_Adult_Males <- parameter_summary$mean[parameter_summary$Variable == "Natural_Adult_Males"]
Mean_Fs_Adult_Males <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Adult_Males"]

# RESTING_FEMALES (only females)
sd_Ms_Resting_Females  <- parameter_summary$SD[parameter_summary$Variable == "Natural_Resting_Females"]
sd_Fs_Resting_Females  <- parameter_summary$SD[parameter_summary$Variable == "Fishing_Resting_Females"]

Mean_Ms_Resting_Females <- parameter_summary$mean[parameter_summary$Variable == "Natural_Resting_Females"]
Mean_Fs_Resting_Females <- parameter_summary$mean[parameter_summary$Variable == "Fishing_Resting_Females"]


# No fishing mortality ----------------------------------------------------


sink("MatrixModel_NoFs.txt")
cat("
model {

###PRIORS-----------------------------------
  
  
#Stochastic random mortalities for all stages



#set sd/tau for stages and mortalities

 #Neonates
  
tau_Ms_Neo<- 1/(sd_Ms_Neo^2)

 #Juveniles
 
tau_Ms_Juv_Females <- 1 / (sd_Ms_Juv_Females^2)


tau_Ms_Juv_Males <- 1 / (sd_Ms_Juv_Males^2)


#Subadults
 
tau_Ms_Subadult_Females <- 1 / (sd_Ms_Subadult_Females^2)


tau_Ms_Subadult_Males <- 1 / (sd_Ms_Subadult_Males^2)

 
#Adults

tau_Ms_Adult_Females <- 1 / (sd_Ms_Adult_Females^2)


tau_Ms_Adult_Males <- 1 / (sd_Ms_Adult_Males^2)


#Resting females

tau_Ms_Resting_Females <- 1 / (sd_Ms_Resting_Females^2)


#Stage specific mortality priors using parameterization.

# Neonates
M[1]  ~ dnorm(Mean_Ms_Neo,  tau_Ms_Neo)


# Juvenile_Males
M[2]  ~ dnorm(Mean_Ms_Juv_Males,  tau_Ms_Juv_Males)


# Subadult_Males
M[3]  ~ dnorm(Mean_Ms_Subadult_Males, tau_Ms_Subadult_Males)


# Adult_Males
M[4]  ~ dnorm(Mean_Ms_Adult_Males, tau_Ms_Adult_Males)


# Juvenile_Females
M[5]  ~ dnorm(Mean_Ms_Juv_Females,  tau_Ms_Juv_Females)


# Subadult_Females
M[6]  ~ dnorm(Mean_Ms_Subadult_Females, tau_Ms_Subadult_Females)


# Adult_Females
M[7]  ~ dnorm(Mean_Ms_Adult_Females, tau_Ms_Adult_Females)

# Resting_Females
M[8]  ~ dnorm(Mean_Ms_Resting_Females, tau_Ms_Resting_Females)



#Total mortalities for each stage

for (i in 1:8) {

 Z[i]<- M[i]

}



###Stochastic fecundity


#Still diffuse
f ~ dunif(16, 41)

# Survival and transition calculations --------------------------------------

# Male survival for stages 1–4
for (i in 1:4) {
  Sa_m[i] <- exp(-Z[i])
}

# Female survival for stages 1 and 5–8
Sa_f[1] <- exp(-Z[1])  # same shared stage for neaonates
for (i in 2:5) {
  Sa_f[i] <- exp(-Z[i + 3])  # maps to Z[5:8]
}

##Probability of growing to the next step for males y_im[i]

  for (i in 1:4) {


    y_im[i] <- ifelse(i == 1, 
                      Sa_m[i], 
                      ((Sa_m[i]^T_im[i]) - (Sa_m[i]^(T_im[i] - 1))) / 
                      (Sa_m[i]^T_im[i] - 1) )
  }

# Transition probabilities G for males

 for (i in 1:4) {
  G_male[i] <- Sa_m[i] * y_im[i]
  p_male[i] <- ifelse(i == 1, 0, Sa_m[i] * (1 - y_im[i]))  # no self-loop for neonate males
  Dp_m[i]   <- 1 - G_male[i] - p_male[i]
}

##Probability of growing to the next step for females y_if[i]

  for (i in 1:5) {

    y_if[i] <- ifelse(i == 1, 
                      Sa_f[i], 
                      ((Sa_f[i]^T_if[i]) - (Sa_f[i]^(T_if[i] - 1))) / 
                      (Sa_f[i]^T_if[i] - 1) )
  }

# Transition probabilities G for females

 for (i in 1:5) {
  G_female[i] <- Sa_f[i] * y_if[i]
  p_female[i] <- ifelse(i == 1, 0, Sa_f[i] * (1 - y_if[i]))  # no self-loop for neonate females
  Dp_f[i]   <- 1 - G_female[i] - p_female[i]
}

#Calculated fecundity only for reproductive stages:

#Males
fec_a[4] <- Sa_m[4] * f * N_t[4] / (N_t[7] + N_t[4]) 

#Females
fec_a[7] <- Sa_f[4] * f * N_t[7] / (N_t[7] + N_t[4]) 


##BUILD Matrix A_m ----------------------------


# Note: p_m/f entries are shifted one step to the right because Jags does not allow to remove the first entry of p (neonates) with p<- p[-1]
 
A_t[1, 4] <- fec_a[4]
A_t[1, 7] <- G_female[4] * fec_a[7]
A_t[2, 1] <- p * G_male[1]
A_t[2, 2] <- p_male[2]
A_t[3, 2] <- G_male[2]
A_t[3, 3] <- p_male[3]
A_t[4, 3] <- G_male[3]
A_t[4, 4] <- p_male[4]
A_t[5, 1] <- (1 - p) * G_male[1]
A_t[5, 5] <- p_female[2]
A_t[6, 5] <- G_female[1]
A_t[6, 6] <- p_female[3]
A_t[7, 6] <- G_female[2]
A_t[7, 8] <- G_female[4]
A_t[8, 7] <- G_female[5]


}  # end model
")
sink()


## Parameter definition ---------------------------------------------------

#Parameters

#Create prior parameter list

win.data <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  Mean_Ms_Neo = Mean_Ms_Neo,

  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,

  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,

  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females
  
)

#Define output to save un model object
params <- c("M", "Sa_m", "y_im", "G_male", "p_male", 
            "Sa_f", "y_if", "G_female", "p_female", 
            "A_t", "Z", "f" )

#Run model
out_matrix_noFs <- jags(data = win.data,
                   inits = NULL,
                   parameters.to.save = params,
                   model.file = "MatrixModel_NoFs.txt",
                   n.chains = 2,
                   n.iter = 10000,
                   n.burnin = 1000,
                   n.thin = 1,
                   DIC=F)


#MCMC vis to quickly check traveplots and densities
MCMCtrace(out_matrix_noFs, 
          params = c("M"), 
          pdf = FALSE)


#Store Sample Distributions of mortalities

out_matrix_noFs$BUGSoutput$sims.list$G_female

# Total mortality
Z_all <- out_matrix_noFs$BUGSoutput$sims.list$Z  
#Natural mortality
M_all <- out_matrix_noFs$BUGSoutput$sims.list$M  
#Fishing mortality
Fs_all <- out_matrix_noFs$BUGSoutput$sims.list$Fs  



# Calculate the mean mortality across stages for each draw
mean_Z <- rowMeans(Z_all)
mean_M <- rowMeans(M_all)
# mean_Fs <- rowMeans(Fs_all) #No fishing mortality


#Extract matrix entries

# Extract all iterations of monitored parameters
A_raw <- out_matrix_noFs$BUGSoutput$sims.matrix

# Get the number of MCMC samples
n_iter <- nrow(A_raw)

# Initialize a list to store matrices per iteration
A_t_list <- vector("list", n_iter)

# Loop over iterations and fill in the matrix
for (i in 1:n_iter) {
  A <- matrix(0, nrow = 8, ncol = 8)
  
  A[1, 4] <- A_raw[i, "A_t[1,4]"]
  A[1, 7] <- A_raw[i, "A_t[1,7]"]
  A[2, 1] <- A_raw[i, "A_t[2,1]"]
  A[2, 2] <- A_raw[i, "A_t[2,2]"]
  A[3, 2] <- A_raw[i, "A_t[3,2]"]
  A[3, 3] <- A_raw[i, "A_t[3,3]"]
  A[4, 3] <- A_raw[i, "A_t[4,3]"]
  A[4, 4] <- A_raw[i, "A_t[4,4]"]
  A[5, 1] <- A_raw[i, "A_t[5,1]"]
  A[5, 5] <- A_raw[i, "A_t[5,5]"]
  A[6, 5] <- A_raw[i, "A_t[6,5]"]
  A[6, 6] <- A_raw[i, "A_t[6,6]"]
  A[7, 6] <- A_raw[i, "A_t[7,6]"]
  A[7, 8] <- A_raw[i, "A_t[7,8]"]
  A[8, 7] <- A_raw[i, "A_t[8,7]"]
  A_t_list[[i]] <- A
}


#Mean Matrix
# Compute mean matrix from list of matrices
A_t_mean <- Reduce("+", A_t_list) / n_iter

#Lambdas

# Calculate population growth rate (lambda)

lambdas <- sapply(A_t_list, function(A) {
  eigvals <- eigen(A)$values
  Re(eigvals[which.max(Re(eigvals))])  # Get the dominant real part
})

x<-sapply(A_t_list, function(A) {
  eigvals <- eigen(A)$values
})

head(x)

#Mean lambda
# Compute mean matrix from list of matrices
mean_lambda <- mean(lambdas)
sd(lambdas)
#n lambdas, just to chek its one per itteration

length(lambdas)

#Lambdas with a value larger than 1

length(lambdas[lambdas>1])

#Proportion


length(lambdas[lambdas>1])/length(lambdas)

mean(lambdas)
sd(lambdas)
##Append mortalities and lambda to matrix list

A_t_list$M<-mean_M
A_t_list$Z<-mean_Z
A_t_list$fecundity<-mean_f
A_t_list$lambda<-lambdas

# 1. Neonates
A_t_list$M_Neonates         <- M_all[, 1]
A_t_list$Z_Neonates         <- Z_all[, 1]

# 2. Juvenile Males
A_t_list$M_Juvenile_Males   <- M_all[, 2]
A_t_list$Z_Juvenile_Males   <- Z_all[, 2]

# 3. Subadult Males
A_t_list$M_Subadult_Males   <- M_all[, 3]
A_t_list$Z_Subadult_Males   <- Z_all[, 3]

# 4. Adult Males
A_t_list$M_Adult_Males      <- M_all[, 4]
A_t_list$Z_Adult_Males      <- Z_all[, 4]

# 5. Juvenile Females
A_t_list$M_Juvenile_Females <- M_all[, 5]
A_t_list$Z_Juvenile_Females <- Z_all[, 5]

# 6. Subadult Females
A_t_list$M_Subadult_Females <- M_all[, 6]
A_t_list$Z_Subadult_Females <- Z_all[, 6]

# 7. Adult Females
A_t_list$M_Adult_Females    <- M_all[, 7]
A_t_list$Z_Adult_Females    <- Z_all[, 7]

# 8. Resting Females
A_t_list$M_Resting_Females  <- M_all[, 8]
A_t_list$Z_Resting_Females  <- Z_all[, 8]



## Iterative population projection -----------------------------------------

# PARAMETERS

n_iter <- length(A_t_list$lambda) #Extract number of model itterations
n_years <-  36  #Years to project population
n_stages <- length(N_t)  #Number of lifestages
N_t <- c(10000, 1000, 500, 400, 1000, 500, 200, 200) #Initial population


# Preallocate empty list to store results
pop_summary_list <- vector("list", n_iter)


# LOOP: project total population over time for each iteration
for (i in 1:n_iter) {
  N_temp <- N_t  # use a temporary variable!
  pop_iter <- numeric(n_years + 1)
  pop_iter[1] <- sum(N_temp)
  pop_sub <- numeric(n_years + 1)
  pop_sub[1] <- sum(N_temp[c(3,4,6,7,8)])
  
  for (t in 1:n_years) {
    N_temp <- A_t_list[[i]] %*% N_temp
    pop_iter[t + 1] <- sum(N_temp)
    pop_sub[t + 1] <- sum(N_temp[c(3,4,6,7,8)])
  }
  
  
  # STORE RESULTS 
  pop_summary_list[[i]] <- data.frame(
    Iteration = i,
    Year = 0:n_years,
    
    # Total populations
    Population = pop_iter,
    SubPop = pop_sub,
    Natural = A_t_list$M[i],
      Total = A_t_list$Z[i],
    
    # Natural mortality by stage
    Natural_Neonates         = A_t_list$M_Neonates[i],
    Natural_Juvenile_Males   = A_t_list$M_Juvenile_Males[i],
    Natural_Subadult_Males   = A_t_list$M_Subadult_Males[i],
    Natural_Adult_Males      = A_t_list$M_Adult_Males[i],
    Natural_Juvenile_Females = A_t_list$M_Juvenile_Females[i],
    Natural_Subadult_Females = A_t_list$M_Subadult_Females[i],
    Natural_Adult_Females    = A_t_list$M_Adult_Females[i],
    Natural_Resting_Females  = A_t_list$M_Resting_Females[i],
    
    
    # Total mortality by stage
    Total_Neonates           = A_t_list$Z_Neonates[i],
    Total_Juvenile_Males     = A_t_list$Z_Juvenile_Males[i],
    Total_Subadult_Males     = A_t_list$Z_Subadult_Males[i],
    Total_Adult_Males        = A_t_list$Z_Adult_Males[i],
    Total_Juvenile_Females   = A_t_list$Z_Juvenile_Females[i],
    Total_Subadult_Females   = A_t_list$Z_Subadult_Females[i],
    Total_Adult_Females      = A_t_list$Z_Adult_Females[i],
    Total_Resting_Females    = A_t_list$Z_Resting_Females[i],
    
    # Vital rates
    Fecundity = A_t_list$fecundity[i],
    lambda    = A_t_list$lambda[i]
  )
}

# Combine all iterations into one data frame
pop_proj_df <- do.call(rbind, pop_summary_list)

pop_summary <- pop_proj_df %>%
  group_by(Year) %>%
  summarise(
    mean = mean(Population),
    lower = quantile(Population, 0.025),
    upper = quantile(Population, 0.975)
  )


summary(pop_summary)



# Variations in Fishing mortality -----------------------------------------------

sink("MatrixModel_scenario.txt")
cat("
model {

###PRIORS-----------------------------------
  
  
#Stochastic random mortalities for all stages



#set sd/tau for stages and mortalities

 #Neonates
  
tau_Ms_Neo<- 1/(sd_Ms_Neo^2)
tau_Fs_Neo<- 1/(sd_Fs_Neo^2)

 #Juveniles
 
tau_Ms_Juv_Females <- 1 / (sd_Ms_Juv_Females^2)
tau_Fs_Juv_Females <- 1 / (sd_Fs_Juv_Females^2)

tau_Ms_Juv_Males <- 1 / (sd_Ms_Juv_Males^2)
tau_Fs_Juv_Males <- 1 / (sd_Fs_Juv_Males^2)

#Subadults
 
tau_Ms_Subadult_Females <- 1 / (sd_Ms_Subadult_Females^2)
tau_Fs_Subadult_Females <- 1 / (sd_Fs_Subadult_Females^2)

tau_Ms_Subadult_Males <- 1 / (sd_Ms_Subadult_Males^2)
tau_Fs_Subadult_Males <- 1 / (sd_Fs_Subadult_Males^2)
 
#Adults

tau_Ms_Adult_Females <- 1 / (sd_Ms_Adult_Females^2)
tau_Fs_Adult_Females <- 1 / (sd_Fs_Adult_Females^2)

tau_Ms_Adult_Males <- 1 / (sd_Ms_Adult_Males^2)
tau_Fs_Adult_Males <- 1 / (sd_Fs_Adult_Males^2)

#Resting females

tau_Ms_Resting_Females <- 1 / (sd_Ms_Resting_Females^2)
tau_Fs_Resting_Females <- 1 / (sd_Fs_Resting_Females^2)


#Stage specific mortality priors using parameterization.

# Neonates
M[1]  ~ dnorm(Mean_Ms_Neo,  tau_Ms_Neo)
Fs[1] ~ dnorm(Mean_Fs_Neo,tau_Fs_Neo)

# Juvenile_Males
M[2]  ~ dnorm(Mean_Ms_Juv_Males,  tau_Ms_Juv_Males)
Fs[2] ~ dnorm(Mean_Fs_Juv_Males,  tau_Fs_Juv_Males)

# Subadult_Males
M[3]  ~ dnorm(Mean_Ms_Subadult_Males, tau_Ms_Subadult_Males)
Fs[3] ~ dnorm(Mean_Fs_Subadult_Males, tau_Fs_Subadult_Males)

# Adult_Males
M[4]  ~ dnorm(Mean_Ms_Adult_Males, tau_Ms_Adult_Males)
Fs[4] ~ dnorm(Mean_Fs_Adult_Males, tau_Fs_Adult_Males)

# Juvenile_Females
M[5]  ~ dnorm(Mean_Ms_Juv_Females,  tau_Ms_Juv_Females)
Fs[5] ~ dnorm(Mean_Fs_Juv_Females,  tau_Fs_Juv_Females)

# Subadult_Females
M[6]  ~ dnorm(Mean_Ms_Subadult_Females, tau_Ms_Subadult_Females)
Fs[6] ~ dnorm(Mean_Fs_Subadult_Females, tau_Fs_Subadult_Females)

# Adult_Females
M[7]  ~ dnorm(Mean_Ms_Adult_Females, tau_Ms_Adult_Females)
Fs[7] ~ dnorm(Mean_Fs_Adult_Females, tau_Fs_Adult_Females)

# Resting_Females
M[8]  ~ dnorm(Mean_Ms_Resting_Females, tau_Ms_Resting_Females)
Fs[8] ~ dnorm(Mean_Fs_Resting_Females, tau_Fs_Resting_Females)


#Total mortalities for each stage

for (i in 1:8) {

 Z[i]<- M[i]+Fs[i]

}



###Stochastic fecundity


#Still diffuse
f ~ dunif(16, 41)

# Survival and transition calculations --------------------------------------

# Male survival for stages 1–4
for (i in 1:4) {
  Sa_m[i] <- exp(-Z[i])
}

# Female survival for stages 1 and 5–8
Sa_f[1] <- exp(-Z[1])  # same shared stage for neaonates
for (i in 2:5) {
  Sa_f[i] <- exp(-Z[i + 3])  # maps to Z[5:8]
}

##Probability of growing to the next step for males y_im[i]

  for (i in 1:4) {


    y_im[i] <- ifelse(i == 1, 
                      Sa_m[i], 
                      ((Sa_m[i]^T_im[i]) - (Sa_m[i]^(T_im[i] - 1))) / 
                      (Sa_m[i]^T_im[i] - 1) )
  }

# Transition probabilities G for males

 for (i in 1:4) {
  G_male[i] <- Sa_m[i] * y_im[i]
  p_male[i] <- ifelse(i == 1, 0, Sa_m[i] * (1 - y_im[i]))  # no self-loop for neonate males
  Dp_m[i]   <- 1 - G_male[i] - p_male[i]
}

##Probability of growing to the next step for females y_if[i]

  for (i in 1:5) {

    y_if[i] <- ifelse(i == 1, 
                      Sa_f[i], 
                      ((Sa_f[i]^T_if[i]) - (Sa_f[i]^(T_if[i] - 1))) / 
                      (Sa_f[i]^T_if[i] - 1) )
  }

# Transition probabilities G for females

 for (i in 1:5) {
  G_female[i] <- Sa_f[i] * y_if[i]
  p_female[i] <- ifelse(i == 1, 0, Sa_f[i] * (1 - y_if[i]))  # no self-loop for neonate females
  Dp_f[i]   <- 1 - G_female[i] - p_female[i]
}

#Calculated fecundity only for reproductive stages:

#Males
fec_a[4] <- Sa_m[4] * f * N_t[4] / (N_t[7] + N_t[4]) 

#Females
fec_a[7] <- Sa_f[4] * f * N_t[7] / (N_t[7] + N_t[4]) 


##BUILD Matrix A_m ----------------------------


# Note: p_m/f entries are shifted one step to the right because Jags does not allow to remove the first entry of p (neonates) with p<- p[-1]
 
A_t[1, 4] <- fec_a[4]
A_t[1, 7] <- G_female[4] * fec_a[7]
A_t[2, 1] <- p * G_male[1]
A_t[2, 2] <- p_male[2]
A_t[3, 2] <- G_male[2]
A_t[3, 3] <- p_male[3]
A_t[4, 3] <- G_male[3]
A_t[4, 4] <- p_male[4]
A_t[5, 1] <- (1 - p) * G_male[1]
A_t[5, 5] <- p_female[2]
A_t[6, 5] <- G_female[1]
A_t[6, 6] <- p_female[3]
A_t[7, 6] <- G_female[2]
A_t[7, 8] <- G_female[4]
A_t[8, 7] <- G_female[5]


}  # end model
")
sink()


## Define scenario ---------------------------------------------------------


#Create prior parameter list
# Here you can set how much you want to change mortality values. No need to alter jags code.


#Current ETP fishing mortality

{win.data_current <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females
  
)
}

#Increased neonate/juvenile ETP fishing mortality

{win.data_juvenile_fishing <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*2,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females*2,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males*2,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females
  
)
}


#Extreme neonate/juvenile ETP fishing mortality

{win.data_juvenile_extreme <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*3,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females*3,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males*3,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females
  
)
}

# Overall 25% reduction

{win.data_25 <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo*.75,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*.75,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females*.75,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females*.75,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males*.75,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males*.75,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females*.75,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females*.75,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males*.75,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males*.75,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females*.75,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females*.75,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males*.75,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males*.75,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females*.75,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females*.75
  
)
}



# Overall 50% reduction

{win.data_50 <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo*.5,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*.5,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females*.5,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females*.5,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males*.5,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males*.5,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females*.5,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females*.5,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males*.5,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males*.5,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females*.5,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females*.5,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males*.5,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males*.5,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females*.5,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females*.5
  
)
}

# 75 Reduction of higher elasticity

{win.data_elastic <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females*.25,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females*.25,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males*.25,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males*.25,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females*.25,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females*.25,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males*.25,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males*.25,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females*.25,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females*.25
  
)
}

# 75 Reduction of neonate/juvenile mortality
{win.data_juvenile_protection <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo*0.25,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*0.25,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females*0.25,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females*0.25,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males*0.25,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males*0.25,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females
  
)
}

# 90 Reduction of neonate fishing mortality
{win.data_nursey_protection <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo*0.10,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*0.10,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females
  
)
}

# 75 Reduction of higher elasticity + extreme juvenile mortality

{win.data_elastic_extreme <- list(
  # Population and transition data
  N_t  = N_t,
  T_im = T_im,
  T_if = T_if,
  p    = p,
  
  # NEONATES
  sd_Ms_Neo   = sd_Ms_Neo,
  sd_Fs_Neo   = sd_Fs_Neo,
  Mean_Ms_Neo = Mean_Ms_Neo,
  Mean_Fs_Neo = Mean_Fs_Neo*3,
  
  # JUVENILES — Females
  sd_Ms_Juv_Females   = sd_Ms_Juv_Females,
  sd_Fs_Juv_Females   = sd_Fs_Juv_Females,
  Mean_Ms_Juv_Females = Mean_Ms_Juv_Females,
  Mean_Fs_Juv_Females = Mean_Fs_Juv_Females*3,
  
  # JUVENILES — Males
  sd_Ms_Juv_Males   = sd_Ms_Juv_Males,
  sd_Fs_Juv_Males   = sd_Fs_Juv_Males,
  Mean_Ms_Juv_Males = Mean_Ms_Juv_Males,
  Mean_Fs_Juv_Males = Mean_Fs_Juv_Males*3,
  
  # SUBADULTS — Females
  sd_Ms_Subadult_Females   = sd_Ms_Subadult_Females,
  sd_Fs_Subadult_Females   = sd_Fs_Subadult_Females*.25,
  Mean_Ms_Subadult_Females = Mean_Ms_Subadult_Females,
  Mean_Fs_Subadult_Females = Mean_Fs_Subadult_Females*.25,
  
  # SUBADULTS — Males
  sd_Ms_Subadult_Males   = sd_Ms_Subadult_Males,
  sd_Fs_Subadult_Males   = sd_Fs_Subadult_Males*.25,
  Mean_Ms_Subadult_Males = Mean_Ms_Subadult_Males,
  Mean_Fs_Subadult_Males = Mean_Fs_Subadult_Males*.25,
  
  # ADULTS — Females
  sd_Ms_Adult_Females   = sd_Ms_Adult_Females,
  sd_Fs_Adult_Females   = sd_Fs_Adult_Females*.25,
  Mean_Ms_Adult_Females = Mean_Ms_Adult_Females,
  Mean_Fs_Adult_Females = Mean_Fs_Adult_Females*.25,
  
  # ADULTS — Males
  sd_Ms_Adult_Males   = sd_Ms_Adult_Males,
  sd_Fs_Adult_Males   = sd_Fs_Adult_Males*.25,
  Mean_Ms_Adult_Males = Mean_Ms_Adult_Males,
  Mean_Fs_Adult_Males = Mean_Fs_Adult_Males*.25,
  
  # RESTING — Females only
  sd_Ms_Resting_Females   = sd_Ms_Resting_Females,
  sd_Fs_Resting_Females   = sd_Fs_Resting_Females*.25,
  Mean_Ms_Resting_Females = Mean_Ms_Resting_Females,
  Mean_Fs_Resting_Females = Mean_Fs_Resting_Females*.25
  
)
}



#Define output to save un model object
params <- c("M", "Sa_m", "y_im", "G_male", "p_male", 
            "Sa_f", "y_if", "G_female", "p_female", 
            "A_t", "Fs", "Z", "f" )

#Run model

##**IMPORTANT**##
##change the first line (data) to project a different scenario#
#Store outputs as needed. Im sure theres a way to build a complex function or loop to do it in one go. But I think its better to monitor each scenario step by step

out_matrix <- jags(
                    data = win.data_current,
                   #data = win.data_25,
                   #data = win.data_50,
                   #data = win.data_elastic,
                   #data = win.data_elastic_extreme,
                   #data = win.data_juvenile_protection, 
                   #data = win.data_juvenile_fishing,
                   #data = win.data_juvenile_extreme,
                   #data = win.data_nursey_protection,
                   inits = NULL,
                   parameters.to.save = params,
                   model.file = "MatrixModel_scenario.txt",
                   n.chains = 2,
                   n.iter = 10000,
                   n.burnin = 1000,
                   n.thin = 1,
                   DIC=F)


#MCMC vis to quickly check traveplots and densities
MCMCtrace(out_matrix, 
          params = c("M"), 
          pdf = FALSE)


## Extract summaries from simmulations -------------------------------------

#Store Sample Distributions of mortalities


# Total mortality
Z_all <- out_matrix$BUGSoutput$sims.list$Z  
#Natural mortality
M_all <- out_matrix$BUGSoutput$sims.list$M  
#Fishing mortality
Fs_all <- out_matrix$BUGSoutput$sims.list$Fs  



# Calculate the mean mortality across stages for each draw
mean_Z <- rowMeans(Z_all)
mean_M <- rowMeans(M_all)
mean_Fs <- rowMeans(Fs_all)


# Overall fecundity
f_all <- out_matrix$BUGSoutput$sims.list$f  

# Vectorize fecundity. single value for each draw
mean_f <- rowMeans(f_all)



#Extract matrix entries

# Extract all iterations of monitored parameters
A_raw <- out_matrix$BUGSoutput$sims.matrix

# Get the number of MCMC samples
n_iter <- nrow(A_raw)

# Initialize a list to store matrices per iteration
A_t_list <- vector("list", n_iter)

# Loop over iterations and fill in the matrix
for (i in 1:n_iter) {
  A <- matrix(0, nrow = 8, ncol = 8)
  
  A[1, 4] <- A_raw[i, "A_t[1,4]"]
  A[1, 7] <- A_raw[i, "A_t[1,7]"]
  A[2, 1] <- A_raw[i, "A_t[2,1]"]
  A[2, 2] <- A_raw[i, "A_t[2,2]"]
  A[3, 2] <- A_raw[i, "A_t[3,2]"]
  A[3, 3] <- A_raw[i, "A_t[3,3]"]
  A[4, 3] <- A_raw[i, "A_t[4,3]"]
  A[4, 4] <- A_raw[i, "A_t[4,4]"]
  A[5, 1] <- A_raw[i, "A_t[5,1]"]
  A[5, 5] <- A_raw[i, "A_t[5,5]"]
  A[6, 5] <- A_raw[i, "A_t[6,5]"]
  A[6, 6] <- A_raw[i, "A_t[6,6]"]
  A[7, 6] <- A_raw[i, "A_t[7,6]"]
  A[7, 8] <- A_raw[i, "A_t[7,8]"]
  A[8, 7] <- A_raw[i, "A_t[8,7]"]
  A_t_list[[i]] <- A
}


#Mean Matrix
# Compute mean matrix from list of matrices
A_t_mean <- Reduce("+", A_t_list) / n_iter

#Lambdas

# Calculate population growth rate (lambda)

lambdas <- sapply(A_t_list, function(A) {
  eigvals <- eigen(A)$values
  Re(eigvals[which.max(Re(eigvals))])  # Get the dominant real part
})

x<-sapply(A_t_list, function(A) {
  eigvals <- eigen(A)$values
})

#Mean Mlambda
# Compute mean matrix from list of matrices
mean_lambda <- mean(lambdas)
sd(lambdas)
#n lambdas, just to chek its one per itteration

length(lambdas)

#Lambdas with a value larger than 1

length(lambdas[lambdas>1])

#Proportion


length(lambdas[lambdas>1])/length(lambdas)

mean(lambdas)
sd(lambdas)
##Append mortalities and lambda to matrix list

A_t_list$M<-mean_M
A_t_list$Fs<-mean_Fs
A_t_list$Z<-mean_Z
A_t_list$fecundity<-mean_f
A_t_list$lambda<-lambdas

# 1. Neonates
A_t_list$M_Neonates         <- M_all[, 1]
A_t_list$Fs_Neonates        <- Fs_all[, 1]
A_t_list$Z_Neonates         <- Z_all[, 1]

# 2. Juvenile Males
A_t_list$M_Juvenile_Males   <- M_all[, 2]
A_t_list$Fs_Juvenile_Males  <- Fs_all[, 2]
A_t_list$Z_Juvenile_Males   <- Z_all[, 2]

# 3. Subadult Males
A_t_list$M_Subadult_Males   <- M_all[, 3]
A_t_list$Fs_Subadult_Males  <- Fs_all[, 3]
A_t_list$Z_Subadult_Males   <- Z_all[, 3]

# 4. Adult Males
A_t_list$M_Adult_Males      <- M_all[, 4]
A_t_list$Fs_Adult_Males     <- Fs_all[, 4]
A_t_list$Z_Adult_Males      <- Z_all[, 4]

# 5. Juvenile Females
A_t_list$M_Juvenile_Females <- M_all[, 5]
A_t_list$Fs_Juvenile_Females<- Fs_all[, 5]
A_t_list$Z_Juvenile_Females <- Z_all[, 5]

# 6. Subadult Females
A_t_list$M_Subadult_Females <- M_all[, 6]
A_t_list$Fs_Subadult_Females<- Fs_all[, 6]
A_t_list$Z_Subadult_Females <- Z_all[, 6]

# 7. Adult Females
A_t_list$M_Adult_Females    <- M_all[, 7]
A_t_list$Fs_Adult_Females   <- Fs_all[, 7]
A_t_list$Z_Adult_Females    <- Z_all[, 7]

# 8. Resting Females
A_t_list$M_Resting_Females  <- M_all[, 8]
A_t_list$Fs_Resting_Females <- Fs_all[, 8]
A_t_list$Z_Resting_Females  <- Z_all[, 8]


## Project population (with burn in) ----------------------------------------------------------

# PARAMETERS

n_iter <- length(A_t_list$lambda) #Extract number of model itterations
n_years_b <- 36  #Years to project population
n_stages <- length(N_t)  #Number of lifestages
N_t <- c(10000, 1000, 500, 400, 1000, 500, 200, 200) #Initial population


# Preallocate empty list to store results
popb_summary_list <- vector("list", n_iter)


# LOOP: project total population over time for each iteration
for (i in 1:n_iter) {
  N_temp <- N_t  # use a temporary variable!
  popb_iter <- numeric(n_years_b + 1)
  popb_iter[1] <- sum(N_temp)
  popb_sub <- numeric(n_years_b + 1)
  popb_sub[1] <- sum(N_temp[c(3,4,6,7,8)])
  
  for (t in 1:n_years_b) {
    N_temp <- A_t_list[[i]] %*% N_temp
    popb_iter[t + 1] <- sum(N_temp)
    popb_sub[t + 1] <- sum(N_temp[c(3,4,6,7,8)])
  }
  
  
  # STORE RESULTS 
  popb_summary_list[[i]] <- data.frame(
    Iteration = i,
    Year = 0:n_years_b,
    
    # Total populations
    Population = popb_iter,
    SubPop = popb_sub,
    
    Natural = A_t_list$M[i],
    Fishing= A_t_list$Fs[i],
    Total = A_t_list$Z[i],
    
    # Natural mortality by stage
    Natural_Neonates         = A_t_list$M_Neonates[i],
    Natural_Juvenile_Males   = A_t_list$M_Juvenile_Males[i],
    Natural_Subadult_Males   = A_t_list$M_Subadult_Males[i],
    Natural_Adult_Males      = A_t_list$M_Adult_Males[i],
    Natural_Juvenile_Females = A_t_list$M_Juvenile_Females[i],
    Natural_Subadult_Females = A_t_list$M_Subadult_Females[i],
    Natural_Adult_Females    = A_t_list$M_Adult_Females[i],
    Natural_Resting_Females  = A_t_list$M_Resting_Females[i],
    
    # Fishing mortality by stage
    Fishing_Neonates         = A_t_list$Fs_Neonates[i],
    Fishing_Juvenile_Males   = A_t_list$Fs_Juvenile_Males[i],
    Fishing_Subadult_Males   = A_t_list$Fs_Subadult_Males[i],
    Fishing_Adult_Males      = A_t_list$Fs_Adult_Males[i],
    Fishing_Juvenile_Females = A_t_list$Fs_Juvenile_Females[i],
    Fishing_Subadult_Females = A_t_list$Fs_Subadult_Females[i],
    Fishing_Adult_Females    = A_t_list$Fs_Adult_Females[i],
    Fishing_Resting_Females  = A_t_list$Fs_Resting_Females[i],
    
    # Total mortality by stage
    Total_Neonates           = A_t_list$Z_Neonates[i],
    Total_Juvenile_Males     = A_t_list$Z_Juvenile_Males[i],
    Total_Subadult_Males     = A_t_list$Z_Subadult_Males[i],
    Total_Adult_Males        = A_t_list$Z_Adult_Males[i],
    Total_Juvenile_Females   = A_t_list$Z_Juvenile_Females[i],
    Total_Subadult_Females   = A_t_list$Z_Subadult_Females[i],
    Total_Adult_Females      = A_t_list$Z_Adult_Females[i],
    Total_Resting_Females    = A_t_list$Z_Resting_Females[i],
    
    # Vital rates
    Fecundity = A_t_list$fecundity[i],
    lambda    = A_t_list$lambda[i]
  )
}

# Combine all iterations into one data frame
popb_proj_df <- do.call(rbind, popb_summary_list)

popb_summary <- popb_proj_df %>%
  group_by(Year) %>%
  summarise(
    mean = mean(Population),
    lower = quantile(Population, 0.025),
    upper = quantile(Population, 0.975)
  )



#Parameterize population


#Long format of filtered estimates
sim_pop_long <- pop_proj_df%>%
  dplyr::select(-Population, -SubPop,-Year) %>%
  pivot_longer(
    cols = -c(Iteration),  # all columns to pivot
    names_to = "Variable",   # name for the new column holding original column names
    values_to = "Estimate"   # name for the new column holding the values
  )


# group and get minuimum and maximum parameter value per estimate

parameter_summary<-sim_pop_long%>%group_by(Variable)%>%
  summarize(mean=mean(Estimate),
            minimum=min(Estimate),
            maximum=max(Estimate),
            SD=sd(Estimate))%>%
  ungroup()

print(parameter_summary)

#write.csv(parameter_summary,"ParameterEstimates.csv")


## Spaghetti plots -----


#ease visualization trimming extreme values
popb_proj_dftrim<-popb_proj_df[popb_proj_df$Population<10000 & popb_proj_df$Population> 0,]

#this one is for lambda critical only
popb_proj_dftrim2<-popb_proj_df[popb_proj_df$lambda<1.025 & popb_proj_df$lambda> 0.975,]


popb_summary_trim <- popb_proj_dftrim2 %>%
  group_by(Year) %>%
  summarise(
    mean = mean(Population),
    lower = quantile(Population, 0.025),
    upper = quantile(Population, 0.975)
  )


#Violin plot -------------------------------------------------------------

 # 
 # current_violin_df<-popb_proj_df|>
 #   group_by(Iteration)|>
 #   summarise(Growth_Rate=first(lambda))

# elastic_violin_df<-popb_proj_df|>
#   group_by(Iteration)|>
#   summarise(Growth_Rate=first(lambda))

 # elastic_extreme_violin_df<-popb_proj_df|>
 #  group_by(Iteration)|>
 #  summarise(Growth_Rate=first(lambda))

# juveniles_protected_violin_df<-popb_proj_df|>
#   group_by(Iteration)|>
#   summarise(Growth_Rate=first(lambda))

 # juveniles_fishing_violin_df<-popb_proj_df|>
 #   group_by(Iteration)|>
 #   summarise(Growth_Rate=first(lambda))

# juveniles_extreme_violin_df<-popb_proj_df|>
#   group_by(Iteration)|>
#   summarise(Growth_Rate=first(lambda))


 # nursery_protected_violin_df<-popb_proj_df|>
 #   group_by(Iteration)|>
 #   summarise(Growth_Rate=first(lambda))
 # 
 # uniform_smallprotection_violin_df<-popb_proj_df|>
 #   group_by(Iteration)|>
 #   summarise(Growth_Rate=first(lambda))

# uniform_protection_violin_df<-popb_proj_df|>
#   group_by(Iteration)|>
#   summarise(Growth_Rate=first(lambda))


 
 nofishing_violin_df<-popb_proj_df|>
   group_by(Iteration)|>
   summarise(Growth_Rate=first(lambda))




#shades or ablines

v_plot_full_Scenario <- #current_violin_df|>
                        #elastic_violin_df|>
                        #elastic_extreme_violin_df|>
                        #juveniles_protected_violin_df|>
                        #juveniles_fishing_violin_df|>
                        #juveniles_extreme_violin_df|>
                        #nursery_protected_violin_df|>
                        #uniform_protection_violin_df|>
                         uniform_smallprotection_violin_df|>
  ggplot(aes(x = "", y = Growth_Rate)) +
  geom_violin(fill = viridis(12, option = "mako")[10]) +
  geom_boxplot(width = 0.15,
               color = viridis(12, option = "mako")[1],  
               fill  = viridis(12, option = "mako")[8],
               linewidth = 1.1) +
  geom_point(
    stat = "summary",
    fun = "mean",
    color = viridis(12, option = "mako")[1],  # dark cyan
    fill  = viridis(12, option = "mako")[6],
    shape = 21,
    size = 10
  ) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "red",
             linewidth = 1)+
  theme_classic(base_size = 18) +
  theme() +
  labs(y = "Lambda", x="Hammerheads")


print (v_plot_full_Scenario)

#here you can change the name and such if you want to save individual plots. Let a couple examples below.

ggsave("ViolinPlot_Full_Current.jpg", plot = v_plot_full_Scenario,
       width = 12, height = 10, dpi = 300)

  #For example. Saving a plot for the overall 50% reduction df
# ggsave("ViolinPlot_Full_50Reduced.jpg", plot = v_plot_full_Reduced_Fishing,
#        width = 12, height = 10, dpi = 300)



# Combined Violin ----------------------------------------------------------

#This creates Figure 5. Merging all violin DFs into a single object and ploting

current_violin_df$Model <- "Current"
elastic_violin_df$Model <- "Adults_protected"
elastic_extreme_violin_df$Model <- "Extreme Elastic"
juveniles_protected_violin_df$Model <- "Juveniles_protected"
juveniles_fishing_violin_df$Model <- "2xJuvenile_fishing"
juveniles_extreme_violin_df$Model <-"3xJuvenile_fishing"
nursery_protected_violin_df$Model <- "Neonates_protected"
uniform_smallprotection_violin_df$Model <- "25%_reduction"
uniform_protection_violin_df$Model <- "50%_reduction"
nof_violin_df$Model <- "No_fishing"


Violin_df<-rbind(current_violin_df,
                 elastic_violin_df,
                 elastic_extreme_violin_df,
                 juveniles_protected_violin_df,
                 juveniles_fishing_violin_df,
                 juveniles_extreme_violin_df,
                 nursery_protected_violin_df,
                 uniform_smallprotection_violin_df,
                 uniform_protection_violin_df,
                 nof_violin_df)

#Arrange 

Violin_df$Model <- factor(
  Violin_df$Model,
  levels = c(
    "Current",
    "2xJuvenile_fishing",
    "3xJuvenile_fishing",
    "Extreme Elastic",
    "Neonates_protected",
    "Juveniles_protected",
    "25%_reduction",
    "50%_reduction",
    "Adults_protected",
    "No_fishing"
  )
)




v_plot_Scenarios <-ggplot(Violin_df, aes(x = Model, y = Growth_Rate)) +
  geom_violin(aes(fill =Model)) +
  scale_fill_viridis_d(option = "mako",begin = 0.3,alpha=0.5) +   # mako palette for stages
  geom_boxplot(width = 0.15,
               color = viridis(12, option = "mako")[1],  
               fill  = viridis(12, option = "mako")[8],
               linewidth = 1.1) +
  geom_point(
    stat = "summary",
    fun = "mean",
    color = viridis(12, option = "mako")[1],  # dark cyan
    fill  = viridis(12, option = "mako")[6],
    shape = 21,
    size = 10
  ) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "red",
             linewidth = 1)+
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom",
        axis.text.x.bottom = element_text(size=12)) +
  labs(y = "Lambda", x="Model")


#visualize final plot

plot(v_plot_Scenarios)


#save
ggsave("ViolinPlot_Model_Scenarios.jpg", plot = v_plot_Scenarios,
       width = 16, height = 12, dpi = 300)

#estimate mean lambdas per scenario. This is a suppelmentary table at the time of writing. 
summary_Violin_df<-Violin_df|>group_by(Model)|>
  summarise(mean_lambda = mean(Growth_Rate),
            sd_lambda=sd(Growth_Rate),
            percent = sum(Growth_Rate >= 1) / n(),
            n=n())

write.csv(summary_Violin_df,"growth_by_model.csv")



# Single unbiased model ---------------------------------------------------
#Mortalities
M<-c(0.7,0.38,0.38,0.19,0.19,0.19,0.19,0.19)
F_low<-c(0.1,0.15,0.15,0.15,0.15,0.15,0.15,0.15)
F_mean<-c(0.1,0.25,0.25,0.25,0.25,0.25,0.25,0.25)
F_high<-c(0.1,0.4,0.4,0.4,0.4,0.4,0.4,0.4)
F_juvs<-c(0.1,0.0,0.0,0.25,0.0,0.0,0.25,0.25)
F_adults<-c(0.1,0.25,0.25,0.0,0.25,0.25,0.0,0.0)

f<-19

#previously created parameters
T_im
T_if
N_t
p
Sa_f<-c()

scenarios <- list(
  nf     = list(F = rep(0, 8),  label = "No fishing"),
  low    = list(F = F_low,      label = "Low fishing"),
  mean   = list(F = F_mean,     label = "Mean fishing"),
  high   = list(F = F_high,     label = "High fishing"),
  juvs   = list(F = F_juvs,     label = "Juveniles protected"),
  adults = list(F = F_adults,   label = "Adults protected")
)

#Matrix building function

build_matrix <- function(M, F, f, T_im, T_if, N_t, p) {
  
  Z <- M + F
  
  # Survival
  Sa_m <- exp(-Z[1:4])
  Sa_f <- c(exp(-Z[1]), exp(-Z[5:8]))
  
  # Growth probs
  y_im <- sapply(1:4, function(i)
    if (i == 1) Sa_m[i]
    else ((Sa_m[i]^T_im[i]) - (Sa_m[i]^(T_im[i]-1))) / (Sa_m[i]^T_im[i] - 1)
  )
  
  y_if <- sapply(1:5, function(i)
    if (i == 1) Sa_f[i]
    else ((Sa_f[i]^T_if[i]) - (Sa_f[i]^(T_if[i]-1))) / (Sa_f[i]^T_if[i] - 1)
  )
  
  # Transitions
  G_male   <- Sa_m * y_im
  p_male   <- pmax(0, Sa_m - y_im)
  
  G_female <- Sa_f * y_if
  p_female <- pmax(0, Sa_f - y_if)
  
  # Fecundity
  fec <- rep(0, 8)
  fec[4] <- Sa_m[4] * f * N_t[4] / (N_t[7] + N_t[4])
  fec[7] <- Sa_f[4] * f * N_t[7] / (N_t[7] + N_t[4])
  
  # Matrix
  A <- matrix(0, 8, 8)
  
  A[1,4] <- fec[4]
  A[1,7] <- G_female[4] * fec[7]
  A[2,1] <- p * G_male[1]
  A[2,2] <- p_male[2]
  A[3,2] <- G_male[2]
  A[3,3] <- p_male[3]
  A[4,3] <- G_male[3]
  A[4,4] <- p_male[4]
  A[5,1] <- (1-p) * G_male[1]
  A[5,5] <- p_female[2]
  A[6,5] <- G_female[1]
  A[6,6] <- p_female[3]
  A[7,6] <- G_female[2]
  A[7,8] <- G_female[4]
  A[8,7] <- G_female[5]
  
  list(
    A = A,
    lambda = max(Re(eigen(A)$values),
    Sa_f=Sa_f)
  )
}


#Model matrix per scenario

model_results <- lapply(scenarios, function(s) {
  build_matrix(M, s$F, f, T_im, T_if, N_t, p)
})



model_results$nf$lambda
model_results$low$lambda
model_results$mean$lambda
model_results$high$lambda
model_results$juvs$lambda
model_results$adults$lambda


#Project scenarios

N0 <- c(10000, 1000, 500, 400, 1000, 500, 200, 200)
n_years<-35

project_population <- function(A, N0, n_years, scenario_name) {
  
  n_stages <- length(N0)
  
  N_mat <- matrix(NA, nrow = n_stages, ncol = n_years + 1)
  N_mat[, 1] <- N0
  
  for (t in 1:n_years) {
    N_mat[, t + 1] <- A %*% N_mat[, t]
  }
  
  pop_total <- colSums(N_mat)
  
  data.frame(
    Year = 0:n_years,
    Population = pop_total,
    Scenario = scenario_name
  )
}

pop_list <- mapply(function(res, s) {
  project_population(res$A, N0, n_years, s$label)
}, model_results, scenarios, SIMPLIFY = FALSE)

pop_all_long <- dplyr::bind_rows(pop_list)

lambda_df <- data.frame(
  Scenario = c("No fishing", "Low fishing", "Mean fishing", "High fishing","Adults protected", "Juveniles protected"),
  lambda = c(1.07, 0.89, 0.82, 0.71, 0.93,0.95))

pop_all_long <- pop_all_long |>
  dplyr::left_join(lambda_df, by = "Scenario") |>
  dplyr::mutate(
    Scenario_label = paste0(
      Scenario, " (\u03BB = ", sprintf("%.2f", lambda), ")"
    )
  )


pop_plot <- ggplot(
  pop_all_long,
  aes(x = Year, y = Population, color = Scenario_label)) +
  geom_line(linewidth = 2,aes(linetype = Scenario)) +
  guides(linetype = "none",color = guide_legend(override.aes = list(linewidth = 4)))+
  scale_color_viridis_d(
    option = "mako",
    end = 0.85,
    name = "Scenario") +
  labs(
    x = "Year",
    y = "Total population size") +
  theme_minimal(base_size = 25) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 20),
    legend.text  = element_text(size = 15),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"))

ggsave("Supplementary_scenarios.jpg", plot = 
         pop_plot,
       width = 12, height = 10, dpi = 300)

