################################################################################
# Monte-Carlo Simulation Study to Verify Bootstrap
# Authors: Ties Bos, Maciej Korek and Joshua Smilde
################################################################################
# Used packages
library("dplyr")
library("stats")
library("extraDistr")
library("tsoutliers")
library("goft")

# Import the data:
scan.data <- read.csv("ScanRecords.csv", header=TRUE)

# select the data of type 2 patients:
type_2_ind <- scan.data$PatientType == "Type 2"
type_2_data <- scan.data[type_2_ind, ]

# Assuming the duration of type 2 patients follow a gamma distribution, we 
# require the two parameter estimates.
mean_duration_2 <- mean(type_2_data$Duration)
sd_duration_2 <- sd(type_2_data$Duration)

beta_param <- sd_duration_2^2/mean_duration_2
alpha_param <- mean_duration_2/beta_param

# Assuming the time between two consecutive calls for type 2 patients is normally
# distributed, we need the mean and standard deviation estimates:
time_patient2 <- rep(NA, nrow(type_2_data)-1)
for(i in 2:nrow(type_2_data)){
  time <- type_2_data$Time[i]-type_2_data$Time[i-1]
  if(time < 0){
    time_patient2[i-1] <- time+9
  } else {
    time_patient2[i-1] <- time
  }
}

mean_time_2 <- mean(time_patient2)
sd_time_2 <- sd(time_patient2)

# We use a Monte-Carlo approach to verify our Bootstrap approach for type 2.
# We focus on the established confidence intervals, surrounding means and 
# standard deviation.

# For starters, we make a function of our bootstrap, with input:
#   - A sample for the duration X1
#   - A sample for the time between two calls X2
#   - Bootstrap sample size B
#   - Level of the test alpha


boot_type2 <- function(X1, X2, B, alpha){
  n1 <- length(X1)
  n2 <- length(X2)
  
  mean.X1 <- mean(X1)
  sd.X1 <- sd(X1)
  
  mean.X2 <- mean(X2)
  sd.X2 <- sd(X2)
  
  # We will store the statistics in the vectors:
  # T-statistics for CI surrounding the mean:
    T_stat_dur_2 <- rep(NA, B)
    T_stat_time_2 <- rep(NA,B)
  # The statistics for the CI surrounding the variances:
    Q_stat_dur_2 <- rep(NA,B)
    Q_stat_time_2 <- rep(NA,B)
    
  # We run the bootstrap:
    for(b in 1:B){
      J1 <- sample(1:n1, replace=TRUE)
      X1.star <- X1[J1]
      mean.X1.star <- mean(X1.star)
      sd.X1.star <- sd(X1.star)
      T_stat_dur_2[b] <- sqrt(n1)*(mean.X1.star-mean.X1)/sd.X1.star
      Q_stat_dur_2[b] <- (n1-1)*sd.X1.star^2/(sd.X1^2)
      
      J2 <- sample(1:n2, replace=TRUE)
      X2.star <- X2[J2]
      mean.X2.star <- mean(X2.star)
      sd.X2.star <- sd(X2.star)
      T_stat_time_2[b] <- sqrt(n2)*(mean.X2.star-mean.X2)/sd.X2.star
      Q_stat_time_2[b] <- (n2-1)*sd.X2.star^2/(sd.X2^2)
    }
    
    # Critical values
    cv_T_dur_2 <- quantile(T_stat_dur_2, probs=c(alpha/2, 1-alpha/2))
    cv_T_time_2 <- quantile(T_stat_time_2, probs=c(alpha/2, 1-alpha/2))
    
    cv_Q_dur_2 <- quantile(Q_stat_dur_2, probs=c(alpha/2, 1-alpha/2))
    cv_Q_time_2 <- quantile(Q_stat_time_2, probs=c(alpha/2, 1-alpha/2))
    
    # two-sided 1-alpha-confidence-intervals
    CIs <- matrix(NA, nrow=4, ncol=2, 
                  dimnames = list(c("CI Mean Duration", "CI Variance Duration",
                                    "CI Mean Time", "CI Variance Time"), 
                                  c("Lowerbound", "Upperbound")))
    CIs[1,1] <- mean.X1-cv_T_dur_2[2]*sd.X1/sqrt(n1)
    CIs[1,2] <- mean.X1-cv_T_dur_2[1]*sd.X1/sqrt(n1)
    CIs[2,1] <- (n1-1)*sd.X1^2/cv_Q_dur_2[2]
    CIs[2,2] <- (n1-1)*sd.X1^2/cv_Q_dur_2[1]
    
    CIs[3,1] <- mean.X2-cv_T_time_2[2]*sd.X2/sqrt(n2)
    CIs[3,2] <- mean.X2-cv_T_time_2[1]*sd.X2/sqrt(n2)
    CIs[4,1] <- (n2-1)*sd.X2^2/cv_Q_time_2[2]
    CIs[4,2] <- (n2-1)*sd.X2^2/cv_Q_time_2[1]
    
    return(CIs)
}

# Monte-Carlo Simulations based on these Confidence-interval

# Parametric approach
#   If we assume the time between two calls is normal and the duration is gamma
#   distributed (with the specified parameter):
#     The number of MC-simulation iterations
no.sims <- 900
#     The number of bootstrap draws
B <- 1000
#     The rejection level of interest:
alpha <- 0.05

# We choose a variety of sample sizes if we want to make a plot:
# n <- seq(10, 500, 50)
# or, more relevant, we use the sample size to find more accurate information
# about the dataset at hand:
n <- nrow(type_2)


# Save rejection probabilities for all n in a matrix:
Rej_prob_matrix <- matrix(NA, nrow=0, ncol=5,
                          dimnames = list(c(), c("Sample Size","Rej. Prob Mean Duration", 
                                                 "Rej. Prob. Variance Duration",
                                                 "Rej. Prob. Mean Time", 
                                                 "Rej. Prob. Variance Time")))


for(i in 1:length(n)){
  n_i <- n[i]
  #     Store rejections in the following matrix:
  Reject.mat <- matrix(NA, nrow=no.sims, ncol=4, 
                       dimnames = list(c(), c("Reject Mean Duration", 
                                              "Reject Variance Duration",
                                              "Reject Mean Time", 
                                              "Reject Variance Time")))
  for(i in 1:no.sims){
    X_dur <- rgamma(n_i, alpha_param, 1/beta_param)
    
    X_time <- rnorm(n_i-1, mean_time_2, sd_time_2)
    
    Conf_Ints <- boot_type2(X_dur, X_time, B, alpha)
    
    Reject.mat[i,1] <- 1*(mean_duration_2 < Conf_Ints[1,1]) + 1*(mean_duration_2> Conf_Ints[1,2])
    Reject.mat[i,2] <- 1*(sd_duration_2^2 < Conf_Ints[2,1]) + 1*(sd_duration_2^2 > Conf_Ints[2,2])
    Reject.mat[i,3] <- 1*(mean_time_2 < Conf_Ints[3,1]) + 1*(mean_time_2 > Conf_Ints[3,2])
    Reject.mat[i,4] <- 1*(sd_time_2^2 < Conf_Ints[4,1]) + 1*(sd_time_2^2 > Conf_Ints[4,2])
  }
  
  rej_prob <- c(n_i, colMeans(Reject.mat))
  Rej_prob_matrix <- rbind(Rej_prob_matrix, rej_prob)
  print(n_i) # Check progress
}

# We can make nice plots (if n was chosen as a vector of values):
plot(n, Rej_prob_matrix[,2], type="l", col="black", 
     ylim=c(min(Rej_prob_matrix[, 2:5]), max(Rej_prob_matrix[, 2:5])),
     ylab="Rejection Probabilities", xlab="Sample Size",
     main="Rejection Probabilities from Bootstrapped 95%-Confidence Intervals")
lines(n, Rej_prob_matrix[,3], col="black", lty=3)
lines(n, Rej_prob_matrix[,4], col="black", lty=5)
lines(n, Rej_prob_matrix[,5], col="black", lty=6)
abline(a=alpha, b=0, col="red")
legend("topright", c("Mean Duration", "Variance Duration", "Mean Time", 
                     "Variance Time", "Rejection Level"),
       bty="n", lty=c(1, 3, 5, 6, 1), col=c("black", "black", "black", "black", "red"))


# Non-Parametric Approach 
#   --> THE RESULTS OF THIS METHOD ARE ACTUALLY USED IN THE REPORT (SECTION 2.4)
# Alternatively, we can take on a more non-parametric approach by sampling from the 
# empirical distribution function:
#     The number of MC-simulations
no.sims <- 900
#     The number of Bootstrap draws
B <- 1000
# Rejection level of the test:
alpha <- 0.05

# We choose a variety of sample sizes if we want to make a plot:
# n <- seq(10, 500, 50)
# or, more relevant, we use the sample size to find more accurate information
# about the dataset at hand:
n <- nrow(type_2_data)
# (Both the vector or this single value can be used!)

# Save rejection probabilities for all n in a matrix:
Rej_prob_matrix_2 <- matrix(NA, nrow=0, ncol=5,
                          dimnames = list(c(), c("Sample Size","Rej. Prob Mean Duration", 
                                                 "Rej. Prob. Variance Duration",
                                                 "Rej. Prob. Mean Time", 
                                                 "Rej. Prob. Variance Time")))

# We allow that n is a vector of sample sizes for a nice plot
for(i in 1:length(n)){
  n_i <- n[i]
  #     Store rejections in the following matrix:
  Reject.mat <- matrix(NA, nrow=no.sims, ncol=4, 
                       dimnames = list(c(), c("Reject Mean Duration", 
                                              "Reject Variance Duration",
                                              "Reject Mean Time", 
                                              "Reject Variance Time")))
  for(i in 1:no.sims){
    J_dur <- sample(1:length(type_2_data$Duration), n_i, replace=TRUE)
    X_dur <- type_2_data$Duration[J_dur]
    
    J_time <- sample(1:length(time_patient2), n_i-1, replace=TRUE)
    X_time <- time_patient2[J_time]
    
    Conf_Ints <- boot_type2(X_dur, X_time, B, alpha)
    
    Reject.mat[i,1] <- 1*(mean_duration_2 < Conf_Ints[1,1]) + 1*(mean_duration_2> Conf_Ints[1,2])
    Reject.mat[i,2] <- 1*(sd_duration_2^2 < Conf_Ints[2,1]) + 1*(sd_duration_2^2 > Conf_Ints[2,2])
    Reject.mat[i,3] <- 1*(mean_time_2 < Conf_Ints[3,1]) + 1*(mean_time_2 > Conf_Ints[3,2])
    Reject.mat[i,4] <- 1*(sd_time_2^2 < Conf_Ints[4,1]) + 1*(sd_time_2^2 > Conf_Ints[4,2])
  }
  
  rej_prob <- c(n_i, colMeans(Reject.mat))
  Rej_prob_matrix_2 <- rbind(Rej_prob_matrix_2, rej_prob)
  print(n_i)
}

# We can make nice plots (if n was chosen as an array of values):
plot(n, Rej_prob_matrix_2[,2], type="l", col="black", 
     ylim=c(min(Rej_prob_matrix_2[, 2:5]), max(Rej_prob_matrix_2[, 2:5])),
     ylab="Rejection Probabilities", xlab="Sample Size",
     main="Rejection Probabilities from Bootstrapped 95%-Confidence Intervals")
lines(n, Rej_prob_matrix_2[,3], col="black", lty=3)
lines(n, Rej_prob_matrix_2[,4], col="black", lty=5)
lines(n, Rej_prob_matrix_2[,5], col="black", lty=6)
abline(a=alpha, b=0, col="red")
legend("topright", c("Mean Duration", "Variance Duration", "Mean Time", 
                     "Variance Time", "Rejection Level"),
       bty="n", lty=c(1, 3, 5, 6, 1), col=c("black", "black", "black", "black", "red"))














