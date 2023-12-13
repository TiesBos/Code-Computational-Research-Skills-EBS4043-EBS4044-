################################################################################
# Statistical Analysis Part - Computational Research Skills
# Authors: Ties Bos, Maciej Korek and Joshua Smilde
#
################################################################################

############ Packages ##########################################################
library("dplyr")
library("stats")
library("extraDistr")
library("tsoutliers")
library("goft")

# Set a seed:
set.seed(1)

############ Import The Dataset ################################################
scan.data <- read.csv("ScanRecords.csv", header=TRUE)

# We check whether the data has been imported nicely:
  head(scan.data)

############ Key Statistics (2.1) ##############################################
  summary(scan.data)
  
# For the average number of calls, we need the number of registered days in the set:
  dates <- scan.data$Date
# Select only the dates that are not duplicated
  unique_dates <- dates[!duplicated(dates)]
  n_dates <- length(unique_dates)
  avg_no_calls_general <- nrow(scan.data)/n_dates
  
# Split up the dataset in two parts for the types of patients:
  type_1_ind <- scan.data$PatientType == "Type 1"
  type_2_ind <- scan.data$PatientType == "Type 2"
# We select the rows of the dataset of type 1 and type 2 patients:
  type_1_data <- scan.data[type_1_ind, ]
  type_2_data <- scan.data[type_2_ind, ]
# No of type 1 patients:
  n_1 <- nrow(type_1_data)
  n_2 <- nrow(type_2_data)
  
# Then, the summary statistics per group:
  summary(type_1_data)
  summary(type_2_data)
# Average number of calls per day for both groups:
  avg_no_calls_1 <- n_1/n_dates
  avg_no_calls_2 <- n_2/n_dates

############ Duration Dynamics (2.2) ###########################################
# Sample means:
  mean_duration_1 <- mean(type_1_data$Duration)
  mean_duration_2 <- mean(type_2_data$Duration)
# Sample Standard Deviations:
  sd_duration_1 <- sd(type_1_data$Duration)
  sd_duration_2 <- sd(type_2_data$Duration)

# We plot the histograms alongside the fitted distributions of the duration.
  # Type 1 Patients:
  hist(type_1_data$Duration, col="deepskyblue4", prob=TRUE, 
       main="MRI Scan Duration of Type 1 Patients", xlab="MRI Scan Duration", 
       breaks=20)
  # plot the density of a normal distribution with sample mean and std. dev.:
  x_1 <- seq(min(type_1_data$Duration), max(type_1_data$Duration), 0.01)
  lines(x_1, dnorm(x_1, mean=mean_duration_1, sd=sd_duration_1), col="#CC0000", 
        lwd=2.4)
  
  # Type 2 Patients:
  hist(type_2_data$Duration, breaks=10, probability=TRUE, col="deepskyblue4", 
       main="MRI Scan Duration of Type 2 Patients", 
       xlab= "MRI Scan Duration")
  # We can also plot the density of the corresponding normal distribution:
  x_2 <- c(0.2, seq(min(type_2_data$Duration), max(type_2_data$Duration), 0.0001), 1.2)
  lines(x_2, dnorm(x_2, mean=mean_duration_2, sd=sd_duration_2), col="#CC0000", lwd=2.4)
  # --> distribution looks slightly skewed compared to the fitted normal distribution.
  # A Gamma distribution seems like a better fit:
  beta_param <- sd_duration_2^2/mean_duration_2
  alpha_param <- mean_duration_2/beta_param
  lines(x_2, dgamma(x_2, shape=alpha_param, rate=1/beta_param), type="l",
        col="#FF9900", lwd=2)

# Next, we plot the empirical cumulative distribution function of the duration
# with the fitted CDF for both types of patients.
  # the ECDFs:
  ecdf_duration1 <- ecdf(type_1_data$Duration)
  ecdf_duration2 <- ecdf(type_2_data$Duration)
  
  # Plot Type 1 Patients:
  plot(x_1, ecdf_duration1(x_1), type="l", col="deepskyblue4", lwd=2, xlab="", 
       ylab="(Empirical) CDF", main="Scan Duration ECDF of Type 1 Patients")
  # we include the CDF of the fitted normal distribution:
  lines(x_1, pnorm(x_1, mean_duration_1, sd_duration_1), col="#CC0000", lwd=2)
  
  # Plot Type 2 Patients:
  plot(x_2, ecdf_duration2(x_2), type="l", ylab="(Empirical) CDF", xlab="",
       main="Scan Duration ECDF of Type 2 Patients", col="deepskyblue4", lwd=2)
  # we can also include the cdf of the corresponding normal distribution:
  lines(x_2, pnorm(x_2, mean_duration_2, sd_duration_2), type="l", col="#CC0000", lwd=2)
  # and the CDF of the fitted gamma distribution:
  lines(x_2, pgamma(x_2, shape=alpha_param, rate=1/beta_param), type="l",
        col="#FF9900", lwd=2)
  
# Lastly, we can do more appropriate test:
  # A Jarque-Bera Test to test for normality:
  JarqueBera.test(type_1_data$Duration)
  # --> fail to reject the null of normality for duration of type 1 patients.
  JarqueBera.test(type_2_data$Duration)
  # --> reject null of normality for the duration of type 2 patients (mostly due 
  #     to skewnessin the distribution)
  gamma_test(type_2_data$Duration)
  # --> fail to reject the null of a gamma distribution.
  
# We apply a non-parametric bootstrap to measure uncertainty around the distribution
# of the duration of patients and to retrieve information of our distribution. 
# More specifically, for patients of type 1 we are interested in:
#  - Uncertainty around the mean -> apply bootstrap to construct a confidence interval.
#  - Uncertainty around the variance (/sd) -> apply bootstrap to construct a confidence interval.
#  - Quantiles + bootstrapped percentile intervals:
     quants_dur <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975)
     quant.point.1 <- quantile(type_1_data$Duration, quants_dur)

# The bootstrap for statistics around patients of type 1:
  # Step 1: Number of draws
     B <- 100000
  # Step 2: Vectors and Matrices to store the statistics:
  #   - T-statistic for mean:
     T_dur_1 <- rep(NA, B)
  #   - Statistic for variance:
     Q_dur_1 <- rep(NA, B)
  #   - Store values for quantiles
     Quant.matrix.1 <- matrix(NA, nrow=0, ncol=length(quants_dur))
  
  # Step 3: Run the boostrap:
     for(b in 1:B){
       J <- sample.int(n_1, size = n_1, replace = TRUE)
       X.star <- type_1_data$Duration[J]
       mean.star <- mean(X.star)
       sd.star <- sd(X.star)
       
       T.star <- sqrt(n_1)*(mean.star-mean_duration_1)/sd.star
       T_dur_1[b] <- T.star
       
       Q.star <- (n_1-1)*sd.star^2/(sd_duration_1^2)
       Q_dur_1[b] <- Q.star
       
       Quant.matrix.1 <- rbind(Quant.matrix.1, (quantile(X.star, probs=quants_dur)-quant.point.1))
     }
    
  # Step 4: find the critical values
     T.cv.dur.1 <- quantile(T_dur_1, probs=c(0.025, 0.975))
     Q.cv.dur.1 <- quantile(Q_dur_1, probs=c(0.025, 0.975))
     
     quant.dur.cv.1_up <- rep(NA, length(quants_dur))
     quant.dur.cv.1_down <- rep(NA, length(quants_dur))
     for(i in 1:length(quants_dur)){
       quant.dur.cv.1_down[i] <- quantile(Quant.matrix.1[,i], probs=0.025)
       quant.dur.cv.1_up[i] <- quantile(Quant.matrix.1[,i], probs=0.975)
     }
    
  # Step 5: the 95%-confidence-/percentile-intervals:
     CI_mean_dur_1 <- c(mean_duration_1-T.cv.dur.1[2]*sd_duration_1/sqrt(n_1), 
                        mean_duration_1-T.cv.dur.1[1]*sd_duration_1/sqrt(n_1))
     CI_var_dur_1 <- c((n_1-1)*sd_duration_1^2/Q.cv.dur.1[2], 
                       (n_1-1)*sd_duration_1^2/Q.cv.dur.1[1])
     CI_sd_dur_1 <- sqrt(CI_var_dur_1)
     
     # Compare to the theoretical confidence intervals:
     T_CI_mean_1 <- c(mean_duration_1 + qt(0.025, n_1-1)*sd_duration_1/sqrt(n_1), 
                      mean_duration_1 + qt(0.975, n_1-1)*sd_duration_1/sqrt(n_1))
     T_CI_var_1 <- c((n_1-1)*sd_duration_1^2/qchisq(0.975, n_1-1), 
                     (n_1-1)*sd_duration_1^2/qchisq(0.025, n_1-1))
     
    
     PI_quant_dur_1 <- cbind(quant.point.1-quant.dur.cv.1_up, quant.point.1-quant.dur.cv.1_down)

# The bootstrap for patients of type 2:
     # Similarly, we can use bootstrap to find uncertainty around the distribution
     # of the duration of type 2 patients.
     quant.point.2 <- quantile(type_2_data$Duration, quants_dur)
     # Step 1: Number of draws
     B <- 100000
     # Step 2: Vectors and Matrices to store the statistics:
     #   - T-statistic for mean:
     T_dur_2 <- rep(NA, B)
     #   - Statistic for variance:
     Q_dur_2 <- rep(NA, B)
     #   - Store values for quantiles
     Quant.matrix.2 <- matrix(NA, nrow=0, ncol=length(quants_dur))
     
     # Step 3: Run the boostrap:
     for(b in 1:B){
       J <- sample.int(n_2, size = n_2, replace = TRUE)
       X.star <- type_2_data$Duration[J]
       mean.star <- mean(X.star)
       sd.star <- sd(X.star)
       
       T.star <- sqrt(n_2)*(mean.star-mean_duration_2)/sd.star
       T_dur_2[b] <- T.star
       
       Q.star <- (n_2-1)*sd.star^2/(sd_duration_2^2)
       Q_dur_2[b] <- Q.star
       
       Quant.matrix.2 <- rbind(Quant.matrix.2, (quantile(X.star, probs=quants_dur)-quant.point.2))
     }
     
     # Step 4: find the critical values
     T.cv.dur.2 <- quantile(T_dur_2, probs=c(0.025, 0.975))
     Q.cv.dur.2 <- quantile(Q_dur_2, probs=c(0.025, 0.975))
     
     quant.dur.cv.2_up <- rep(NA, length(quants_dur))
     quant.dur.cv.2_down <- rep(NA, length(quants_dur))
     for(i in 1:length(quants_dur)){
       quant.dur.cv.2_down[i] <- quantile(Quant.matrix.2[,i], probs=0.025)
       quant.dur.cv.2_up[i] <- quantile(Quant.matrix.2[,i], probs=0.975)
     }
     
     # Step 5: the 95%-confidence-/percentile-intervals:
     CI_mean_dur_2 <- c(mean_duration_2-T.cv.dur.2[2]*sd_duration_2/sqrt(n_2), 
                        mean_duration_2-T.cv.dur.2[1]*sd_duration_2/sqrt(n_2))
     CI_var_dur_2 <- c((n_2-1)*sd_duration_2^2/Q.cv.dur.2[2], 
                       (n_2-1)*sd_duration_2^2/Q.cv.dur.2[1])
     CI_sd_dur_2 <- sqrt(CI_var_dur_2)
     PI_quant_dur_2 <- cbind(quant.point.2-quant.dur.cv.2_up, quant.point.2-quant.dur.cv.2_down)
     
     
############ Time between calls Dynamics (2.2) ###########################################
# We find the time between consecutive calls for both patient types:
  time_patient1 <- rep(NA, n_1-1)
  time_patient2 <- rep(NA, n_2-1)
 
# For type 1 patients:
   for(i in 2:n_1){
     time <- type_1_data$Time[i]-type_1_data$Time[i-1]
     if(time < 0){
       time_patient1[i-1] <- time+9
     } else {
       time_patient1[i-1] <- time
     }
   }
  
  # For type 2 patients:
  for(i in 2:n_2){
    time <- type_2_data$Time[i]-type_2_data$Time[i-1]
    if(time < 0){
      time_patient2[i-1] <- time+9
    } else {
      time_patient2[i-1] <- time
    }
  }
  
  average_time_1 <- mean(time_patient1)
  sd_time_1 <- sd(time_patient1)
  average_time_2 <- mean(time_patient2)
  sd_time_2 <- sd(time_patient2)
  
# We visualize the distributions by a histogram
  # For patients of type 1:
  hist(time_patient1, col="deepskyblue4", breaks=15, xlab="", 
       main="Time Between Calls of Type 1 Patients", probability = TRUE)
  # Add the line of the pdf of an exponential pdf:
  x_3 <- seq(0, 3.5, 0.05)
  lines(x_3, dexp(x_3, rate=1/average_time_1), col="#9E9D24", lwd=2.4)
  
  # and for patients of type 2:
  hist(time_patient2, col="deepskyblue4", breaks=12, xlab="", 
       main="Time Between Calls of Type 2 Patients", probability = TRUE)
  # Add the line of the pdf of an exponential pdf:
  x_4 <- seq(0, 2, 0.05)
  lines(x_4, dnorm(x_4, mean=average_time_2, sd=sd_time_2), col="#CC0000", lwd=2.4)

# or by visualising the ECDFs:
  # For type 1 patients:
  ecdf_time1 <- ecdf(time_patient1)
  # take a grid:
  x_5 <- seq(min(time_patient1), max(time_patient1), 0.0001)
  plot(x_5, pexp(x_5, rate=1/average_time_1), col="#9E9D24", lwd=2, type="l", 
       ylab="(Empirical) CDF", xlab="",
       main="ECDF of Time between Type 1 Patients' calls")
  lines(x_5, ecdf_time1(x_5), col="deepskyblue4", lwd=2)
  
  # For type 2 patients:
  ecdf_time2 <- ecdf(time_patient2)
  x_6 <- seq(min(time_patient2), max(time_patient2), 0.0001)
  plot(x_6, pnorm(x_6, average_time_2, sd_time_2), col="#CC0000", lwd=2, type="l", 
       ylab="(Empirical) CDF", xlab="",
       main="ECDF of Time between Type 2 Patients' calls")
  lines(x_6, ecdf_time2(x_6), col="deepskyblue4", lwd=2)
  
# Formal Test:
  # Test for exponential distribution for time between two calls for type 1 patients:
  # --> Kolmogorov–Smirnov test
  ks.test(time_patient1, pexp(length(time_patient1), rate=1/average_time_1))
  # we fail to reject the null that time_patient1 follows an exponential distribution  
  
  # Test for normality for the distribution of the time between two calls for type 2 patients:
  # --> Jarque-Bera Test
  JarqueBera.test(time_patient2)
  # --> fails to reject the null of normality.
  
  # We apply bootstrap for uncertainty surrounding:
  # - Patients of type 1:
  #   + Uncertainty around mean, using Confidence Intervals
  #   + Quantiles + uncertainty (Percentile Intervals) at the following percentages:
  quants <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95, 0.975)
  sample_quants1 <- quantile(time_patient1, probs=quants)
  sample_quants2 <- quantile(time_patient2, probs=quants)
  # - Patients of type 2:
  #   + uncertainty around mean, using Confidence Intervals
  #   + uncertainty around variance, using Confidence Intervals
  #   + Quantiles + percentile intervals for the same percentages.
  
  # Step 1: Choose number of sample draws
  B <- 10000
  # Step 2: Define arrays to store statistics:
  T_time_1 <- rep(NA, B)
  Quantile_matrix.1 <- matrix(data=NA, nrow=0, ncol=length(quants))
  
  T_time_2 <- rep(NA, B)
  Q_time_2 <- rep(NA, B)
  Quantile_matrix.2 <- matrix(data=NA, nrow=0, ncol=length(quants))
  
  # Step 3: Run the bootstrap
  for(b in 1:B){
    #X.star.1 <- rexp((n_1-1), rate=1/average_time_1)
    J.1 <- sample(1:(n_1-1), replace=TRUE)
    X.star.1 <- time_patient1[J.1]
    mean.star.1 <- mean(X.star.1)
    sd.star.1 <- sd(X.star.1)
    T.star.1 <- sqrt(n_1-1)*(mean.star.1-average_time_1)/sd.star.1
    T_time_1[b] <- T.star.1
    Quantile_matrix.1 <- rbind(Quantile_matrix.1, (quantile(X.star.1, probs=quants)-sample_quants1))
    
    J.2 <- sample(1:(n_2-1), replace=TRUE)
    X.star.2 <- time_patient2[J.2]
    mean.star.2 <- mean(X.star.2)
    sd.star.2 <- sd(X.star.2)
    T.star.2 <- sqrt(n_2-1)*(mean.star.2-average_time_2)/sd.star.2
    T_time_2[b] <- T.star.2
    
    Q.star.2 <- (n_2-2)*sd.star.2^2/(sd_time_2^2)
    Q_time_2[b] <- Q.star.2
    Quantile_matrix.2 <- rbind(Quantile_matrix.2, (quantile(X.star.2, probs=quants)-sample_quants2))
  }
  
  # Step 4: find critical values for the T-statistics and Q-statistic:
  cvs_T_1 <- quantile(T_time_1, probs=c(0.025, 0.975))
  cvs_T_2 <- quantile(T_time_2, probs=c(0.025, 0.975))
  cvs_Q_2 <- quantile(Q_time_2, probs=c(0.025, 0.975))
  
  quant.time.cv.1_up <- rep(NA, length(quants))
  quant.time.cv.1_down <- rep(NA, length(quants))
  for(i in 1:length(quants)){
    quant.time.cv.1_down[i] <- quantile(Quantile_matrix.1[,i], probs=0.025)
    quant.time.cv.1_up[i] <- quantile(Quantile_matrix.1[,i], probs=0.975)
  }
  
  quant.time.cv.2_up <- rep(NA, length(quants))
  quant.time.cv.2_down <- rep(NA, length(quants))
  for(i in 1:length(quants)){
    quant.time.cv.2_down[i] <- quantile(Quantile_matrix.2[,i], probs=0.025)
    quant.time.cv.2_up[i] <- quantile(Quantile_matrix.2[,i], probs=0.975)
  }
  
  # Step 5: Construct Confidence and Percentile Intervals:
  CI_mean_time_1 <- c(average_time_1-cvs_T_1[2]*sd_time_1/sqrt(length(time_patient1)), 
                      average_time_1-cvs_T_1[1]*sd_time_1/sqrt(length(time_patient1)))
  CI_mean_time_2 <- c(average_time_2-cvs_T_2[2]*sd_time_2/sqrt(length(time_patient2)), 
                      average_time_2-cvs_T_2[1]*sd_time_2/sqrt(length(time_patient2)))
  CI_var_time_2 <-  c((length(time_patient2)-1)*sd_time_2^2/cvs_Q_2[2], 
                      (length(time_patient2)-1)*sd_time_2^2/cvs_Q_2[1])
  CI_sd_time_2 <- sqrt(CI_var_time_2)
  
  PI_quant_time_2 <- cbind(sample_quants2-quant.time.cv.2_up, 
                           sample_quants2-quant.time.cv.2_down)
  PI_quant_time_1 <- cbind(sample_quants1-quant.time.cv.1_up, 
                           sample_quants1-quant.time.cv.1_down)
  
  











