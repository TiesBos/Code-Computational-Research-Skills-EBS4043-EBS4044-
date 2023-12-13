################################################################################
# Statistical Analysis Computational Research Skills
# Authors: Ties Bos
################################################################################
# packages used:
library("dplyr")
library("stats")
library("extraDistr")
library("tsoutliers")

############ Import the Dataset ################################################
scan.data <- read.csv("ScanRecords.csv", header=TRUE)

# We check whether the data has been imported nicely:
head(scan.data)
summary(scan.data)

############ Type 1 Patients ###################################################
  # We split take the subset of type I patients:
  type_1_ind <- scan.data$PatientType == "Type 1"
  # We select the rows of the dataset of type 1 patients:
  type_1_data <- scan.data[type_1_ind, ]
  # No of type 1 patients:
  n_1 <- nrow(type_1_data)
  
  #---------- Duration of Type 1 Patients --------------------------------------
    # As it is given that the duration of type 1 patients are iid and normally
    # distributed, our main concern regards uncertainty surrounding the mean and  
    # variance of this distribution.
  
    # Sample mean:
    mean_duration_1 <- mean(type_1_data$Duration)
    # Sample sd:
    sd_duration_1 <- sd(type_1_data$Duration)
    
    # Plotting the histogram of the durations:
    hist(type_1_data$Duration, col="deepskyblue4", prob=TRUE, 
         main="MRI Scan Duration of Type 1 Patients", xlab="MRI Scan Duration", breaks=20)
    # Choose a grid:
    x <- seq(min(type_1_data$Duration), max(type_1_data$Duration), 0.01)
    lines(x, dnorm(x, mean=mean_duration_1, sd=sd_duration_1), col="#CC0000", lwd=2.4)
    
    # Plotting the ECDF:
    ecdf_duration1 <- ecdf(type_1_data$Duration)
    plot(x, ecdf_duration1(x), type="l", col="deepskyblue4", lwd=2, xlab="", 
         ylab="(Empirical) CDF", main="Scan Duration ECDF of Type 1 Patients")
    lines(x, pnorm(x, mean_duration_1, sd_duration_1), col="#CC0000", lwd=2)
    
    # To obtain valid confidence intervals surrounding our parameters, we need 
    # the correct critical values. 
    
      #   [UNCERTAINTY SURROUNDING THE MEAN]
      # Option 1: We use parametric bootstrapping to do so:
      # 1) The number of sample draws:
      B <- 100000
      # 2) We use the following vector to store the T-statistics:
      T_boot_1 <- rep(NA, B)
      # 3) Run the bootstrap:
      for(b in 1:B){
        X.star <- rnorm(n_1, mean_duration_1, sd=sd_duration_1)
        mean_star <- mean(X.star)
        sd_star <- sd(X.star)
        T_star <- sqrt(n_1)*(mean_star-mean_duration_1)/sd_star
        T_boot_1[b] <- T_star
      }
      # 4) Order the test-statistics;
      ordered_T_1 <- T_boot_1[order(T_boot_1)]
      # 5) Find the quantiles:
      alpha_1 <- 0.05
      cv_lower <- ordered_T_1[(alpha_1/2)*B]
      cv_up <- ordered_T_1[(1-alpha_1/2)*B]
      # 6) the confidence intervals:
      CI_mean_1_low <- mean_duration_1 - cv_up*sd_duration_1/sqrt(n_1)
      CI_mean_1_up <- mean_duration_1 - cv_lower*sd_duration_1/sqrt(n_1)
      boot_CI_mean_1 <- c(CI_mean_1_low, CI_mean_1_up)
      
      # Option 2: Recognize that as the duration are iid normally distributed,
      # the T-statistic follows a T-distribution with n_1-1 degrees of freedom. 
      # Then,
      T_CI_mean_1 <- c(mean_duration_1 + qt(alpha_1/2, n_1-1)*sd_duration_1/sqrt(n_1), 
                       mean_duration_1 + qt(1-alpha_1/2, n_1-1)*sd_duration_1/sqrt(n_1))
    
      #   [UNCERTAINTY SURROUNDING THE VARIANCE]
      # Option 1: We use parametric bootstrapping to do so:
      # 1) The number of sample draws:
      B <- 100000
      # 2) We use the following vector to store the statistics for tests around variance:
      Q_boot_1 <- rep(NA, B)
      # 3) Run the bootstrap:
      for(b in 1:B){
        X.star <- rnorm(n_1, mean_duration_1, sd=sd_duration_1)
        sd_star <- sd(X.star)
        Q_star <- (n_1-1)*sd_star^2/(sd_duration_1^2)
        Q_boot_1[b] <- Q_star
      }
      # 4) Order the test-statistics;
      ordered_Q_1 <- Q_boot_1[order(Q_boot_1)]
      # 5) Find the quantiles:
      alpha_1 <- 0.05
      Q_cv_lower <- ordered_Q_1[(alpha_1/2)*B]
      Q_cv_up <- ordered_Q_1[(1-alpha_1/2)*B]
      # 6) the confidence intervals:
      CI_var_1_low <- (n_1-1)*sd_duration_1^2/Q_cv_up
      CI_var_1_up <- (n_1-1)*sd_duration_1^2/Q_cv_lower
      boot_CI_var_1 <- c(CI_var_1_low, CI_var_1_up)
      # or, similarly, for the standard deviation:
      boot_CI_sd_1 <- sqrt(boot_CI_var_1)
      
      # Option 2: Recognize that as the duration are iid normally distributed,
      # the test-statistic for the test surrounding the variance follows a chi-
      # squared with n_1-1 degrees of freedom. Then,
      Q_CI_var_1 <- c((n_1-1)*sd_duration_1^2/qchisq(1-alpha_1/2, n_1-1), 
                       (n_1-1)*sd_duration_1^2/qchisq(alpha_1/2, n_1-1))
      # or for the standard deviation:
      Q_CI_sd_1 <- sqrt(Q_CI_var_1)

  #---------- Number of Type 1 Patients ----------------------------------------
    # It is given that the number of daily patients follows a poisson(lambda) 
    # distribution. We are interested in the uncertainty around the lambda, which
    # equals the expected value of the distribution, and therefore around the mean
    # of the distribution. 
    # --> we can find confidence-intervals of lambda with parametric bootstrapping
    #     using the T-statistic.
      
      
    # Before we can bootstrap, we aggregate the daily number of type 1 patients:
      dates <- scan.data$Date
      # Select only the dates that are not duplicated
      unique_dates <- dates[!duplicated(dates)]
      n_dates <- length(unique_dates)
      daily_1 <- rep(0, n_dates)
      for(i in 1:length(unique_dates)){
        date_i <- unique_dates[i]
        daily_1[i] <- sum(type_1_data$Date==date_i)
      }
      daily_type_1 <- data.frame(unique_dates, daily_1)
      colnames(daily_type_1) <- c("Date", "# Type 1 Patients")
    
      # [ Parametric Bootstrapping]
    # Sample Mean number of daily patients of type 1:
      mean_daily_1 <- mean(daily_type_1$`# Type 1 Patients`)
    # Sample SD of number of daily patients of type 1:
      sd_daily_1 <- sd(daily_type_1$`# Type 1 Patients`)
    # Next, we use parametric with the parameter estimate being the sample mean:
    # 1) set the number of sample draws:
      B <- 10000
    # 2) Define an array to store the T-statistics:
      T_daily_1_param <- rep(NA,B)
    # 3) Run the parametric bootstrap:
      for(b in 1:B){
        X.star <- rpois(n_dates, mean_daily_1)
        sd.star <- sd(X.star)
        mean.star <- mean(X.star)
        T.star <- sqrt(n_dates)*(mean.star-mean_daily_1)/sd.star
        T_daily_1_param[b] <- T.star
      }
    # 4) We order the T-statistics
      ordered_T_1_param <- T_daily_1_param[order(T_daily_1_param)]
    # 5) Find the quantiles of this ordered vector, i.e. the critical values
      alpha_1_param <- 0.05
      daily_1_cv_up <- ordered_T_1_param[B*(1-alpha_1_param/2)]
      daily_1_cv_low <- ordered_T_1_param[B*alpha_1_param/2]
    # 6) We construct the bootstrapped confidence intervals
      CI_lambda_param_low <- mean_daily_1-daily_1_cv_up*sd_daily_1/sqrt(n_dates)
      CI_lambda_param_up <- mean_daily_1-daily_1_cv_low*sd_daily_1/sqrt(n_dates)
      CI_lambda_param <- c(CI_lambda_param_low, CI_lambda_param_up)
      
    #   [Non-Parametric Approach]
    # Next, we can also use a non-parametric approach:
    # 1) set the number of sample draws:
    B <- 10000
    # 2) Define an array to store the T-statistics:
    T_daily_1_non_param <- rep(NA,B)
    # 3) Run the parametric bootstrap:
    for(b in 1:B){
      J <- sample.int(n_dates, size = n_dates, replace = TRUE)
      X.star <- daily_type_1$`# Type 1 Patients`[J]
      sd.star <- sd(X.star)
      mean.star <- mean(X.star)
      T.star <- sqrt(n_dates)*(mean.star-mean_daily_1)/sd.star
      T_daily_1_non_param[b] <- T.star
    }
    # 4) We order the T-statistics
    ordered_T_1_non_param <- T_daily_1_non_param[order(T_daily_1_non_param)]
    # 5) Find the quantiles of this ordered vector, i.e. the critical values
    alpha_1_np <- 0.05
    daily_1_cv_up_np <- ordered_T_1_non_param[B*(1-alpha_1_np/2)]
    daily_1_cv_low_np <- ordered_T_1_non_param[B*alpha_1_np/2]
    # 6) We construct the bootstrapped confidence intervals
    CI_lambda_np_low <- mean_daily_1-daily_1_cv_up_np*sd_daily_1/sqrt(n_dates)
    CI_lambda_np_up <- mean_daily_1-daily_1_cv_low_np*sd_daily_1/sqrt(n_dates)
    CI_lambda_np <- c(CI_lambda_np_low, CI_lambda_np_up)  
      
  
############ Type 2 Patients ###################################################
  # We split take the subset of type 2 patients:
  type_2_ind <- scan.data$PatientType == "Type 2"
  # We select the rows of the dataset of type 2 patients:
  type_2_data <- scan.data[type_2_ind, ]
  # Number of of type 2 patients:
  n_2 <- nrow(type_2_data)
  
  # -- Duration of the scans ---------------------------------------------------
  # The average duration of the scans:
  mean_duration_2 <- mean(type_2_data$Duration)
  # The average standard deviation of the duration of the scan:
  sd_duration_2 <- sd(type_2_data$Duration)
  
  # To gain some insights as to how the duration is distributed for type 2 patients
  # we start by constructing the EDF function:
  ecdf_duration2 <- ecdf(type_2_data$Duration)
  # take a grid:
  x <- seq(min(type_2_data$Duration), max(type_2_data$Duration), 0.0001)
  plot(x, ecdf_duration2(x), type="l", ylab="(Empirical) CDF", xlab="",
       main="Scan Duration ECDF of Type 2 Patients", col="deepskyblue4", lwd=2)
  
  # we can also include the cdf of the corresponding normal distribution:
  lines(x, pnorm(x, mean_duration_2, sd_duration_2), type="l", col="#CC0000", lwd=2)
  # --> normal distribution seems like an appropriate fit. 
  lines(x, pgamma(x, shape=alpha_param, rate=1/beta_param), type="l",
        col="#FF9900", lwd=2)
  
  # Looking at the histogram:
  hist(type_2_data$Duration, breaks=10, probability=TRUE, col="deepskyblue4", 
       main="MRI Scan Duration of Type 2 Patients", 
       xlab= "MRI Scan Duration")
  # We can also plot the density of the corresponding normal distribution:
  lines(x, dnorm(x, mean=mean_duration_2, sd=sd_duration_2), col="#CC0000", lwd=2.4)
  # --> distribution looks slightly skewed compared to the fitted normal distribution.
  # A Gamma distribution seems like a better fit:
  beta_param <- sd_duration_2^2/mean_duration_2
  alpha_param <- mean_duration_2/beta_param
  lines(c(0.2, x, 1.2), dgamma(c(0.2, x, 1.2), shape=alpha_param, rate=1/beta_param), type="l",
        col="#FF9900", lwd=2)
  
  # We use non-parametric bootstrapping to find out the following specifications of this distribution:
  #   - Mean of duration + uncertainty (i.e. Confidence Intervals)
  #   - SD of durations + uncertainty (i.e. Confidence Intervals)
  #   - probability exceeding a certain thresholds (using equal-tailed percentile intervals):
  Tr1 <- 1 #(one hour)
  Tr2 <- 0.75 #(45 minutes)
  Tr3 <- 50/60 #(50 minutes)
  Tr4 <- 0.5
  #     We use the estimated probability based on the ECDF:
  prob_Tr1_hat <- ecdf_duration2(Tr1)
  prob_Tr2_hat <- ecdf_duration2(Tr2)
  prob_Tr3_hat <- ecdf_duration2(Tr3)
  prob_Tr4_hat <- ecdf_duration2(Tr4)
  
  # Non-Parametric Bootstrap:
  # 1) Set the number of sample draws B:
  B <- 100000
  # 2) Vectors to store statistics:
  # Mean:
  T_boot_2 <- rep(NA, B)
  # SD:
  Q_boot_2 <- rep(NA, B)
  # Prob crossing a threshold
  Tr1_boot_2 <- rep(NA, B)
  Tr2_boot_2 <- rep(NA, B)
  Tr3_boot_2 <- rep(NA, B)
  Tr4_boot_2 <- rep(NA, B)
  
  # 3) Run the bootstrap:
  for(b in 1:B){
    J <- sample.int(n_2, size = n_2, replace = TRUE)
    X.star <- type_2_data$Duration[J]
    mean.star <- mean(X.star)
    sd.star <- sd(X.star)
    T.star <- sqrt(n_2)*(mean.star-mean_duration_2)/sd.star
    T_boot_2[b] <- T.star
    
    Q_boot_2[b] <- (n_2-1)*sd.star^2/(sd_duration_2^2)
    
    ecdf_star <- ecdf(X.star)
    Tr1_boot_2[b] <- prob_Tr1_hat-ecdf_star(Tr1)
    Tr2_boot_2[b] <- prob_Tr2_hat-ecdf_star(Tr2)
    Tr3_boot_2[b] <- prob_Tr3_hat-ecdf_star(Tr3)
    Tr4_boot_2[b] <- prob_Tr4_hat-ecdf_star(Tr4)
  }

  # 4) Order the statistics:
  ordered_Q_stat <- Q_boot_2[order(Q_boot_2)]
  ordered_T_2 <- T_boot_2[order(T_boot_2)]
  Tr1_boot_2_ord <- Tr1_boot_2[order(Tr1_boot_2)]
  Tr2_boot_2_ord <- Tr2_boot_2[order(Tr2_boot_2)]
  Tr3_boot_2_ord <- Tr3_boot_2[order(Tr3_boot_2)]
  Tr4_boot_2_ord <- Tr4_boot_2[order(Tr4_boot_2)]
  
  # 5) constructing the intervals:
  # The equal-tailed 95%-percentile interval surrounding the probability of exceeding Tr:
  alpha_3 <- 0.05
  PI_prob_Tr1 <- c(prob_Tr1_hat-quantile(Tr1_boot_2_ord, probs=c(0.025, 0.975))[2],
                   prob_Tr1_hat-quantile(Tr1_boot_2_ord, probs=c(0.025, 0.975))[1])
  PI_prob_Tr2 <- c(prob_Tr2_hat-quantile(Tr2_boot_2_ord, probs=c(0.025, 0.975))[2],
                   prob_Tr2_hat-quantile(Tr2_boot_2_ord, probs=c(0.025, 0.975))[1])
  PI_prob_Tr3 <- c(prob_Tr3_hat-quantile(Tr3_boot_2_ord, probs=c(0.025, 0.975))[2],
                   prob_Tr3_hat-quantile(Tr3_boot_2_ord, probs=c(0.025, 0.975))[1])
  PI_prob_Tr4 <- c(prob_Tr4_hat-quantile(Tr4_boot_2_ord, probs=c(0.025, 0.975))[2],
                   prob_Tr4_hat-quantile(Tr4_boot_2_ord, probs=c(0.025, 0.975))[1])
  
  # The equal-tailed confidence interval surrounding the mean duration:
  #cv_T_dur_2_up <- ordered_T_2[(1-alpha_3/2)*B]
  #cv_T_dur_2_down <- ordered_T_2[(alpha_3/2)*B]
  cv_T_dur_2 <- quantile(ordered_T_2, probs=c(0.025, 0.975))
  # --> the critical values being suspiciously close to those if the durations follow
  #     a normal distribution would be another argument to approximate using a 
  #     normal distribution.
  CI_mean_dur_2 <- c(mean_duration_2-cv_T_dur_2[2]*sd_duration_2/sqrt(n_2), 
                     mean_duration_2-cv_T_dur_2[1]*sd_duration_2/sqrt(n_2))
  
  # The confidence interval surrounding the variance:
  cv_Q_dur_2 <- quantile(ordered_Q_stat, probs=c(0.025, 0.975))
  CI_Q_dur_2 <- c((n_2-1)*sd_duration_2^2/cv_Q_dur_2[2], 
                  (n_2-1)*sd_duration_2^2/cv_Q_dur_2[1])
  CI_SD_dur_2 <- sqrt(CI_Q_dur_2)
  
  # -- Daily number of scans ---------------------------------------------------
  # Again, we start by aggregating the daily demand for type 2 patients:
  daily_2 <- rep(0, n_dates)
  for(i in 1:n_dates){
    date_i <- unique_dates[i]
    daily_2[i] <- sum(type_2_data$Date==date_i)
  }
  daily_type_2 <- data.frame(unique_dates, daily_2)
  colnames(daily_type_2) <- c("Date", "# Type 2 Patients")
  
  mean_daily_2 <- mean(daily_type_2$`# Type 2 Patients`)
  sd_daily_2 <- sd(daily_type_2$`# Type 2 Patients`)
  
  # We are interested in:
  #   - uncertainty surrounding the mean
  #   - probability of crossing a threshold Tr_daily
  # --> Use non-parametric bootstrap to construct percentile intervals:
  
  # 1) Set the number of Bootstrap sample draws:
    B <- 100000
  # 2) construct arrays to store the test-statistics:
    # Uncertainty around mean (T-statistics):
    T_daily_2 <- rep(NA, B)
    # Uncertainty around threshold (Call statistic Q)
    Tr <- 10
    prob_Tr_hat <- rep(NA, B)
    
  # 3) Run Bootstrap
    for(b in 1:B){
      J <- sample(1:n_dates, replace=TRUE)
      X.star <- daily_2[J]
      mean.star <- mean(X.star)
      sd.star <- sd(X.star)
      
      T.star <- sqrt(n_dates)*(mean.star-mean_daily_2)/sd.star
      T_daily_2[b] <- T.star
      
      ecdf_star <- ecdf(X.star)
      prob_Tr_hat[b] <- ecdf_star(Tr)
    }
  
  # 4) order the statistics:
    prob_Tr_hat <- prob_Tr_hat[order(prob_Tr_hat)]
    T_daily_2 <- T_daily_2[order(T_daily_2)]
    
  # 5) Find the quantiles:
    cv_T_daily_2 <- quantile(T_daily_2, probs=c(0.025, 0.975))
    cv_prob_Tr_hat <- quantile(prob_Tr_hat,probs=c(0.025, 0.975))
    
  # 6) Create percentile/confidence intervals:
    # Confidence interval around the mean:
    mean_daily_2_PI <- c(mean_daily_2 - cv_T_daily_2[2]*sd_daily_2/sqrt(n_dates), 
                         mean_daily_2- cv_T_daily_2[1]*sd_daily_2/sqrt(n_dates))
    
    # Percentile interval around the probability of exceeding the threshold Tr:
    prob_Tr_daily_2_PI <- cv_prob_Tr_hat
    
    # Approximating the distribution
    ecdf_daily_2 <- ecdf(daily_type_2$`# Type 2 Patients`)
    x <- seq(min(daily_type_2$`# Type 2 Patients`), max(daily_type_2$`# Type 2 Patients`), 1)
    plot(ecdf_daily_2)  
    hist(daily_type_2$`# Type 2 Patients`, prob=TRUE)
    lines(x, pgeom(x, 0.02), col="red")
    
#============= TIME BETWEEN PATIENTS ===========================================
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
  
  average_time_1 <- mean(time_patient1)
  sd_time_1 <- sd(time_patient1)
  
  # We plot the histogram of the time between arrivals:
  hist(time_patient1, col="deepskyblue4", breaks=15, xlab="", 
       main="Time Between Calls of Type 1 Patients", probability = TRUE)
  # Add the line of the pdf of an exponential pdf:
  x <- seq(0, 3.5, 0.05)
  lines(x, dexp(x, rate=1/average_time_1), col="#9E9D24", lwd=2.4)
  
  # Similarly, we plot the ECDF:
  ecdf_time1 <- ecdf(time_patient1)
  # take a grid:
  x <- seq(min(time_patient1), max(time_patient1), 0.0001)
  plot(x, ecdf_time1(x), type="l", ylab="(Empirical) CDF", xlab="",
       main="ECDF of Time between Type 1 Patients' calls", col="deepskyblue4", 
       lwd=2)
  # and the CDF of the fitted exponential distribution:
  lines(x, pexp(x, rate=1/average_time_1), col="#9E9D24", lwd=2)
  
  # Alternatively:
  x <- seq(min(time_patient1), max(time_patient1), 0.0001)
  plot(x, pexp(x, rate=1/average_time_1), col="#9E9D24", lwd=2, type="l", 
       ylab="(Empirical) CDF", xlab="",
       main="ECDF of Time between Type 1 Patients' calls")
  lines(x, ecdf_time1(x), col="deepskyblue4", lwd=2)
  
  
  # For type 2 patients:
  for(i in 2:n_2){
    time <- type_2_data$Time[i]-type_2_data$Time[i-1]
    if(time < 0){
      time_patient2[i-1] <- time+9
    } else {
      time_patient2[i-1] <- time
    }
  }
  
  average_time_2 <- mean(time_patient2)
  sd_time_2 <- sd(time_patient2)
    
  # We plot the histogram of the time between arrivals:
  hist(time_patient2, col="deepskyblue4", breaks=12, xlab="", 
       main="Time Between Calls of Type 2 Patients", probability = TRUE)
  # Add the line of the pdf of an exponential pdf:
  x <- seq(0, 2, 0.05)
  lines(x, dnorm(x, mean=average_time_2, sd=sd_time_2), col="#CC0000", lwd=2.4)
  # Formal test for normality:
  JarqueBera.test(time_patient2)
  # --> fails to reject the null of normality.
  
  # Similarly, we plot the ECDF:
  ecdf_time2 <- ecdf(time_patient2)
  # take a grid:
  x <- seq(min(time_patient2), max(time_patient2), 0.0001)
  plot(x, ecdf_time2(x), type="l", ylab="(Empirical) CDF", xlab="",
       main="ECDF of Time between Type 2 Patients' calls", col="deepskyblue4", 
       lwd=2)
  # and the CDF of the fitted exponential distribution:
  lines(x, pnorm(x, average_time_2, sd_time_2), col="#CC0000", lwd=2) 
  # --> Fits very well.  
  
  # Alternatively:
  x <- seq(min(time_patient2), max(time_patient2), 0.0001)
  plot(x, pnorm(x, average_time_2, sd_time_2), col="#CC0000", lwd=2, type="l", 
       ylab="(Empirical) CDF", xlab="",
       main="ECDF of Time between Type 2 Patients' calls")
  lines(x, ecdf_time2(x), col="deepskyblue4", lwd=2)
  
  
  
      
  # We apply bootstrap for uncertainty surrounding:
  # - Patients of type 1:
  #   + Uncertainty around mean, using Confidence Intervals
  #   + Quantiles + uncertainty (Percentile Intervals) at the following percentages:
        quants <- c(0.025, 0.05, 0.5, 0.75, 0.8, 0.85, 0.90, 0.95, 0.975)
        sample_quants1 <- quantile(time_patient1, probs=quants)
        sample_quants2 <- quantile(time_patient2, probs=quants)
  # - Patients of type 2:
  #   + uncertainty around mean, using Confidence Intervals
  #   + uncertainty around variance, using Confidence Intervals
  #   + Quantiles + percentile intervals for the same percentages.
  
  # Step 1: Choose number of sample draws
  B <- 100000
  # Step 2: Define arrays to store statistics:
  T_time_1 <- rep(NA, B)
  Quantile_matrix.1 <- matrix(data=NA, nrow=0, ncol=length(quants))
  
  T_time_2 <- rep(NA, B)
  Q_time_2 <- rep(NA, B)
  Quantile_matrix.2 <- matrix(data=NA, nrow=0, ncol=length(quants))
  
  # Step 3: Run the bootstrap
  for(b in 1:B){
    X.star.1 <- rexp((n_1-1), rate=1/average_time_1)
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
  
  cv_quant1_1 <- quantile(Quantile_matrix.1[,1], probs=c(0.025, 0.975))
  cv_quant2_1 <- quantile(Quantile_matrix.1[,2], probs=c(0.025, 0.975))
  cv_quant3_1 <- quantile(Quantile_matrix.1[,3], probs=c(0.025, 0.975))
  cv_quant4_1 <- quantile(Quantile_matrix.1[,4], probs=c(0.025, 0.975))
  cv_quant5_1 <- quantile(Quantile_matrix.1[,5], probs=c(0.025, 0.975))
  cv_quant6_1 <- quantile(Quantile_matrix.1[,6], probs=c(0.025, 0.975))
  cv_quant7_1 <- quantile(Quantile_matrix.1[,7], probs=c(0.025, 0.975))
  cv_quant8_1 <- quantile(Quantile_matrix.1[,8], probs=c(0.025, 0.975))
  cv_quant9_1 <- quantile(Quantile_matrix.1[,9], probs=c(0.025, 0.975))
  
  cv_quant1_2 <- quantile(Quantile_matrix.2[,1], probs=c(0.025, 0.975))
  cv_quant2_2 <- quantile(Quantile_matrix.2[,2], probs=c(0.025, 0.975))
  cv_quant3_2 <- quantile(Quantile_matrix.2[,3], probs=c(0.025, 0.975))
  cv_quant4_2 <- quantile(Quantile_matrix.2[,4], probs=c(0.025, 0.975))
  cv_quant5_2 <- quantile(Quantile_matrix.2[,5], probs=c(0.025, 0.975))
  cv_quant6_2 <- quantile(Quantile_matrix.2[,6], probs=c(0.025, 0.975))
  cv_quant7_2 <- quantile(Quantile_matrix.2[,7], probs=c(0.025, 0.975))
  cv_quant8_2 <- quantile(Quantile_matrix.2[,8], probs=c(0.025, 0.975))
  cv_quant9_2 <- quantile(Quantile_matrix.2[,9], probs=c(0.025, 0.975))
  
  # Step 5: construct confidence and percentile intervals
  CI_param_1 <- c(average_time_1 - cvs_T_1[2]*sd_time_1/sqrt(n_1-1), 
                  average_time_1 - cvs_T_1[1]*sd_time_1/sqrt(n_1-1))
  
  CI_mean_time_2 <- c(average_time_2 - cvs_T_2[2]*sd_time_2/sqrt(n_2-1), 
                      average_time_2 - cvs_T_2[1]*sd_time_2/sqrt(n_2-1))
  CI_var_time_2 <- c((n_2-2)*sd_time_2^2/cvs_Q_2[2], (n_2-2)*sd_time_2^2/cvs_Q_2[1])
  CI_sd_time_2 <- sqrt(CI_var_time_2)
  
  PI_quant1_1 <- c(sample_quants1[1] - cv_quant1_1[2], 
                   sample_quants1[1] - cv_quant1_1[1])
  PI_quant2_1 <- c(sample_quants1[2] - cv_quant2_1[2], 
                   sample_quants1[2] - cv_quant2_1[1])
  PI_quant3_1 <- c(sample_quants1[3] - cv_quant3_1[2], 
                   sample_quants1[3] - cv_quant3_1[1])
  PI_quant4_1 <- c(sample_quants1[4] - cv_quant4_1[2], 
                   sample_quants1[4] - cv_quant4_1[1])
  PI_quant5_1 <- c(sample_quants1[5] - cv_quant5_1[2], 
                   sample_quants1[5] - cv_quant5_1[1])
  PI_quant6_1 <- c(sample_quants1[6] - cv_quant6_1[2], 
                   sample_quants1[6] - cv_quant6_1[1])
  PI_quant7_1 <- c(sample_quants1[7] - cv_quant7_1[2], 
                   sample_quants1[7] - cv_quant7_1[1])
  PI_quant8_1 <- c(sample_quants1[8] - cv_quant8_1[2], 
                   sample_quants1[8] - cv_quant8_1[1])
  PI_quant9_1 <- c(sample_quants1[9] - cv_quant9_1[2], 
                   sample_quants1[9] - cv_quant9_1[1])
  
  PI_quant1_2 <- c(sample_quants2[1] - cv_quant1_2[2], 
                   sample_quants2[1] - cv_quant1_2[1])
  PI_quant2_2 <- c(sample_quants2[2] - cv_quant2_2[2], 
                   sample_quants2[2] - cv_quant2_2[1])
  PI_quant3_2 <- c(sample_quants2[3] - cv_quant3_2[2], 
                   sample_quants2[3] - cv_quant3_2[1])
  PI_quant4_2 <- c(sample_quants2[4] - cv_quant4_2[2], 
                   sample_quants2[4] - cv_quant4_2[1])
  PI_quant5_2 <- c(sample_quants2[5] - cv_quant5_2[2], 
                   sample_quants2[5] - cv_quant5_2[1])
  PI_quant6_2 <- c(sample_quants2[6] - cv_quant6_2[2], 
                   sample_quants2[6] - cv_quant6_2[1])
  PI_quant7_2 <- c(sample_quants2[7] - cv_quant7_2[2], 
                   sample_quants2[7] - cv_quant7_2[1])
  PI_quant8_2 <- c(sample_quants2[8] - cv_quant8_2[2], 
                   sample_quants2[8] - cv_quant8_2[1])
  PI_quant9_2 <- c(sample_quants2[9] - cv_quant9_2[2], 
                   sample_quants2[9] - cv_quant9_2[1])

  
  
  
  quantile(Quantile_matrix.1[,1], probs=c(0.975, 0.025))
  
  
  
  
  
  