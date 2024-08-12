#Pooling functions and other functions
#Author: Mary Ann Binuya
#Last updated: December 21, 2023

#References:
  #Pooling: https://bookdown.org/mwheymans/bookmi/rubins-rules.html; https://bookdown.org/mwheymans/bookmi/pooling-methods-for-categorical-variables.html
  #Transformations for pooling:
    ##Marshall, et al, 2009 (https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-9-57)
    ##Debray et al, 2017 (https://www.bmj.com/content/356/bmj.i6460)
  #Calibration curves: https://cran.r-project.org/web/packages/CalibrationCurves/CalibrationCurves.pdf
  #Log-log transformation for predicted probabilities in calibration curves: Austin, et al., 2020

#1. Pool estimates and produce 95% confidence intervals:
  #est <- array of (performance) estimates
  #se <- array of standard errors
  pool_estimates <- function(est, se, logit_trans = FALSE, log_trans = FALSE, exp_trans = FALSE){
    RR_se <- function(est, se){
        m <- length(est)
        w_var <- mean(se^2) # within variance
        b_var <- var(est) # between variance
        t_var <- w_var + (1 + (1/m)) * b_var # total variance
        se_total <- sqrt(t_var) # total se
        r <- (1 + 1 / m) * (b_var / w_var)
        df <- (m - 1) * (1 + (1/r))^2 # degrees of freedom
        t <- qt(0.975, df) # inverse cdf of t dist (~N dist for large df)
        res <- c(se_total, t)
        return(res)
      }
      
      est <- unlist(est)
      se <- unlist(se)
      
      if(logit_trans){ #if logit transformation requested
        est_log <- log(est/(1-est)) 
        se_log <- se / (est * (1-est)) #for formula, see supplementary page 9 of Debray et al, 2017
        se_total <- RR_se(est_log, se_log) #RR pooling after logit transformation
        
        inv.est <- round(exp(mean(est_log))/(1+exp(mean(est_log))), 2) #back-transform
        inv.est.u <- round(exp(mean(est_log) + (se_total[2]*se_total[1])) /
          (1 + exp(mean(est_log) + (se_total[2]*se_total[1]))), 2)
        inv.est.l <- round(exp(mean(est_log) - (se_total[2]*se_total[1])) /
          (1 + exp(mean(est_log) - (se_total[2]*se_total[1]))), 2)
        res <- paste0(inv.est, " (", inv.est.l, ", ", inv.est.u, ")")
        
      } else if (exp_trans) { #if exp transformation requested for coefs (i.e., beta to hazard ratios)
        se_total <- RR_se(est, se) #RR pooling before exp transformation
        
        exp.est <- round(exp(mean(est)), 1)
        exp.est.u <- round(exp(mean(est) + (se_total[2]*se_total[1])), 1)
        exp.est.l <- round(exp(mean(est) - (se_total[2]*se_total[1])), 1)
      
        res <- paste0(exp.est, " (", exp.est.l, ", ", exp.est.u, ")")
        
      } else if (log_trans) { #if log transformation requested (e.g., for O/E ratios)
        est_log <- log(est) 
        se_log <- se #for O/E ratio, "se" here is se(obs risk)/obs risk (see supplementary page 9 of Debray et al, 2017); otherwise se(obs risk) / exp_risk
        se_total <- RR_se(est_log, se_log) # RR pooling after log transformation
        
        inv.est <- round(exp(mean(est_log)), 1) #back-transform
        inv.est.u <- round(exp(mean(est_log) + (se_total[2]*se_total[1])), 1)
        inv.est.l <- round(exp(mean(est_log) - (se_total[2]*se_total[1])), 1)
        res <- paste0(inv.est, " (", inv.est.l, ", ", inv.est.u, ")")
        
      } else {
        mean.est <- round(mean(est), 1)
        se_total <- RR_se(est, se)
        mean.est.u <- round(mean(est) + (se_total[2]*se_total[1]), 1)
        mean.est.l <- round(mean(est) - (se_total[2]*se_total[1]), 1)
        
        res <- paste0(mean.est, " (", mean.est.l, ", ", mean.est.u, ")")
      }
      return(res)
  }

#2. Calculated pooled coefficients and hazard ratios
  #data <- long form of the imputed dataset without the original dataset
  #impvar <- name of variable distinguishing the imputed datasets, in quotation marks
  #nimp <- number of imputed datasets
  #formula <- survival formula (e.g., Surv(time, event) ~ var 1 + var2)
  #exp_trans <- exponential transformation (e.g., betas -> HR)

  pool_coefs <- function(data, impvar, nimp, formula, exp_trans = TRUE) { #default gives hazard ratios
    coef_f <- se_f <- list()
    
    for (i in 1:nimp) { #function to calculate coefs, se, and HR
      data_compl <- data[data[impvar] == i, ]
      
      f <- coxph(formula, data = data_compl)
      f_sum <- summary(f)
      
      coef_f[[i]] <- f_sum$coefficients[, "coef"] #betas
      se_f[[i]] <- f_sum$coefficients[, "se(coef)"]
    }  
    
    perf_stats <- list(
      coef = do.call("rbind", coef_f),
      se = do.call("rbind", se_f)
    )
    
    pp_var <- function(df1, df2) { #pool for each covariate
      n <- ncol(df1)
      pooled_coef <- data.frame(matrix(nrow = n, ncol = 1))
      
      for (j in 1:n) {
        pooled_coef[j,1] <- pool_estimates(df1[, j], df2[, j], exp_trans = exp_trans)
      }
      
      rownames(pooled_coef) <- colnames(df1)
      return(pooled_coef)
    }
    
    # Pool
    res <- pp_var(perf_stats$coef, perf_stats$se)
    colnames(res) <- "pooled_HRs"
    return(res)
  }
  
#3. Calculate pooled Harrell's C and time-dependent AUC:
  #tvar <- time to event variable, in quotation marks
  #statvar <- event variable, in quotation marks
  #tmax <- evaluate until tmax time
  
  #NOTE: The C-statistic measures the ability of a model to correctly rank individuals according to their risk of experiencing an event (ratio of concordant to discordant pairs divided by total evaluable pairs). It ranges from 0.5 to 1. A value of 0.5 indicates a model that performs no better than random chance, while a value of 1 represents a perfect discriminatory model.
  #NOTE: The Area under the AUC Curve (AUC) measures the ability of the model to discriminate between positive and negative outcomes.
  #TLDR; The closer the C/AUC to 1 the better. C/AUC=0.5 means model is no better than random chance.
  
  pool_auc <- function(data, tvar, statvar, impvar, nimp, formula, tmax, logit_trans = FALSE) {
    perf_stats <- matrix(NA, nimp, 2)
    
    for (i in 1:nimp) {
      data_compl <- data[data[impvar] == i, ]
      f <- coxph(formula, data = data_compl)
     
      # Uno"s AUC up to tmax
      auc_temp <- timeROC(
        T = data_compl[, tvar],
        delta = data_compl[, statvar],
        cause = 1, #modify event indicator as preferred
        marker = predict(f, newdata = data_compl), #predict.coxph gives LP by default; here ,marker assumes larger values = higher risk of events, without loss of generality
        weighting = "marginal", 
        times = tmax,
        iid = TRUE) #very slow with large data
      auc_stats <- auc_temp$AUC[2]
      se_auc_stats <- auc_temp$inference$vect_sd_1[2]
      
      perf_stats[i, ] <- c(auc_stats, se_auc_stats)
    }
    
      # Pool
      auc_res <- pool_estimates(perf_stats[, 1], perf_stats[, 2], logit_trans = logit_trans)
      
      res <- list(pooled_AUC = auc_res)
    
      return(res)
  }
  
  pool_c <- function(data, tvar, statvar, impvar, nimp, formula, tmax, logit_trans = FALSE) {
    perf_stats <- matrix(NA, nimp, 2)
    
    for (i in 1:nimp) {
      data_compl <- data[data[impvar] == i, ]
      f <- coxph(formula, data = data_compl)
      
      # Harell"s C up to tmax
      c_temp <- concordance(f, ymax = tmax)
      c_stats <- c_temp$concordance
      se_c_stats <- sqrt(c_temp$var)
      
      perf_stats[i, ] <- c(c_stats, se_c_stats)
    }
    
    # Pool
    c_res <- pool_estimates(perf_stats[, 1], perf_stats[, 2], logit_trans = logit_trans)

    res <- list(pooled_C = c_res)
    
    return(res)
  }

#4. Calculate pooled calibration metrics:
  #NOTE: The calibration slope represents the change in predicted risk per unit change in observed risk (i.e., calibration-in-the-large). A calibration slope equal to 1 indicates perfect calibration, where the predicted risks are in line with the observed risks. Values greater than 1 suggest overestimation, while values less than 1 indicate underestimation of the event risk by the model.
  #NOTE: The O/E ratio compares the observed number of events to the expected number of events based on a prediction model (i.e., mean calibration). A value greater than 1 indicates higher observed events than expected, suggesting overprediction by the model, while a value less than 1 suggests lower observed events than expected, indicating underprediction. The O/E ratio provides insights into the accuracy and performance of the prediction model in estimating event probabilities.
  #TLDR; Calibration slope near 1 is represents no under/overestimation of risks. O/E ratio near represents good calibration on overage.
  
  pool_cal <- function(data, tvar, statvar, impvar, nimp, formula, tmax, log_trans = FALSE) { #log transformation applied to O/E ratio only
    perf_stats <- matrix(NA, nimp, 4)
    risk_list <- vector("list", nimp)

    for (i in 1:nimp) {
      data_compl <- data[data[impvar] == i, ]
      dat <- data_compl
      dat$subject <- seq_len(nrow(dat))
      f <- coxph(formula, data = dat)

      # Calibration slope
      lp <- f$linear.predictors #LP"s are automatically mean centered, where mean(categorical variable)=reference value
      cs_formula <- as.formula(paste("Surv(", tvar, ",", statvar, ") ~ ", paste("lp")))
      cs_fit <- coxph(cs_formula, data = dat)
      cs_stats <- cs_fit$coefficients
      se_cs_stats <- sqrt(diag(vcov(cs_fit)))
      
      # Observed risk at tmax
      surv_fit <- survfit(Surv(dat[, tvar], dat[, statvar]) ~ 1)
      obs_risk <- 1 - summary(surv_fit, times = tmax)$surv
      obs_se <- summary(surv_fit, times = tmax)$std.err
      obs_lrisk <- 1 - summary(surv_fit, times = tmax)$upper
      obs_urisk <- 1 - summary(surv_fit, times = tmax)$lower
      
      # Predicted risk at tmax
      bhaz <- basehaz(f, centered = TRUE) #basehaz centered by default
      bhaz_t <- max(bhaz$hazard[bhaz$time <= tmax]) #baseline cumhaz at time t
      bsurv_t <- exp(-bhaz_t) #baseline surv at time t: where S0 = exp(-H0 or cumhaz)
      dat$predrisk <- 1 - bsurv_t^(exp(lp)) #S(t) = S0^exp(LP), pred risk = 1-S(t) 
      exp_risk <- mean(dat$predrisk)
      
      # O/E ratio
      OE_stats <- obs_risk / exp_risk
      se_OE_stats <- obs_se/obs_risk #if log trans is required (see Debray 2017 formula)
      #se_OE_stats <- obs_se / exp_risk #if log trans not required
      
      risk_list[[i]] <- data.frame(subject = dat$subject,
                                   predrisk = dat$predrisk)
      
      perf_stats[i, ] <- c(cs_stats, se_cs_stats, OE_stats, se_OE_stats)
    }
    
    # Pool
    CS_res <- pool_estimates(perf_stats[, 1], perf_stats[, 2], log_trans = FALSE) #no transformation applied
    OE_res <- pool_estimates(perf_stats[, 3], perf_stats[, 4], log_trans = log_trans)
    
    risk_res  <- do.call(rbind, risk_list) %>%
      group_by(subject) %>%
      dplyr::summarize(
        predicted_risk = mean(predrisk))
    
    risk_df <- cbind(risk_res, time = data_compl[,tvar], event = data_compl[,statvar])
    
    res <- list(pooled_CalSlope = CS_res, pooled_OEratio = OE_res, pooled_risks = risk_df)
    
    return(res)
  }

#5. Create a calibration plot:
  #data <- data derived from pool_cal function
  #times <- calibration at t years
  #main <- plot title
  #C <- C-index
  #calslope <- calibration slope
  #OEratio <- OE ratio
  #size_legend <- size of texts
  #size_bintext <- size of texts next to histograms
  #line_bins <- position of histograms
  #triangles <- quantile means
  #g <- number of groups to define quantiles

  calplot_fn <- function(data, tmax, main = "", C, calslope, OEratio, limit = 0.5, size_lab = 1, size_legend = 0.45, size_bintext = 0.5, line_bins = 0, triangles = FALSE, g) {
    df <- data
    
    # Predicted risks
    df$x <- df$predicted_risk #predicted probabilities at tmax
    
    df$x.ll <- log(-log(1 - df$x)) #complementary log-log transformation of the predicted survival; improves linearity and lessens # of knots needed for rcs per Austin, et al, 2020
    
    model <- cph(Surv(time, event) ~ rcs(x.ll, 5), data = df, x = TRUE, y = TRUE, surv = TRUE)
    
    # Observed proportions
    xx <- seq(quantile(df$x, prob = 0.01), quantile(df$x, prob = 0.99), length = 100)
    xx.ll <- log(-log(1 - xx))
    xx.ll.df <- data.frame(x.ll = xx.ll)
    
    y <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$surv
    y.lower <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$upper
    y.upper <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$lower
    
    # Plot parameters
    xlim <- c(0, limit + 0.01)
    ylim <- c(-0.04, limit + 0.01)
    xlab <- "Predicted probability"
    ylab <- "Observed proportion"
    
    # Plot
    par(las = 1, xaxs = "i", yaxs = "i")
    plot(0, 0, type = "n",
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, cex.main = 1,
         cex.lab = size_lab, cex.axis = size_lab, lwd = 2)
    polygon(c(xx, rev(xx), xx[1]),
            c(y.lower, rev(y.upper), y.lower[1]),
            border = NA, density = 50, angle = -20, col = "gray")
    abline(coef = c(0,1), lty = 1, col = "gray") #diagonal line
    lines(xx, y, lwd = 2, col = "black")
    
    # Triangles
    if (triangles) {
      q <- Hmisc::cut2(df$x, levels.mean = TRUE, g=g) #group predicted risks
      means <- as.double(levels(q))
      y1 <- 1-survest(model, newdata=df, times = tmax)$surv
      prop <- tapply(y1, q, mean) #mean Observed proportions
      points(means, prop, pch=17, cex=1, col="black") #triangles
    }
    
    # Histograms:
    #rect(0,-0.025,1,0,lwd=0,border=NA)
    lim <- c(min(df$x), max(df$x))
    bins <- seq(lim[1], lim[2], length=101) #bins
    
    f0	<- table(cut(df$x[df$event == 0], bins)) #non-events
    f1	<- table(cut(df$x[df$event == 1], bins)) #events
    
    bins0 <- (bins[-101])[f0 > 0] #x bins, non-events
    bins1 <- (bins[-101])[f1 > 0] #x bins, events
    
    pcty1 <- as.numeric(f1/sum(f1))*0.1 #y bin, fix height multiplier for now
    pctx1 <- rep(0, length(pcty1))
    for (i in 1:length(pcty1)) { #histograms
      if (pcty1[i]>0) 
        rect(bins1[i], line_bins,
             bins1[i+1], line_bins + pcty1[i],
             col="gray", border=NA)
    }
    
    pcty0 <- as.numeric(f0/sum(f0))*0.1 #y bin, fix height multiplier for now
    pctx0 <- rep(0, length(pcty0))
    for (i in 1:length(pcty0)) { #histograms
      if (pcty0[i]>0)
        rect(bins0[i], line_bins,
             bins0[i+1], line_bins - pcty0[i],
             col="gray", lwd = 0.5, border=NA)
    }
    
    abline(h = line_bins, lty = 3)
    text(x = limit + 0.01, y = (line_bins + 0.005),
         labels = "Events",
         cex = size_bintext, pos = 2, col = "darkgray")
    text(x = limit + 0.01, y = (line_bins - 0.005),
         labels = "Non-events",
         cex = size_bintext, pos = 2, col = "darkgray")
    
    #Texts
    text(x=0, y=(limit), labels=paste("Discrimination"), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.01), labels=paste("...C: ", C, sep=""), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.02), labels=paste("Calibration "), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.03), labels=paste("...Slope: ", calslope, sep=""), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.04), labels=paste("...O/E ratio: ", OEratio, sep=""), cex=size_legend, pos=4)
  }
  
#6. Calculate Net Benefit (NB):
  #In survival settings, we can manually calculate NB over a range of thresholds using this formula (Vickers et al, 2008):
  #NB = TP/n -  w x FP/n, where:
  #w= pt/(1-pt), where that pt = threshold probability
  #TP = [1 - (S(t) | x = 1)] * P(x = 1) * n,
  #FP = (S(t) | x = 1) * P(x = 1) * n,
  #(x = 1) = if a patient has a predicted risk from the model b	% pt
  #S(t) = the Kaplan-Meier survival probability at time t
  #n = sample size
  
  #data <- data derived from pool_cal function
  
  nb_fn_TP <- function(data, tmax, thresholdmax = 1) {
    
    thresholds <- seq(0.01, thresholdmax, by = 0.001) #use by=0.001 due to low risks and more homogenous case mix; otherwise use 0.01
    
    NB_f <- lapply(thresholds, function(pt) {
      
      #NB treat all
      m_all <- 1 - summary(survfit(Surv(time, event) ~ 1, data = data), times = tmax)$surv #BC deaths
      NB_all <- m_all - (1-m_all) * (pt/(1-pt))
      
      #NB for model
      prop_pred <- nrow(subset(data, predicted_risk >= pt))/nrow(data)  #proportion predicted high risk from model; alternatively, mean(predrisk>=pt)
      
      surv <- try(
        summary(survfit(Surv(time, event) ~ 1, data = data[data$predicted_risk >= pt, ]), 
                times = tmax), silent = TRUE)
      
      if (class(surv) == "try-error") {
        TP <- 0
        FP <- 0
        NB <- 0 #no observations above threshold
      } else {
        m_model <- 1 - surv$surv
        TP <- m_model * prop_pred
        FP <- (1 - m_model) * prop_pred
        NB <- TP - FP * (pt/(1-pt))
      }
      
      NBres <- data.frame("threshold" = pt,
                          "NB_all" = NB_all,
                          "TP_all" = m_all,
                          "FP_all" = (1-m_all),
                          "NB" = NB,
                          "TP" = TP,
                          "FP" = FP)
    })
    
    #Bind results
    NB_res <- do.call(rbind, NB_f)
    
    return(NB_res)
  }
  
  nb_fn_TN <- function(data, tmax, thresholdmax = 1) { #note: labels can be confusing, refer to formula in manuscript
    
    thresholds <- seq(0.01, thresholdmax, by = 0.001)
    
    NB_f <- lapply(thresholds, function(pt) {
      
      #NB treat none
      m_none <- 1 - summary(survfit(Surv(time, event) ~ 1, data = data), times = tmax)$surv #BC deaths
      NB_none <- (1-m_none) - m_none * ((1-pt)/pt) #(1-m_none) is just surv prob estimate at tmax, weight is reversed
      
      #NB for model
      prop_pred <- nrow(subset(data, predicted_risk < pt))/nrow(data)  #proportion predicted low risks/below threshold
      
      surv <- try(
        summary(survfit(Surv(time, event) ~ 1, data = data[data$predicted_risk < pt, ]), 
                times = tmax), silent = TRUE)
      
      if (class(surv) == "try-error") {
        TN <- 1
        FN <- 1
        NB <- 1 #no observations below threshold
      } else {
        m_model <- surv$surv
        TN <- m_model * prop_pred
        FN <- (1 - m_model) * prop_pred
        NB <- TN - FN * ((1-pt)/pt)
      }
      
      NBres <- data.frame("threshold" = pt,
                          "NB_none" = NB_none,
                          "TN_none" = (1-m_none),
                          "FN_none" = m_none,
                          "NB" = NB,
                          "TN" = TN,
                          "FN" = FN)
    })
    
    #Bind results
    NB_res <- do.call(rbind, NB_f)
    
    return(NB_res)
  }
  
#7. Calculate pooled Likelihood Ratio (LR) estimates
  
  pool_lrt <- function(data, impvar, nimp, formula, nullformula) {
    perf_stats <- df <- matrix(NA, nimp, 1)
    #anova <- list()
    
      for (i in 1:nimp) {
      data_compl <- data[data[impvar] == i, ]

      f <- coxph(formula, data = data_compl)
      f.null <- coxph(nullformula, data = data_compl)
      
      #Calculate likelihood ratio (LR)
      perf_stats[i,] <- anova(f, f.null)$Chisq[2]
      df[i,] <- anova(f, f.null)$Df[2]
      }
    
    #Pool Chi-square statistics per D2 procedure
    ll.pool <- psfmi::pool_D2(perf_stats, df[1,])
    #alternatively: miceadds::micombine.chisquare(perf_stats, df)
    
    X2_res <- round(ll.pool["D2"], 3) #pooled Chi-square ~ F distribution
    LRT_res <- ll.pool["p"] #p-value

    res <- data.frame(pooled_X2 = X2_res, P_val = LRT_res)
    rownames(res) <- NULL
    return(res)
  }