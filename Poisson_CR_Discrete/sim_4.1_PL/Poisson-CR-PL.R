library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# ====================================================================
# 1(a). Data generation function
# 생존율(Survival Probability) 기반의 이산형 경쟁위험 데이터 생성 함수
# ====================================================================
generate_discrete_competing_risks <- 
  function(n = 500, 
           max_time = 10, 
           beta1 = c(0.5, -0.3), 
           beta2 = c(-0.2, 0.4)){
    
    # 공변량 생성 (Covariates)
    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, 0.5)
    X <- cbind(x1, x2)
    
    # 결과 저장 배열 초기화
    time <- numeric(n)
    status <- numeric(n)
    cause <- numeric(n)
    
    # 시간의 흐름에 따라 증가하는 기저 위험률 설정 (Baseline hazards)
    baseline_hazard1 <- seq(0.05, 0.15, length.out = max_time)
    baseline_hazard2 <- seq(0.03, 0.12, length.out = max_time)
    
    for (i in 1:n) {
      for (j in 1:max_time) {
        # 원인별 위험률(Cause-specific hazards) 계산
        hazard1 <- baseline_hazard1[j] * exp(X[i,] %*% beta1)
        hazard2 <- baseline_hazard2[j] * exp(X[i,] %*% beta2)
        total_hazard <- hazard1 + hazard2
        
        # 연속시간 생존율 구조를 반영한 정확한 조건부 이벤트 발생 식
        if (runif(1) > exp(-total_hazard)) {
          time[i] <- j
          status[i] <- 1
          
          # 이벤트 발생 시 해당 구간 내 원인별 상대 위험 비율로 원인 배정
          if (runif(1) < hazard1/total_hazard) {
            cause[i] <- 1
          } else {
            cause[i] <- 2
          }
          break
        }
      }
      # max_time까지 이벤트가 없으면 해당 시점에 중도절단(Censoring) 처리
      if (time[i] == 0) {
        time[i] <- max_time
        status[i] <- 0
        cause[i] <- 0
      }
    }
    
    # 원본 데이터 프레임 구축 (Wide Format)
    data <- data.frame(
      id = 1:n,
      time = time,
      status = status,
      cause = cause,
      x1 = x1,
      x2 = x2
    )
    
    return(data)
  }

# Generate data
n <- 500
data <- generate_discrete_competing_risks(n = n)

# Check the data
head(data)
table(data$status, data$cause)

# =========================================================================
# 1(b). Data generation function
# Lee et al. (2018) 방식 - 이산 logit 위험률 기반의 경쟁위험 데이터 생성 함수
# =========================================================================

generate_discrete_competing_risks_logit <- function(
    n = 500, max_time = 10, beta1 = c(0.5, -0.3), beta2 = c(-0.2, 0.4)
){
  
  # 공변량 생성 (Covariates)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  X <- cbind(x1, x2)
  
  # 결과 저장 배열 초기화
  time <- numeric(n)
  status <- numeric(n)
  cause <- numeric(n)
  
  # 시간의 흐름에 따라 증가하는 기저 위험률 설정 (Baseline hazards)
  baseline_hazard1 <- seq(0.05, 0.15, length.out = max_time)
  baseline_hazard2 <- seq(0.03, 0.12, length.out = max_time)
  alpha1 <- qlogis(baseline_hazard1)   # logit(p) = log(p/(1-p))
  alpha2 <- qlogis(baseline_hazard2)
  
  for (i in 1:n) {
    for (j in 1:max_time) {
      # 원인별 이산 위험률(Cause-specific discrete hazards)을 logit 모형으로 계산
      lambda1 <- plogis(alpha1[j] + sum(X[i,] * beta1))
      lambda2 <- plogis(alpha2[j] + sum(X[i,] * beta2))
      
      # 다항 구조: u < lambda1 -> 원인 1, lambda1 <= u < lambda1+lambda2 -> 원인 2
      u <- runif(1)
      if (u < lambda1) {
        time[i] <- j; status[i] <- 1; cause[i] <- 1
        break
      } else if (u < lambda1 + lambda2) {
        time[i] <- j; status[i] <- 1; cause[i] <- 2
        break
      }
    }
    # max_time까지 이벤트가 없으면 해당 시점에 중도절단(Censoring) 처리
    if (time[i] == 0) {
      time[i] <- max_time
      status[i] <- 0
      cause[i] <- 0
    }
  }
  
  # 원본 데이터 프레임 구축 (Wide Format)
  data <- data.frame(
    id = 1:n,
    time = time,
    status = status,
    cause = cause,
    x1 = x1,
    x2 = x2
  )
  
  return(data)
}

# Generate data
n <- 500
data <- generate_discrete_competing_risks_logit(n = n)

# Check the data
head(data)
table(data$status, data$cause)

# ============================================================
# 2. Long format expansion function
# 분석을 위한 롱 포맷 데이터 변환 함수 (Counting Process Format)
# ============================================================
create_expanded_data <- function(data) {
  expanded_list <- list()
  row_idx <- 1
  
  for (i in 1:nrow(data)) {
    subject <- data[i, ]
    event_time <- subject$time
    event_status <- subject$status
    event_cause <- subject$cause
    
    # 각 개체가 거쳐간 시간 구간(time_interval)만큼 행을 생성
    for (j in 1:event_time) {
      # 해당 구간이 최종 시점이면서 원인 1 또는 2를 겪었는지 여부
      ev1 <- ifelse(j == event_time & event_cause == 1, 1, 0)
      ev2 <- ifelse(j == event_time & event_cause == 2, 1, 0)
      
      expanded_list[[row_idx]] <- data.frame(
        id = subject$id,
        time = event_time,      # 최종 관측 시간
        time_interval = j,      # 현재 지나가고 있는 구간 인덱스
        x1 = subject$x1,
        x2 = subject$x2,
        event1 = ev1,
        event2 = ev2,
        # 논문의 delta^c 정의: 해당 구간에서 원인 1, 2 모두 발생하지 않고 생존해 낸 경우
        censored = 1 - ev1 - ev2
      )
      row_idx <- row_idx + 1
    }
  }
  
  # 리스트를 하나의 데이터 프레임으로 결합 (rbind 반복보다 속도가 빠름)
  return(do.call(rbind, expanded_list))
}

# ============================================================
# 3. Proposed partial likelihood estimation
# R 기본 문법만을 사용하여 롱 데이터로 부분 우도를 추정하는 함수
# ============================================================
partial_likelihood_cr_pure_long <- function(long_data, cause_type = 1) {
  
  # 해당 원인의 이벤트 타깃 변수명 지정 ('event1' 또는 'event2')
  event_col <- paste0("event", cause_type)
  
  # 오직 long_data 내에서 해당 원인의 이벤트가 실제 발생한 시점(time_interval)들만 추출
  event_rows <- long_data[long_data[, event_col] == 1, ]
  event_times <- sort(unique(event_rows$time_interval))
  
  # 최적화를 위한 초기 설정
  beta <- c(0, 0)  # 공변량 x1, x2의 초기 계수 값
  max_iter <- 20
  tolerance <- 1e-4
  
  for(iter in 1:max_iter) {
    score <- c(0, 0)
    information <- matrix(0, nrow = 2, ncol = 2)
    
    for(t in event_times) {
      
      # 현재 분석 시점 t에 해당하는 행들만 슬라이싱 (Risk set 후보군)
      current_set <- long_data[long_data$time_interval == t, ]
      
      # 1) 앞 항 계산 대상: 현재 시점 t에서 해당 원인 이벤트가 발생한 개체들 (delta_kil == 1)
      events_df <- current_set[current_set[, event_col] == 1, ]
      
      # 2) 뒷 항(위험집합) 계산 대상: 해당 구간 t를 성공적으로 생존해 낸 개체들 (delta_jl^c == 1)
      at_risk_df <- current_set[current_set$censored == 1, ]
      
      n_events <- nrow(events_df)
      
      if(n_events > 0 && nrow(at_risk_df) > 0) {
        
        # 디자인 행렬 추출 (Matrix 변환)
        x_event <- as.matrix(events_df[, c("x1", "x2")])
        x_risk  <- as.matrix(at_risk_df[, c("x1", "x2")])
        
        # 선형 예측량 및 웨이트(exp_eta) 계산
        eta_risk <- x_risk %*% beta
        exp_eta  <- as.vector(exp(eta_risk))
        
        sum_exp_eta   <- sum(exp_eta)                            # 위험집합 분모
        sum_x_exp_eta <- colSums(x_risk * exp_eta)               # 위험집합 분자 (공변량 가중합)
        
        # 정보 행렬(Information Matrix)을 위한 xx' * exp(eta) 누적 계산
        xx_exp_eta <- matrix(0, 2, 2)
        for(i in 1:nrow(x_risk)) {
          xx <- x_risk[i,] %*% t(x_risk[i,])
          xx_exp_eta <- xx_exp_eta + xx * exp_eta[i]
        }
        
        # 스코어 함수 및 정보 행렬 업데이트 (R 기본 연산)
        score <- score + colSums(x_event) - n_events * (sum_x_exp_eta / sum_exp_eta)
        
        E_xx <- xx_exp_eta / sum_exp_eta
        V_x  <- (sum_x_exp_eta %*% t(sum_x_exp_eta)) / (sum_exp_eta^2)
        information <- information + n_events * (E_xx - V_x)
      }
    }
    
    # Newton-Raphson 방식을 통한 beta 업데이트
    delta <- solve(information) %*% score
    beta <- beta + delta
    
    # 수렴 조건 검사
    if(sqrt(sum(delta^2)) < tolerance) {
      break
    }
  }
  
  # 표준오차 및 검정 통계량 산출
  se <- sqrt(diag(solve(information)))
  result_df <- data.frame(
    coefficient = c("x1", "x2"),
    true_beta   = if(cause_type == 1) c(1, -1) else c(-1, 1),
    estimate    = beta,
    std_error   = se,
    z_value     = beta/se,
    p_value     = 2 * pnorm(-abs(beta/se))
  )
  
  return(list(estimates = result_df, beta = beta, information = information))
}

# ============================================================
# [추가] Cluster-robust sandwich SE 계산 함수
# 원본 partial_likelihood_cr_pure_long()의 결과(fit)를 받아
# 개체 단위 score 기여 U_i를 모아 I^{-1} V I^{-1}를 계산
# ============================================================
compute_robust_se <- function(long_data, cause_type, fit) {
  
  event_col <- paste0("event", cause_type)
  beta <- as.vector(fit$beta)
  
  event_rows <- long_data[long_data[, event_col] == 1, ]
  event_times <- sort(unique(event_rows$time_interval))
  
  n_id <- max(long_data$id)
  U <- matrix(0, nrow = n_id, ncol = 2)   # 개체별 score 기여
  
  for(t in event_times) {
    current_set <- long_data[long_data$time_interval == t, ]
    events_df   <- current_set[current_set[, event_col] == 1, ]
    at_risk_df  <- current_set[current_set$censored == 1, ]
    
    n_events <- nrow(events_df)
    if(n_events > 0 && nrow(at_risk_df) > 0) {
      x_event <- as.matrix(events_df[, c("x1", "x2")])
      x_risk  <- as.matrix(at_risk_df[, c("x1", "x2")])
      exp_eta <- as.vector(exp(x_risk %*% beta))
      
      sum_exp_eta <- sum(exp_eta)
      xbar <- colSums(x_risk * exp_eta) / sum_exp_eta
      
      # 이벤트 항: delta_kil * (x_i - xbar_l)
      for(j in 1:n_events) {
        i_id <- events_df$id[j]
        U[i_id, ] <- U[i_id, ] + (x_event[j, ] - xbar)
      }
      
      # 위험집합 항 (Breslow형 보정): - d_l * w_il / S0_l * (x_i - xbar_l)
      wgt <- n_events * exp_eta / sum_exp_eta
      for(j in 1:nrow(x_risk)) {
        i_id <- at_risk_df$id[j]
        U[i_id, ] <- U[i_id, ] - wgt[j] * (x_risk[j, ] - xbar)
      }
    }
  }
  
  V <- t(U) %*% U
  I_inv <- solve(fit$information)
  se_robust <- sqrt(diag(I_inv %*% V %*% I_inv))
  
  return(se_robust)
}


# ============================================================
# [추가] Meir & Gorfine (2025) two-step estimator
# ============================================================

gorfine_alpha_step <- function(data, beta_hat, cause,
                               covariates = c("x1", "x2"),
                               times = NULL,
                               lower = -50, upper = 50) {

  Z   <- as.matrix(data[, covariates, drop = FALSE])
  eta <- as.vector(Z %*% beta_hat)

  if (is.null(times)) times <- sort(unique(data$time))

  alpha <- rep(NA_real_, length(times))

  for (m in seq_along(times)) {
    tt      <- times[m]
    at_risk <- data$time >= tt                              # I(X_i >= t)
    Yt      <- sum(at_risk)                                 # Y.(t)
    Njt     <- sum(data$time == tt & data$cause == cause)   # N_j(t)

    # N_j(t) = 0 이면 alpha -> -Inf, N_j(t) = Y.(t) 이면 alpha -> +Inf : 식별 불가
    if (Yt == 0 || Njt == 0 || Njt >= Yt) next

    g <- function(a) sum(plogis(a + eta[at_risk])) - Njt    # a 에 대해 순증가

    alpha[m] <- tryCatch(
      uniroot(g, lower = lower, upper = upper, tol = 1e-10)$root,
      error = function(e) NA_real_
    )
  }

  data.frame(time = times, alpha = alpha)
}


gorfine_two_step <- function(data, long_data,
                             causes = c(1, 2),
                             covariates = c("x1", "x2"),
                             ties = "exact",           # "exact" = 식 (4) 그대로
                             estimate_alpha = FALSE) {

  stopifnot(all(covariates %in% names(long_data)))
  stopifnot("time_interval" %in% names(long_data))

  out <- list()

  for (k in causes) {
    ev  <- paste0("event", k)
    if (!ev %in% names(long_data))
      stop(sprintf("long_data 에 '%s' 열이 없습니다.", ev))

    # ---- Step 1: conditional logistic, strata = 시점 ----
    fml <- as.formula(
      paste0(ev, " ~ ", paste(covariates, collapse = " + "),
             " + strata(time_interval)")
    )
    fit <- clogit(fml, data = long_data, method = ties)

    beta_hat <- coef(fit)
    se_hat   <- sqrt(diag(vcov(fit)))
    names(beta_hat) <- names(se_hat) <- covariates

    # ---- Step 2: alpha_jt (식 5) ----
    alpha_hat <- NULL
    if (estimate_alpha) {
      alpha_hat <- gorfine_alpha_step(data, beta_hat, cause = k,
                                      covariates = covariates)
    }

    out[[paste0("cause", k)]] <- list(
      beta  = beta_hat,
      se    = se_hat,
      alpha = alpha_hat,
      fit   = fit
    )
  }

  out
}


# ============================================================
# 5. Simulation function  (기존 구조 그대로 + Method 4: gorfine)
# ============================================================

simulate_comparison <- function(n = 300, max_time = 10, n_sim = 100,
                                beta1_true = c(0.5, -0.3),
                                beta2_true = c(-0.2, 0.4),
                                dgp_fun = generate_discrete_competing_risks_log,
                                gorfine_ties = "exact",
                                verbose_every = 20) {

  results <- list()

  for (sim in 1:n_sim) {
    data      <- dgp_fun(n, max_time, beta1_true, beta2_true)
    long_data <- create_expanded_data(data)

    # ---- Method 1: Logistic (collapsed likelihood, Lee et al. 2018) ----
    logistic_results <- tryCatch({
      fit1_logistic <- glm(event1 ~ x1 + x2 + factor(time_interval),
                           family = binomial, data = long_data)
      fit2_logistic <- glm(event2 ~ x1 + x2 + factor(time_interval),
                           family = binomial, data = long_data)
      list(
        cause1 = summary(fit1_logistic)$coefficients[2:3, ],
        cause2 = summary(fit2_logistic)$coefficients[2:3, ]
      )
    }, error = function(e) NULL)

    # ---- Method 2: Traditional Cox ----
    cox_results <- tryCatch({
      fit1_cox <- coxph(Surv(time, cause == 1) ~ x1 + x2, data = data,
                        ties = "breslow")
      fit2_cox <- coxph(Surv(time, cause == 2) ~ x1 + x2, data = data,
                        ties = "breslow")
      list(
        cause1 = cbind(coef(fit1_cox), sqrt(diag(vcov(fit1_cox)))),
        cause2 = cbind(coef(fit2_cox), sqrt(diag(vcov(fit2_cox))))
      )
    }, error = function(e) NULL)

    # ---- Method 3: Proposed partial likelihood (naive + sandwich SE) ----
    discrete_results <- tryCatch({
      fit1 <- partial_likelihood_cr_pure_long(long_data, cause_type = 1)
      fit2 <- partial_likelihood_cr_pure_long(long_data, cause_type = 2)
      rob1 <- compute_robust_se(long_data, 1, fit1)
      rob2 <- compute_robust_se(long_data, 2, fit2)
      list(
        cause1 = cbind(fit1$estimates$estimate, fit1$estimates$std_error, rob1),
        cause2 = cbind(fit2$estimates$estimate, fit2$estimates$std_error, rob2)
      )
    }, error = function(e) NULL)

    # ---- Method 4: Meir & Gorfine (2025) two-step ----
    #      표에 쓰는 것은 beta 뿐이므로 루프에서는 alpha 단계를 끔.
    gorfine_results <- tryCatch({
      g <- gorfine_two_step(data, long_data,
                            causes = c(1, 2),
                            ties = gorfine_ties,
                            estimate_alpha = FALSE)
      list(
        cause1 = cbind(g$cause1$beta, g$cause1$se),
        cause2 = cbind(g$cause2$beta, g$cause2$se)
      )
    }, error = function(e) NULL)

    results[[sim]] <- list(
      logistic = logistic_results,
      cox      = cox_results,
      discrete = discrete_results,
      gorfine  = gorfine_results
    )

    if (sim %% verbose_every == 0) cat("Simulation", sim, "completed\n")
  }

  return(results)
}


# ============================================================
# 6. Analysis function  (gorfine 추가)
# ============================================================

analyze_simulation <- function(results,
                               beta1_true = c(0.5, -0.3),
                               beta2_true = c(-0.2, 0.4),
                               methods = c("logistic", "gorfine", "cox",
                                           "discrete_robust")) {

  causes     <- c("cause1", "cause2")
  covariates <- c("x1", "x2")

  true_beta <- list(cause1 = beta1_true, cause2 = beta2_true)

  performance_metrics <- expand.grid(
    method    = methods,
    cause     = causes,
    covariate = covariates,
    stringsAsFactors = FALSE
  )

  performance_metrics$bias      <- NA
  performance_metrics$ASE       <- NA
  performance_metrics$ESE       <- NA
  performance_metrics$CP        <- NA
  performance_metrics$n_success <- 0

  for (i in 1:nrow(performance_metrics)) {
    method        <- performance_metrics$method[i]
    cause         <- performance_metrics$cause[i]
    cov_idx       <- ifelse(performance_metrics$covariate[i] == "x1", 1, 2)
    true_beta_val <- true_beta[[cause]][cov_idx]

    # discrete / discrete_robust 는 같은 적합, SE 열만 다름
    fit_slot <- ifelse(method == "discrete_robust", "discrete", method)
    se_col   <- ifelse(method == "discrete_robust", 3, 2)

    estimates  <- numeric(0)
    std_errors <- numeric(0)

    for (sim in 1:length(results)) {
      if (!is.null(results[[sim]][[fit_slot]])) {
        est <- results[[sim]][[fit_slot]][[cause]][cov_idx, 1]
        se  <- results[[sim]][[fit_slot]][[cause]][cov_idx, se_col]
        estimates  <- c(estimates, est)
        std_errors <- c(std_errors, se)
      }
    }

    if (length(estimates) > 0) {
      performance_metrics$bias[i] <- mean(estimates - true_beta_val)
      performance_metrics$ASE[i]  <- mean(std_errors)
      performance_metrics$ESE[i]  <- sd(estimates)

      ci_lower <- estimates - 1.96 * std_errors
      ci_upper <- estimates + 1.96 * std_errors
      performance_metrics$CP[i] <-
        mean(ci_lower <= true_beta_val & true_beta_val <= ci_upper)

      performance_metrics$n_success[i] <- length(estimates)
    }
  }

  return(performance_metrics)
}

# ====================
# 7. Run simulation
# ====================

cat("Starting simulation...\n")
sim_results <- simulate_comparison(n = 300, max_time = 10, n_sim = 100)

# Analyze results
performance <- analyze_simulation(sim_results)

# Print results
print("Performance Comparison:")
print(performance)

# Create summary table
summary_table <- performance %>%
  mutate(
    Bias = round(bias, 4),
    ASE = round(ASE, 4),
    ESE = round(ESE, 4),
    CP = round(CP, 4),
    `Bias/ESE` = round(abs(bias)/ESE, 4)
  ) %>%
  select(method, cause, covariate, Bias, ASE, ESE, CP, `Bias/ESE`)

print("Summary Table:")
print(summary_table)

# Visualization
ggplot(performance, aes(x = method, y = abs(bias), fill = method)) +
  geom_bar(stat = "identity") +
  facet_grid(cause ~ covariate) +
  labs(title = "Absolute Bias by Method",
       y = "Absolute Bias", x = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot(performance, aes(x = method, y = CP, fill = method)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_grid(cause ~ covariate) +
  labs(title = "Coverage Probability by Method",
       y = "Coverage Probability", x = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))




# ====================
# [추가] M&G Sanity check 
# ====================

check_gorfine <- function(n = 500, max_time = 10, seed = 123,
                          beta1 = c(0.5, -0.3), beta2 = c(-0.2, 0.4),
                          ties = "exact") {

  set.seed(seed)
  data <- generate_discrete_competing_risks(n, max_time, beta1, beta2)
  ld   <- create_expanded_data(data)

  cat("=== (0) 확장자료 구조 확인 ===\n")
  cat("  nrow(long_data) =", nrow(ld),
      " / sum(data$time) =", sum(data$time),
      ifelse(nrow(ld) == sum(data$time), "  [OK]", "  [!! 확인 필요]"), "\n")
  cat("  t 별 위험집합 크기 = Y.(t) 일치 여부: ")
  ok <- all(sapply(1:max_time, function(t)
    sum(ld$time_interval == t) == sum(data$time >= t)))
  cat(ifelse(ok, "OK", "!! 불일치"), "\n\n")

  g <- gorfine_two_step(data, ld, ties = ties, estimate_alpha = TRUE)

  # 비교 대상: Lee et al. collapsed logistic (동일한 logit 링크 모수를 목표)
  f1 <- glm(event1 ~ x1 + x2 + factor(time_interval) - 1,
            family = binomial, data = ld)
  f2 <- glm(event2 ~ x1 + x2 + factor(time_interval) - 1,
            family = binomial, data = ld)

  cat("=== (1) beta: Gorfine two-step vs Lee et al. (2018) ===\n")
  cmp <- data.frame(
    cause     = rep(c(1, 2), each = 2),
    covariate = rep(c("x1", "x2"), 2),
    true      = c(beta1, beta2),
    gorfine   = c(g$cause1$beta, g$cause2$beta),
    gorfine_se= c(g$cause1$se,   g$cause2$se),
    lee       = c(coef(f1)[c("x1", "x2")], coef(f2)[c("x1", "x2")]),
    lee_se    = c(sqrt(diag(vcov(f1)))[c("x1", "x2")],
                  sqrt(diag(vcov(f2)))[c("x1", "x2")])
  )
  print(cmp, row.names = FALSE, digits = 4)
  cat("  -> 두 방법은 같은 logit-link 모수를 추정하므로 값이 매우 가까워야 정상.\n")
  cat("     (논문 Table 1 에서도 two-step 추정치가 Lee et al. 대비\n")
  cat("      절댓값이 아주 조금 작게 나오는 패턴)\n\n")

  cat("=== (2) alpha_1t: Gorfine 식(5) vs Lee et al. GLM ===\n")
  a_lee <- coef(f1)[grep("factor", names(coef(f1)))]
  acmp <- data.frame(
    time    = g$cause1$alpha$time,
    gorfine = g$cause1$alpha$alpha,
    lee     = as.vector(a_lee)[seq_len(nrow(g$cause1$alpha))]
  )
  print(acmp, row.names = FALSE, digits = 4)
  cat("  -> 논문 Figure 1 과 같은 비교. 두 곡선이 거의 겹쳐야 정상.\n")

  invisible(list(data = data, long = ld, gorfine = g, beta_cmp = cmp,
                 alpha_cmp = acmp))
}


benchmark_gorfine_ties <- function(n = 500, max_time = 10, seed = 1) {
  set.seed(seed)
  data <- generate_discrete_competing_risks(n, max_time,
                                            c(0.5, -0.3), c(-0.2, 0.4))
  ld   <- create_expanded_data(data)

  res <- NULL
  for (m in c("exact", "breslow", "efron")) {
    tm <- system.time({
      g <- try(gorfine_two_step(data, ld, ties = m, estimate_alpha = FALSE),
               silent = TRUE)
    })["elapsed"]
    b <- if (inherits(g, "try-error")) c(NA, NA) else g$cause1$beta
    res <- rbind(res, data.frame(ties = m, sec_per_rep = as.numeric(tm),
                                 beta11 = b[1], beta12 = b[2]))
  }
  cat("1회 적합 소요시간 (cause 2개 합산):\n")
  print(res, row.names = FALSE, digits = 4)
  cat("\n  n_sim = 1000 환산:", 
      paste0(round(res$sec_per_rep * 1000 / 60, 1), "분", collapse = " / "), "\n")
  cat("  -> exact 가 너무 느리면 simulate_comparison(gorfine_ties = \"breslow\")\n")
  invisible(res)
}


# ============================================================
# [추가] (선택) alpha 비교 시뮬레이션 — 논문 Figure 1 형태
# ============================================================

simulate_alpha_comparison <- function(n = 500, max_time = 10, n_sim = 100,
                                      beta1_true = c(0.5, -0.3),
                                      beta2_true = c(-0.2, 0.4),
                                      dgp_fun = generate_discrete_competing_risks,
                                      ties = "exact") {

  A <- vector("list", n_sim)

  for (s in 1:n_sim) {
    data <- dgp_fun(n, max_time, beta1_true, beta2_true)
    ld   <- create_expanded_data(data)

    A[[s]] <- tryCatch({
      g <- gorfine_two_step(data, ld, ties = ties, estimate_alpha = TRUE)

      f1 <- glm(event1 ~ x1 + x2 + factor(time_interval) - 1,
                family = binomial, data = ld)
      f2 <- glm(event2 ~ x1 + x2 + factor(time_interval) - 1,
                family = binomial, data = ld)
      lee1 <- as.vector(coef(f1)[grep("factor", names(coef(f1)))])
      lee2 <- as.vector(coef(f2)[grep("factor", names(coef(f2)))])

      rbind(
        data.frame(sim = s, cause = 1, time = g$cause1$alpha$time,
                   method = "gorfine", alpha = g$cause1$alpha$alpha),
        data.frame(sim = s, cause = 2, time = g$cause2$alpha$time,
                   method = "gorfine", alpha = g$cause2$alpha$alpha),
        data.frame(sim = s, cause = 1, time = seq_along(lee1),
                   method = "lee", alpha = lee1),
        data.frame(sim = s, cause = 2, time = seq_along(lee2),
                   method = "lee", alpha = lee2)
      )
    }, error = function(e) NULL)

    if (s %% 20 == 0) cat("alpha sim", s, "completed\n")
  }

  A <- do.call(rbind, A)

  summ <- A %>%
    group_by(cause, method, time) %>%
    summarise(mean_alpha = mean(alpha, na.rm = TRUE),
              sd_alpha   = sd(alpha, na.rm = TRUE),
              .groups = "drop")

  p <- ggplot(summ, aes(x = time, y = mean_alpha,
                        colour = method, shape = method)) +
    geom_point(size = 2) + geom_line() +
    facet_wrap(~ cause, labeller = label_both) +
    labs(title = sprintf("Estimated alpha_jt (B = %d)", n_sim),
         x = "Time", y = expression(hat(alpha)[jt])) +
    theme_minimal()

  print(p)
  invisible(list(raw = A, summary = summ, plot = p))
}
