################################################################################
#                    ICLO100 RSI TRADING STRATEGY                              #
#              ICE Gasoil Futures - RSI Momentum/Reversal                      #
################################################################################
# RSI (Relative Strength Index) - Developed by J. Welles Wilder
# 
# Strategy Logic:
# - RSI < 30 → Oversold → BUY (expecting bounce)
# - RSI > 70 → Overbought → SELL (expecting pullback)
# - RSI can also identify trends (persistently high/low)
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(zoo)
})

################################################################################
#                           CONFIGURATION                                       #
################################################################################

DATA_FILE <- "iclo_update2.csv"

# Period definitions
TRAIN_END_DATE    <- as.Date("2025-07-31")
VAL_START_DATE    <- as.Date("2025-08-01")
VAL_END_DATE      <- as.Date("2025-08-31")
TEST_START_DATE   <- as.Date("2025-09-01")
TEST_END_DATE     <- as.Date("2025-09-15")

# RSI Strategy Parameters
RSI_PERIOD        <- 10      # Standard RSI lookback (14 days)
RSI_OVERSOLD      <- 30      # Buy threshold (classic oversold)
RSI_OVERBOUGHT    <- 65      # Sell threshold (classic overbought)
RSI_EXIT_LOW      <- 50      # Exit long when RSI crosses above 50
RSI_EXIT_HIGH     <- 50      # Exit short when RSI crosses below 50

TRANSACTION_COST  <- 0.0010  # 10 bps per trade

# Parameter search grid for optimization
RSI_PERIOD_OPTIONS    <- c(7, 10, 14, 21, 28)
OVERSOLD_OPTIONS      <- c(20, 25, 30, 35)
OVERBOUGHT_OPTIONS    <- c(65, 70, 75, 80)

################################################################################
#                           UTILITY FUNCTIONS                                   #
################################################################################

calc_cagr <- function(returns, periods_per_year = 252) {
  if (length(returns) < 2 || all(is.na(returns))) return(NA_real_)
  equity_curve <- cumprod(1 + returns)
  total_return <- tail(equity_curve, 1)
  n_years <- length(returns) / periods_per_year
  if (n_years <= 0 || !is.finite(total_return)) return(NA_real_)
  (total_return ^ (1/n_years)) - 1
}

calc_ann_vol <- function(returns, periods_per_year = 252) {
  sd(returns, na.rm = TRUE) * sqrt(periods_per_year)
}

calc_sharpe <- function(returns, periods_per_year = 252, rf = 0) {
  ann_ret <- calc_cagr(returns, periods_per_year)
  ann_vol <- calc_ann_vol(returns, periods_per_year)
  if (!is.finite(ann_vol) || ann_vol == 0) return(NA_real_)
  (ann_ret - rf) / ann_vol
}

calc_max_drawdown <- function(returns) {
  if (length(returns) < 2) return(0)
  equity_curve <- cumprod(1 + returns)
  equity_curve <- equity_curve[is.finite(equity_curve) & equity_curve > 0]
  if (length(equity_curve) < 2) return(0)
  peak <- cummax(equity_curve)
  drawdown <- (equity_curve - peak) / peak
  min(drawdown, na.rm = TRUE)
}

calc_win_rate <- function(returns) {
  trades <- returns[returns != 0]
  if (length(trades) == 0) return(NA_real_)
  mean(trades > 0, na.rm = TRUE)
}

count_trades <- function(positions) {
  changes <- diff(c(0, positions))
  sum(changes != 0) / 2
}

calc_avg_holding_period <- function(positions) {
  if (length(positions) < 2 || all(positions == 0)) return(NA_real_)
  rle_result <- rle(positions)
  holding_periods <- rle_result$lengths[rle_result$values != 0]
  if (length(holding_periods) == 0) return(NA_real_)
  mean(holding_periods)
}

count_drawdown_periods <- function(returns) {
  if (length(returns) < 2) return(0)
  equity_curve <- cumprod(1 + returns)
  peak <- cummax(equity_curve)
  in_drawdown <- equity_curve < peak
  sum(diff(c(FALSE, in_drawdown)) == 1)
}

calc_performance_metrics <- function(returns, positions, period_name = "") {
  list(
    Period = period_name,
    Annual_Return = calc_cagr(returns, 252),
    Annual_Volatility = calc_ann_vol(returns, 252),
    Sharpe_Ratio = calc_sharpe(returns, 252),
    Max_Drawdown = calc_max_drawdown(returns),
    Win_Rate = calc_win_rate(returns),
    Number_of_Trades = count_trades(positions),
    Avg_Holding_Period = calc_avg_holding_period(positions),
    Drawdown_Periods = count_drawdown_periods(returns),
    Total_Return = prod(1 + returns, na.rm = TRUE) - 1
  )
}

################################################################################
#                           RSI CALCULATION                                     #
################################################################################

calculate_rsi <- function(prices, n = 14) {
  # Calculate price changes
  price_changes <- c(NA, diff(prices))
  
  # Separate gains and losses
  gains <- ifelse(price_changes > 0, price_changes, 0)
  losses <- ifelse(price_changes < 0, -price_changes, 0)
  
  # Calculate average gains and losses using Wilder's smoothing
  avg_gain <- numeric(length(prices))
  avg_loss <- numeric(length(prices))
  
  # Initial averages (simple moving average for first n periods)
  avg_gain[n] <- mean(gains[1:n], na.rm = TRUE)
  avg_loss[n] <- mean(losses[1:n], na.rm = TRUE)
  
  # Subsequent values use Wilder's smoothing (exponential)
  for (i in (n+1):length(prices)) {
    avg_gain[i] <- (avg_gain[i-1] * (n-1) + gains[i]) / n
    avg_loss[i] <- (avg_loss[i-1] * (n-1) + losses[i]) / n
  }
  
  # Calculate RS (Relative Strength)
  rs <- avg_gain / (avg_loss + 1e-10)  # Add small number to avoid division by zero
  
  # Calculate RSI
  rsi <- 100 - (100 / (1 + rs))
  
  # Set first n-1 values to NA
  rsi[1:(n-1)] <- NA
  
  return(rsi)
}

################################################################################
#                           DATA LOADING                                        #
################################################################################

cat(strrep("=", 80), "\n")
cat("ICLO100 RSI TRADING STRATEGY - ICE Gasoil Futures\n")
cat(strrep("=", 80), "\n\n")

cat("Loading data from:", DATA_FILE, "\n")
df_raw <- read.csv(DATA_FILE, stringsAsFactors = FALSE)
df_raw$assessDate <- as.Date(df_raw$assessDate)

cat("Data range:", as.character(min(df_raw$assessDate)), "to", 
    as.character(max(df_raw$assessDate)), "\n")
cat("Total observations:", nrow(df_raw), "\n\n")

################################################################################
#                      DATA PREPROCESSING                                       #
################################################################################

df <- df_raw %>%
  arrange(assessDate) %>%
  rename(date = assessDate, close = c, high = h, low = l) %>%
  mutate(
    log_return = c(NA, diff(log(close))),
    simple_return = close / lag(close) - 1,
    
    # Calculate RSI
    rsi = calculate_rsi(close, RSI_PERIOD),
    
    # Period classification
    period = case_when(
      date <= TRAIN_END_DATE ~ "Training",
      date >= VAL_START_DATE & date <= VAL_END_DATE ~ "Validation",
      date >= TEST_START_DATE & date <= TEST_END_DATE ~ "Test",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(log_return))

cat("After preprocessing:", nrow(df), "observations with returns\n")
cat("RSI calculation requires", RSI_PERIOD, "days lookback\n\n")

################################################################################
#                        SIGNAL GENERATION                                      #
################################################################################

# RSI Strategy Logic:
# Entry:
#   - Long when RSI < 30 (oversold, expect bounce)
#   - Short when RSI > 70 (overbought, expect pullback)
# Exit:
#   - Exit long when RSI > 50 (momentum shifted)
#   - Exit short when RSI < 50 (momentum shifted)

df <- df %>%
  mutate(position = 0)

current_position <- 0
for (i in 1:nrow(df)) {
  if (is.na(df$rsi[i])) {
    df$position[i] <- 0
    current_position <- 0
    next
  }
  
  rsi_val <- df$rsi[i]
  
  # Entry signals (when not in position)
  if (current_position == 0) {
    if (rsi_val < RSI_OVERSOLD) {
      current_position <- 1    # Enter long (oversold)
    } else if (rsi_val > RSI_OVERBOUGHT) {
      current_position <- -1   # Enter short (overbought)
    }
  }
  # Exit signals (when in position)
  else if (current_position == 1) {
    # Exit long if RSI crosses above midpoint
    if (rsi_val > RSI_EXIT_LOW) {
      current_position <- 0
    }
  }
  else if (current_position == -1) {
    # Exit short if RSI crosses below midpoint
    if (rsi_val < RSI_EXIT_HIGH) {
      current_position <- 0
    }
  }
  
  df$position[i] <- current_position
}

################################################################################
#                        RETURNS CALCULATION                                    #
################################################################################

df <- df %>%
  mutate(
    position_lagged = lag(position, 1, default = 0),
    position_change = abs(position_lagged - lag(position_lagged, 1, default = 0)),
    strategy_return = position_lagged * log_return - 
      TRANSACTION_COST * position_change,
    cum_strategy = cumprod(1 + replace_na(strategy_return, 0)),
    cum_bh = cumprod(1 + replace_na(log_return, 0)),
    running_max = cummax(cum_strategy),
    drawdown = (cum_strategy - running_max) / running_max
  )

################################################################################
#                     DIAGNOSTIC OUTPUT                                         #
################################################################################

cat(strrep("=", 80), "\n")
cat("DIAGNOSTIC INFORMATION\n")
cat(strrep("=", 80), "\n\n")

# RSI statistics by period
for (period_name in c("Training", "Validation", "Test")) {
  df_period <- df %>% filter(period == period_name)
  if (nrow(df_period) > 0) {
    rsi_valid <- df_period$rsi[!is.na(df_period$rsi)]
    cat(sprintf("%s Period (n=%d):\n", period_name, length(rsi_valid)))
    if (length(rsi_valid) > 0) {
      cat(sprintf("  RSI range: [%.1f, %.1f]\n", min(rsi_valid), max(rsi_valid)))
      cat(sprintf("  RSI mean: %.1f, median: %.1f\n", 
                  mean(rsi_valid), median(rsi_valid)))
      cat(sprintf("  Times RSI < %d (oversold): %d (%.1f%%)\n", 
                  RSI_OVERSOLD,
                  sum(rsi_valid < RSI_OVERSOLD),
                  100 * mean(rsi_valid < RSI_OVERSOLD)))
      cat(sprintf("  Times RSI > %d (overbought): %d (%.1f%%)\n", 
                  RSI_OVERBOUGHT,
                  sum(rsi_valid > RSI_OVERBOUGHT),
                  100 * mean(rsi_valid > RSI_OVERBOUGHT)))
      cat(sprintf("  Number of positions taken: %d\n", 
                  sum(df_period$position != 0)))
    }
    cat("\n")
  }
}

# Show test period details
cat("TEST PERIOD DETAILED VIEW:\n")
cat(strrep("-", 80), "\n")
df_test <- df %>% 
  filter(period == "Test") %>%
  select(date, close, log_return, rsi, position, strategy_return) %>%
  as.data.frame()

if (nrow(df_test) > 0) {
  df_test$log_return <- sprintf("%.4f", df_test$log_return)
  df_test$rsi <- sprintf("%.1f", df_test$rsi)
  df_test$strategy_return <- sprintf("%.4f", df_test$strategy_return)
  print(df_test, row.names = FALSE)
} else {
  cat("No data in test period\n")
}
cat("\n")

################################################################################
#                     PERFORMANCE ANALYSIS                                      #
################################################################################

cat(strrep("=", 80), "\n")
cat("PERFORMANCE SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat("RSI Strategy Parameters:\n")
cat(sprintf("  RSI Period: %d days\n", RSI_PERIOD))
cat(sprintf("  Oversold Threshold: %d (buy signal)\n", RSI_OVERSOLD))
cat(sprintf("  Overbought Threshold: %d (sell signal)\n", RSI_OVERBOUGHT))
cat(sprintf("  Exit Thresholds: %d (both directions)\n", RSI_EXIT_LOW))
cat(sprintf("  Transaction Cost: %.2f bps per trade\n\n", TRANSACTION_COST * 10000))

periods <- c("Training", "Validation", "Test")
metrics_list <- list()

for (period in periods) {
  df_period <- df %>% filter(period == !!period)
  
  if (nrow(df_period) > 0) {
    metrics <- calc_performance_metrics(
      df_period$strategy_return,
      df_period$position_lagged,
      period
    )
    metrics_list[[period]] <- metrics
  }
}

metrics_df <- bind_rows(metrics_list)

metrics_formatted <- metrics_df %>%
  mutate(
    Annual_Return = sprintf("%.2f%%", Annual_Return * 100),
    Annual_Volatility = sprintf("%.2f%%", Annual_Volatility * 100),
    Sharpe_Ratio = sprintf("%.3f", Sharpe_Ratio),
    Max_Drawdown = sprintf("%.2f%%", Max_Drawdown * 100),
    Win_Rate = ifelse(is.na(Win_Rate), "NA", sprintf("%.2f%%", Win_Rate * 100)),
    Total_Return = sprintf("%.2f%%", Total_Return * 100),
    Avg_Holding_Period = ifelse(is.na(Avg_Holding_Period), "NA", 
                                sprintf("%.1f days", Avg_Holding_Period))
  )

print(metrics_formatted, row.names = FALSE)



################################################################################
#                    PARAMETER OPTIMIZATION                                     #
################################################################################

cat(strrep("=", 80), "\n")
cat("RSI PARAMETER OPTIMIZATION ON VALIDATION PERIOD\n")
cat(strrep("=", 80), "\n\n")

n_combinations <- length(RSI_PERIOD_OPTIONS) * 
  length(OVERSOLD_OPTIONS) * 
  length(OVERBOUGHT_OPTIONS)
cat("Testing", n_combinations, "parameter combinations...\n\n")

optimization_results <- expand.grid(
  rsi_period = RSI_PERIOD_OPTIONS,
  oversold = OVERSOLD_OPTIONS,
  overbought = OVERBOUGHT_OPTIONS,
  stringsAsFactors = FALSE
) %>%
  filter(overbought > oversold + 20) %>%  # Ensure reasonable spread
  mutate(sharpe = NA_real_, max_dd = NA_real_, 
         win_rate = NA_real_, num_trades = NA_integer_)

for (i in 1:nrow(optimization_results)) {
  rsi_p <- optimization_results$rsi_period[i]
  os <- optimization_results$oversold[i]
  ob <- optimization_results$overbought[i]
  
  df_opt <- df_raw %>%
    arrange(assessDate) %>%
    rename(date = assessDate, close = c) %>%
    mutate(
      log_return = c(NA, diff(log(close))),
      rsi = calculate_rsi(close, rsi_p)
    ) %>%
    filter(!is.na(log_return), 
           date >= VAL_START_DATE, 
           date <= VAL_END_DATE)
  
  current_pos <- 0
  positions <- numeric(nrow(df_opt))
  
  for (j in 1:nrow(df_opt)) {
    rsi_val <- df_opt$rsi[j]
    if (is.na(rsi_val)) {
      positions[j] <- 0
      current_pos <- 0
      next
    }
    
    if (current_pos == 0) {
      if (rsi_val < os) current_pos <- 1
      else if (rsi_val > ob) current_pos <- -1
    } else if (current_pos == 1 && rsi_val > 50) {
      current_pos <- 0
    } else if (current_pos == -1 && rsi_val < 50) {
      current_pos <- 0
    }
    
    positions[j] <- current_pos
  }
  
  df_opt$position <- positions
  df_opt$position_lagged <- lag(positions, 1, default = 0)
  df_opt$position_change <- abs(df_opt$position_lagged - 
                                  lag(df_opt$position_lagged, 1, default = 0))
  df_opt$strategy_return <- df_opt$position_lagged * df_opt$log_return - 
    TRANSACTION_COST * df_opt$position_change
  
  optimization_results$sharpe[i] <- calc_sharpe(df_opt$strategy_return, 252)
  optimization_results$max_dd[i] <- calc_max_drawdown(df_opt$strategy_return)
  optimization_results$win_rate[i] <- calc_win_rate(df_opt$strategy_return)
  optimization_results$num_trades[i] <- count_trades(df_opt$position_lagged)
}

optimization_results <- optimization_results %>%
  arrange(desc(sharpe)) %>%
  mutate(
    sharpe_fmt = sprintf("%.3f", sharpe),
    max_dd_fmt = sprintf("%.2f%%", max_dd * 100),
    win_rate_fmt = ifelse(is.na(win_rate), "NA", sprintf("%.2f%%", win_rate * 100))
  )

cat("Top 10 parameter combinations (ranked by Sharpe ratio on validation):\n\n")
print(optimization_results %>% 
        select(rsi_period, oversold, overbought, sharpe_fmt, max_dd_fmt, 
               win_rate_fmt, num_trades) %>%
        rename(RSI_Period = rsi_period, Oversold = oversold, Overbought = overbought,
               Sharpe = sharpe_fmt, MaxDD = max_dd_fmt, 
               WinRate = win_rate_fmt, Trades = num_trades) %>%
        head(10), 
      row.names = FALSE)

################################################################################
#                        EXPORT RESULTS                                         #
################################################################################

trade_log <- df %>%
  filter(period == "Test", position_lagged != 0) %>%
  select(date, close, rsi, position_lagged, log_return, strategy_return, cum_strategy)

write.csv(trade_log, "rsi_trade_log_test_period.csv", row.names = FALSE)

write.csv(df %>% select(date, close, log_return, rsi, position, 
                        position_lagged, strategy_return, cum_strategy, 
                        drawdown, period),
          "rsi_strategy_full_results.csv", row.names = FALSE)

cat("\n\nResults exported:\n")
cat("  - rsi_trade_log_test_period.csv\n")
cat("  - rsi_strategy_full_results.csv\n")
cat("  - 4 plot PNG files\n\n")

cat(strrep("=", 80), "\n")
cat("RSI STRATEGY ANALYSIS COMPLETE\n")
cat(strrep("=", 80), "\n")
cat("\nRECOMMENDATION: Review training period performance to assess strategy validity.\n")
cat("If training Sharpe is positive, this strategy may be viable.\n")
cat("If training Sharpe is negative, consider alternative parameters or approaches.\n")