# ===============================================================
# Vector-borne model with Omega sensitivity (Baseline + DFC plots)
# ===============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:\\Users\\13522\\Downloads\\Crop ms")

# -------- global constants ----------
YEARS <- 1
harvest_start <- 213
harvest_duration <- 14
CD <- 30

# ===============================================================
# Simulation function (now accepts Omega as argument)
# ===============================================================
run_simulation_with_dif <- function(dif_value, scenario_name = "all", Omega = 1) {
  b_H     <- (30/1000)/365
  beta_H  <- 0.000015
  mu_H    <- (8.6/1000)/365
  r       <- 0.037
  eta_M   <- 2
  beta_M  <- 0.000030
  mu_M    <- 1/13
  delta_H <- 0.000055
  gamma   <- 0.05
  alpha   <- 100 * gamma
  epsilon <- 0.5
  CD_loc  <- CD
  CA      <- 30
  years   <- YEARS
  duration <- 365 * years
  
  # Arrays
  B_IH <- D_IH <- H_IH <- numeric(duration + 1)
  B_DH <- D_DH <- H_DH <- numeric(duration + 1)
  B_SM <- D_SM <- H_SM <- numeric(duration + 1)
  B_IM <- D_IM <- H_IM <- numeric(duration + 1)
  B_NM <- D_NM <- H_NM <- numeric(duration + 1)
  B_SH <- D_SH <- H_SH <- numeric(duration + 1)
  B_RH <- D_RH <- H_RH <- numeric(duration + 1)
  DFC <- DFC_count <- HFC <- HFC_count <- seasonality <- numeric(duration + 1)
  
  # Initial conditions
  B_SH[1] <- D_SH[1] <- H_SH[1] <- 2000
  B_IH[1] <- D_IH[1] <- H_IH[1] <- 1
  B_SM[1] <- D_SM[1] <- H_SM[1] <- 1000
  B_IM[1] <- D_IM[1] <- H_IM[1] <- 0
  B_NM[1] <- D_NM[1] <- H_NM[1] <- 1000
  
  seasonality[1] <- 0.5 * cos(((2 * pi)/365) * (1 + 150 + dif_value)) + 0.5
  current_exceeds_threshold_prev <- 0
  total_control_applied_this_year <- 0
  current_year <- 1
  
  for (t in 1:duration) {
    if (t > 1 && (t-1) %% 365 == 0) {
      current_year <- current_year + 1
      total_control_applied_this_year <- 0
    }
    
    seasonality[t + 1] <- 0.5 * cos(((2 * pi)/365) * ((t + 1) + 150 + dif_value)) + 0.5
    
    # ---------- Baseline ----------
    B_SH[t + 1] <- max(0, B_SH[t] + (B_SH[t] + B_IH[t] + B_RH[t]) * b_H -
                         min(B_SH[t], beta_H * B_IM[t] * B_SH[t]) - mu_H * B_SH[t])
    B_IH[t + 1] <- max(0, B_IH[t] + min(B_SH[t], beta_H * B_IM[t] * B_SH[t]) -
                         r * B_IH[t] - mu_H * B_IH[t] - delta_H * B_IH[t])
    B_RH[t + 1] <- max(0, B_RH[t] + r * B_IH[t] - mu_H * B_RH[t])
    B_DH[t + 1] <- B_DH[t] + delta_H * B_IH[t]
    
    B_SM[t + 1] <- max(0, B_SM[t] + (B_IH[1] + B_SH[1]) * eta_M * seasonality[t] -
                         beta_M * B_IH[t] * B_SM[t] -
                         (mu_M * (1 - seasonality[t])) * B_SM[t])
    B_IM[t + 1] <- max(0, B_IM[t] + beta_M * B_IH[t] * B_SM[t] -
                         (mu_M * (1 - seasonality[t])) * B_IM[t])
    B_NM[t + 1] <- B_SM[t + 1] + B_IM[t + 1]
    
    # ---------- Disease-focused control ----------
    D_SH[t + 1] <- max(0, D_SH[t] + (D_SH[t] + D_IH[t] + D_RH[t]) * b_H -
                         beta_H * D_IM[t] * D_SH[t] - mu_H * D_SH[t])
    D_IH[t + 1] <- max(0, D_IH[t] + beta_H * D_IM[t] * D_SH[t] -
                         r * D_IH[t] - mu_H * D_IH[t] - delta_H * D_IH[t])
    D_RH[t + 1] <- max(0, D_RH[t] + r * D_IH[t] - mu_H * D_RH[t])
    D_DH[t + 1] <- D_DH[t] + delta_H * D_IH[t]
    
    current_health_impact <- gamma * D_IH[t] + alpha * D_DH[t]
    current_exceeds_threshold <- as.numeric(current_health_impact > Omega)
    previous_exceeds_threshold <- if (t == 1) 0 else current_exceeds_threshold_prev
    
    # DFC control logic
    if (DFC_count[t] == 0) {
      if (current_exceeds_threshold > 0 && previous_exceeds_threshold <= 0 &&
          total_control_applied_this_year < 30) {
        DFC_count[t + 1] <- 1
      } else {
        DFC_count[t + 1] <- 0
      }
    } else {
      if (DFC_count[t] < CD && total_control_applied_this_year < 30) {
        DFC_count[t + 1] <- DFC_count[t] + 1
      } else {
        DFC_count[t + 1] <- 0
      }
    }
    if (DFC_count[t + 1] > 0 && DFC_count[t + 1] <= CD &&
        total_control_applied_this_year < 30) {
      DFC[t + 1] <- epsilon
      total_control_applied_this_year <- total_control_applied_this_year + 1
    } else {
      DFC[t + 1] <- 0
    }
    current_exceeds_threshold_prev <- current_exceeds_threshold
    
    D_SM[t + 1] <- max(0, D_SM[t] + (D_IH[1] + D_SH[1]) * eta_M * seasonality[t] -
                         beta_M * D_IH[t] * D_SM[t] -
                         (mu_M * (1 - seasonality[t])) * D_SM[t] - DFC[t] * D_SM[t])
    D_IM[t + 1] <- max(0, D_IM[t] + beta_M * D_IH[t] * D_SM[t] -
                         (mu_M * (1 - seasonality[t])) * D_IM[t] - DFC[t] * D_IM[t])
    D_NM[t + 1] <- D_SM[t + 1] + D_IM[t + 1]
  }
  
  list(
    baseline = list(IH_timeseries = B_IH, DH_timeseries = B_DH),
    disease_control = list(IH_timeseries = D_IH, DH_timeseries = D_DH,
                           DFC_timeseries = DFC),
    seasonality = seasonality
  )
}

# ===============================================================
# Omega comparison plotting (Ω = 1 vs 50)
# ===============================================================
YEARS <- 1
time_days <- 1:(365 * YEARS)
harvest_start <- 213
harvest_duration <- 14
dif_levels <- c(-30, -15, 0, 30)
colors <- c("blue", "purple", "red", "green")

xy_aligned <- function(x, y) {
  n <- min(length(x), length(y))
  list(x = x[seq_len(n)], y = y[seq_len(n)])
}

shade_harvest <- function() {
  usr <- par("usr")
  rect(harvest_start, usr[3],
       harvest_start + harvest_duration - 1, usr[4],
       col = rgb(0.3, 0.3, 0.3, 0.2), border = NA)
}

calc_harvest_infections <- function(IH) {
  idx <- harvest_start:(harvest_start + harvest_duration - 1)
  idx <- idx[idx <= length(IH)]
  sum(IH[idx], na.rm = TRUE)
}

# Run model sets
run_sims_for_omega <- function(Omega_value) {
  sims <- list()
  for (dif in dif_levels) {
    cat("Running Ω =", Omega_value, "dif =", dif, "\n")
    sims[[as.character(dif)]] <- run_simulation_with_dif(dif, Omega = Omega_value)$disease_control
  }
  sims
}

# Baseline and DFC sims
baseline_sims <- lapply(dif_levels, function(dif) {
  run_simulation_with_dif(dif, Omega = Inf)$baseline
})
names(baseline_sims) <- as.character(dif_levels)
dfc_omega1  <- run_sims_for_omega(1)
dfc_omega50 <- run_sims_for_omega(50)

# Table builder
make_table <- function(sims, baseline_sims) {
  do.call(rbind, lapply(dif_levels, function(dif) {
    key <- as.character(dif)
    IH <- sims[[key]]$IH_timeseries
    harvest_inf <- calc_harvest_infections(IH)
    total_inf <- sum(IH, na.rm = TRUE)
    base_total <- sum(baseline_sims[[key]]$IH_timeseries, na.rm = TRUE)
    effort <- sum(sims[[key]]$DFC_timeseries, na.rm = TRUE)
    eff <- if (effort > 0) (base_total - total_inf) / effort else NA
    data.frame(dif, harvest_inf, total_inf, eff)
  }))
}

table_omega1  <- make_table(dfc_omega1, baseline_sims)
table_omega50 <- make_table(dfc_omega50, baseline_sims)

safe_max <- function(list_sims, field) {
  vals <- sapply(list_sims, function(x) max(x[[field]], na.rm = TRUE))
  max(vals[is.finite(vals)], na.rm = TRUE)
}
max_IH <- max(safe_max(dfc_omega1, "IH_timeseries"),
              safe_max(dfc_omega50, "IH_timeseries"))
max_DH <- max(safe_max(dfc_omega1, "DH_timeseries"),
              safe_max(dfc_omega50, "DH_timeseries"))

# ---- Plot ----
png("DFC_Omega1_vs50_comparison_final.png", width=11, height=8.5, units="in", res=300)
par(mfrow=c(2,2), mar=c(4,4,2,2), oma=c(0,0,2,0))

# Panel A
plot(1, type="n", xlim=c(1,365), ylim=c(0,max_IH),
     xlab="Time (days)", ylab="Infected Humans",
     main=expression(paste("Infected Humans (",Omega,"=1)")))
mtext("A", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  xy <- xy_aligned(time_days, dfc_omega1[[key]]$IH_timeseries)
  lines(xy$x, xy$y, col=colors[i], lwd=2)
}
lx <- 0; ly <- 0.95 * max_IH; step <- 0.07 * max_IH
text(lx, ly, "dif   Harvest   Total   Eff", font=2, pos=4)
for (i in seq_len(nrow(table_omega1))) {
  row <- table_omega1[i,]; ry <- ly - i * step
  txt <- sprintf("%3d %9.1f %8.1f %7.3f",
                 row$dif, row$harvest_inf, row$total_inf, row$eff)
  text(lx, ry, txt, col=colors[i], pos=4)
}

# Panel B
plot(1, type="n", xlim=c(1,365), ylim=c(0,max_IH),
     xlab="Time (days)", ylab="",
     main=expression(paste("Infected Humans (",Omega,"=50)")))
mtext("B", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  xy <- xy_aligned(time_days, dfc_omega50[[key]]$IH_timeseries)
  lines(xy$x, xy$y, col=colors[i], lwd=2)
}
lx1 = 210
text(lx1, ly, "dif   Harvest   Total   Eff", font=2, pos=4)
for (i in seq_len(nrow(table_omega50))) {
  row <- table_omega50[i,]; ry <- ly - i * step
  txt <- sprintf("%3d %9.1f %8.1f %7.3f",
                 row$dif, row$harvest_inf, row$total_inf, row$eff)
  text(lx1, ry, txt, col=colors[i], pos=4)
}

# Panel C
plot(1, type="n", xlim=c(1,365), ylim=c(0,max_DH),
     xlab="Time (days)", ylab="Accumulated Deaths",
     main=expression(paste("Accumulated Deaths (",Omega,"=1)")))
mtext("C", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  xy <- xy_aligned(time_days, dfc_omega1[[key]]$DH_timeseries)
  lines(xy$x, xy$y, col=colors[i], lwd=2)
}

# Panel D
plot(1, type="n", xlim=c(1,365), ylim=c(0,max_DH),
     xlab="Time (days)", ylab="",
     main=expression(paste("Accumulated Deaths (",Omega,"=50)")))
mtext("D", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  xy <- xy_aligned(time_days, dfc_omega50[[key]]$DH_timeseries)
  lines(xy$x, xy$y, col=colors[i], lwd=2)
}

mtext("Disease-Focused Control (DFC): Ω Sensitivity and Seasonal Shift",
      outer=TRUE, cex=1.2, font=2)
dev.off()
