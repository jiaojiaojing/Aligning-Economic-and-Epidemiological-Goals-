# ===============================
# 4-panel HFC comparison plot
# Panels:
#   A: Infected Humans (HFC, r = 0.02)  + inline table
#   B: Infected Humans (HFC, r = 0.037) + inline table
#   C: Accumulated Deaths (HFC, r = 0.02)
#   D: Accumulated Deaths (HFC, r = 0.037)
# Colors by dif: -30 (blue), 0 (red), 30 (darkgreen)
# ===============================

# ---- Global controls ----
YEARS <- 1
harvest_start <- 213
harvest_duration <- 14
CD <- 30               # control duration
CA <- 30               # yearly control cap (days)
epsilon <- 0.5         # control intensity
gamma <- 0.05
alpha <- 100 * gamma   # health-impact weight on deaths
Omega_default <- 1     # default DFC threshold (not directly used here but kept for parity)

# ===============================
# Simulation for a given dif (UPDATED to allow r_val and Omega_thr)
# ===============================
run_simulation_with_dif <- function(dif_value, r_val = 0.037, Omega_thr = Omega_default) {
  # Parameters (tweak as needed)
  b_H     <- (30/1000)/365
  beta_H  <- 0.000015
  mu_H    <- (8.6/1000)/365
  r       <- r_val                        # <—— recovery rate is now configurable
  eta_M   <- 2
  beta_M  <- 0.000030
  mu_M    <- 1/13
  delta_H <- 0.000055
  years   <- YEARS
  duration <- 365 * years
  
  # Initialize arrays
  B_SH <- B_IH <- B_RH <- B_DH <- B_new_infH <- numeric(duration + 1)
  B_SM <- B_IM <- B_NM <- numeric(duration + 1)
  
  D_SH <- D_IH <- D_RH <- D_DH <- D_new_infH <- numeric(duration + 1)
  D_SM <- D_IM <- D_NM <- numeric(duration + 1)
  
  H_SH <- H_IH <- H_RH <- H_DH <- H_new_infH <- numeric(duration + 1)
  H_SM <- H_IM <- H_NM <- numeric(duration + 1)
  
  # Control variables
  DFC <- DFC_count <- HFC <- HFC_count <- seasonality <- numeric(duration + 1)
  
  # Initial conditions
  B_SH[1] <- 2000; B_IH[1] <- 1; B_RH[1] <- 0; B_DH[1] <- 0; B_new_infH[1] <- 0
  B_SM[1] <- 1000; B_IM[1] <- 0; B_NM[1] <- 1000
  
  D_SH[1] <- 2000; D_IH[1] <- 1; D_RH[1] <- 0; D_DH[1] <- 0; D_new_infH[1] <- 0
  D_SM[1] <- 1000; D_IM[1] <- 0; D_NM[1] <- 1000
  
  H_SH[1] <- 2000; H_IH[1] <- 1; H_RH[1] <- 0; H_DH[1] <- 0; H_new_infH[1] <- 0
  H_SM[1] <- 1000; H_IM[1] <- 0; H_NM[1] <- 1000
  
  DFC[1] <- 0; DFC_count[1] <- 0; HFC[1] <- 0; HFC_count[1] <- 0
  seasonality[1] <- 0.5 * cos(((2 * pi)/365) * (1 + 150 + dif_value)) + 0.5
  
  current_exceeds_threshold_prev <- 0
  total_control_applied_this_year <- 0
  current_year <- 1
  
  for (t in 1:duration) {
    if (t > 1 && (t-1) %% 365 == 0) {
      current_year <- current_year + 1
      total_control_applied_this_year <- 0
    }
    
    # Seasonality with dif
    seasonality[t + 1] <- 0.5 * cos(((2 * pi)/365) * ((t + 1) + 150 + dif_value)) + 0.5
    
    # ---------- Baseline ----------
    B_SH[t + 1] <- max(0, B_SH[t] + (B_SH[t] + B_IH[t] + B_RH[t]) * b_H -
                         min(B_SH[t], beta_H * B_IM[t] * B_SH[t]) - mu_H * B_SH[t])
    B_IH[t + 1] <- max(0, B_IH[t] + min(B_SH[t], beta_H * B_IM[t] * B_SH[t]) -
                         r * B_IH[t] - mu_H * B_IH[t] - delta_H * B_IH[t])
    B_RH[t + 1] <- max(0, B_RH[t] + min(B_IH[t], r * B_IH[t]) - mu_H * B_RH[t])
    B_DH[t + 1] <- B_DH[t] + delta_H * B_IH[t]
    B_new_infH[t + 1] <- min(B_SH[t], beta_H * B_IM[t] * B_SH[t])
    
    B_new_mosquitoes <- (B_IH[1] + B_SH[1]) * eta_M * seasonality[t]
    B_SM[t + 1] <- max(0, B_SM[t] + B_new_mosquitoes -
                         min(B_SM[t], beta_M * B_IH[t] * B_SM[t]) -
                         (mu_M * (1 - seasonality[t])) * B_SM[t])
    B_IM[t + 1] <- max(0, B_IM[t] + min(B_SM[t], beta_M * B_IH[t] * B_SM[t]) -
                         (mu_M * (1 - seasonality[t])) * B_IM[t])
    B_NM[t + 1] <- B_SM[t + 1] + B_IM[t + 1]
    
    # ---------- Disease-Focused Control (DFC) ----------
    D_SH[t + 1] <- max(0, D_SH[t] + (D_SH[t] + D_IH[t] + D_RH[t]) * b_H -
                         min(D_SH[t], beta_H * D_IM[t] * D_SH[t]) - mu_H * D_SH[t])
    D_IH[t + 1] <- max(0, D_IH[t] + min(D_SH[t], beta_H * D_IM[t] * D_SH[t]) -
                         r * D_IH[t] - mu_H * D_IH[t] - delta_H * D_IH[t])
    D_RH[t + 1] <- max(0, D_RH[t] + min(D_IH[t], r * D_IH[t]) - mu_H * D_RH[t])
    D_DH[t + 1] <- D_DH[t] + delta_H * D_IH[t]
    D_new_infH[t + 1] <- min(D_SH[t], beta_H * D_IM[t] * D_SH[t])
    
    D_new_mosquitoes <- (D_IH[1] + D_SH[1]) * eta_M * seasonality[t]
    D_SM[t + 1] <- max(0, D_SM[t] + D_new_mosquitoes -
                         min(D_SM[t], beta_M * D_IH[t] * D_SM[t]) -
                         (mu_M * (1 - seasonality[t])) * D_SM[t] - DFC[t] * D_SM[t])
    D_IM[t + 1] <- max(0, D_IM[t] + min(D_SM[t], beta_M * D_IH[t] * D_SM[t]) -
                         (mu_M * (1 - seasonality[t])) * D_IM[t] - DFC[t] * D_IM[t])
    D_NM[t + 1] <- D_SM[t + 1] + D_IM[t + 1]
    
    # ---------- Harvest-Focused Control (HFC) ----------
    H_SH[t + 1] <- max(0, H_SH[t] + (H_SH[t] + H_IH[t] + H_RH[t]) * b_H -
                         min(H_SH[t], beta_H * H_IM[t] * H_SH[t]) - mu_H * H_SH[t])
    H_IH[t + 1] <- max(0, H_IH[t] + min(H_SH[t], beta_H * H_IM[t] * H_SH[t]) -
                         r * H_IH[t] - mu_H * H_IH[t] - delta_H * H_IH[t])
    H_RH[t + 1] <- max(0, H_RH[t] + min(H_IH[t], r * H_IH[t]) - mu_H * H_RH[t])
    H_DH[t + 1] <- H_DH[t] + delta_H * H_IH[t]
    H_new_infH[t + 1] <- min(H_SH[t], beta_H * H_IM[t] * H_SH[t])
    
    H_new_mosquitoes <- (H_IH[1] + H_SH[1]) * eta_M * seasonality[t]
    H_SM[t + 1] <- max(0, H_SM[t] + H_new_mosquitoes -
                         min(H_SM[t], beta_M * H_IH[t] * H_SM[t]) -
                         (mu_M * (1 - seasonality[t])) * H_SM[t] - HFC[t] * H_SM[t])
    H_IM[t + 1] <- max(0, H_IM[t] + min(H_SM[t], beta_M * H_IH[t] * H_SM[t]) -
                         (mu_M * (1 - seasonality[t])) * H_IM[t] - HFC[t] * H_IM[t])
    H_NM[t + 1] <- H_SM[t + 1] + H_IM[t + 1]
    
    # ---------- Control logic ----------
    # DFC trigger (uses Omega_thr, though DFC results aren't plotted here)
    current_health_impact <- gamma * D_IH[t] + alpha * D_DH[t]
    current_exceeds_threshold <- as.numeric(current_health_impact > Omega_thr)
    previous_exceeds_threshold <- if (t == 1) 0 else current_exceeds_threshold_prev
    
    # DFC trigger + duration + cap
    if (DFC_count[t] == 0) {
      if (current_exceeds_threshold > 0 && previous_exceeds_threshold <= 0 &&
          total_control_applied_this_year < CA) {
        DFC_count[t + 1] <- 1
      } else {
        DFC_count[t + 1] <- 0
      }
    } else {
      if (DFC_count[t] < CD && total_control_applied_this_year < CA) {
        DFC_count[t + 1] <- DFC_count[t] + 1
      } else {
        DFC_count[t + 1] <- 0
      }
    }
    
    if (DFC_count[t + 1] > 0 && DFC_count[t + 1] <= CD &&
        total_control_applied_this_year < CA) {
      DFC[t + 1] <- epsilon
      total_control_applied_this_year <- total_control_applied_this_year + 1
    } else {
      DFC[t + 1] <- 0
    }
    current_exceeds_threshold_prev <- current_exceeds_threshold
    
    # HFC schedule (CD days ending right before harvest start)
    current_harvest_start <- harvest_start + 365 * (current_year - 1)
    if (HFC_count[t] == 0) {
      if (t + CD == current_harvest_start) {
        HFC_count[t + 1] <- 1
      } else {
        HFC_count[t + 1] <- 0
      }
    } else {
      if (HFC_count[t] < CD) {
        HFC_count[t + 1] <- HFC_count[t] + 1
      } else {
        HFC_count[t + 1] <- 0
      }
    }
    HFC[t + 1] <- if (HFC_count[t + 1] > 0 && HFC_count[t + 1] <= CD) epsilon else 0
  }
  
  list(
    dif = dif_value,
    baseline = list(
      IH_timeseries    = B_IH,
      DH_timeseries    = B_DH,
      NM_timeseries    = B_NM,
      total_infections = sum(B_IH),
      peak_infections  = max(B_IH),
      total_deaths     = B_DH[length(B_DH)],
      final_infected   = B_IH[length(B_IH)],
      new_infH_timeseries = B_new_infH
    ),
    disease_control = list(
      IH_timeseries    = D_IH,
      DH_timeseries    = D_DH,
      NM_timeseries    = D_NM,
      DFC_timeseries   = DFC,
      total_infections = sum(D_IH),
      peak_infections  = max(D_IH),
      total_deaths     = D_DH[length(D_DH)],
      final_infected   = D_IH[length(D_IH)],
      new_infH_timeseries = D_new_infH,
      control_days     = sum(DFC > 0)
    ),
    harvest_control = list(
      IH_timeseries    = H_IH,
      DH_timeseries    = H_DH,
      NM_timeseries    = H_NM,
      HFC_timeseries   = HFC,
      total_infections = sum(H_IH),
      peak_infections  = max(H_IH),
      total_deaths     = H_DH[length(H_DH)],
      final_infected   = H_IH[length(H_IH)],
      new_infH_timeseries = H_new_infH,
      control_days     = sum(HFC > 0)
    ),
    seasonality = seasonality
  )
}

# ===============================
# Helpers
# ===============================
# IH person-days during the harvest window (safe across years)
calculate_harvest_infections <- function(IH_timeseries,
                                         years = YEARS,
                                         hs = harvest_start,
                                         hd = harvest_duration) {
  total <- 0
  for (yr in 1:years) {
    start_day <- hs + 365 * (yr - 1)
    end_day   <- start_day + hd - 1
    idx <- start_day:end_day
    idx <- idx[idx >= 1 & idx <= length(IH_timeseries)]
    if (length(idx) > 0) total <- total + sum(IH_timeseries[idx], na.rm = TRUE)
  }
  total
}

# Harvest shading helper (for base R panels)
shade_harvest <- function() {
  usr <- par("usr")
  rect(harvest_start, usr[3],
       harvest_start + harvest_duration - 1, usr[4],
       col = rgb(0.3, 0.3, 0.3, 0.2), border = NA)
}

# ===============================
# Run HFC packs and build the 4-panel plot
# ===============================
dif_levels   <- c(-30, 0, 30)
curve_colors <- c("blue", "red", "darkgreen")
line_types   <- c(1, 1, 1)              # all solid
time_days    <- 1:(365 * YEARS)

# Helper to run HFC-only sims at a given r and compute table/ylims
run_hfc_pack <- function(r_here) {
  res_list <- list()
  table_df <- data.frame()
  for (dv in dif_levels) {
    sim <- run_simulation_with_dif(dv, r_val = r_here, Omega_thr = Omega_default)
    base_total <- sim$baseline$total_infections
    hfc        <- sim$harvest_control
    
    h_pd    <- calculate_harvest_infections(hfc$IH_timeseries, years = YEARS,
                                            hs = harvest_start, hd = harvest_duration)
    tot_IH  <- hfc$total_infections
    effort  <- sum(hfc$HFC_timeseries, na.rm = TRUE)
    eff     <- if (effort > 0) (base_total - tot_IH) / effort else NA_real_
    
    table_df <- rbind(table_df, data.frame(
      dif = dv,
      harvest_infections = round(h_pd, 1),
      total_infections   = round(tot_IH, 1),
      efficacy           = round(eff, 3)
    ))
    res_list[[as.character(dv)]] <- hfc
  }
  ymax_IH <- max(sapply(res_list, function(x) max(x$IH_timeseries[time_days], na.rm = TRUE)), na.rm = TRUE)
  ymax_DH <- max(sapply(res_list, function(x) max(x$DH_timeseries[time_days], na.rm = TRUE)), na.rm = TRUE)
  if (!is.finite(ymax_IH) || ymax_IH <= 0) ymax_IH <- 1
  if (!is.finite(ymax_DH) || ymax_DH <= 0) ymax_DH <- 1
  list(series = res_list, table = table_df, ymax_IH = ymax_IH, ymax_DH = ymax_DH)
}

pack_r002  <- run_hfc_pack(0.02)
pack_r0037 <- run_hfc_pack(0.037)

# Use common y-lims across rows for clean comparison
ymax_IH_both <- max(pack_r002$ymax_IH, pack_r0037$ymax_IH)
ymax_DH_both <- max(pack_r002$ymax_DH, pack_r0037$ymax_DH)

# ---- Plotting ----
png("HFC_r_comparison_dif_-30_0_30_4panel.png",
    width = 11, height = 8.5, units = "in", res = 300)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))

# Panel A: Infected Humans, r = 0.02
plot(1, type = "n", xlim = c(1, 365), ylim = c(0, ymax_IH_both),
     xlab = "Time (days)", ylab = "Infected Humans",
     main = "Infected Humans (r = 0.02)")
mtext("A", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  lines(time_days, pack_r002$series[[key]]$IH_timeseries[time_days],
        col = curve_colors[i], lty = line_types[i], lwd = 2)
}
# Inline table (colored rows)
legend_x <- 365 * 0.53
legend_y <- ymax_IH_both * 0.95
line_h   <- ymax_IH_both * 0.06
text(legend_x, legend_y, "dif   Harvest       Total        Eff", font = 2, cex = 0.95, pos = 4)
for (i in seq_len(nrow(pack_r002$table))) {
  ry <- legend_y - i * line_h
  txt <- sprintf("%3d %12.1f %12.1f %9.3f",
                 pack_r002$table$dif[i],
                 pack_r002$table$harvest_infections[i],
                 pack_r002$table$total_infections[i],
                 pack_r002$table$efficacy[i])
  text(legend_x, ry, txt, col = curve_colors[i], cex = 0.9, pos = 4)
}

# Panel B: Infected Humans, r = 0.037
plot(1, type = "n", xlim = c(1, 365), ylim = c(0, ymax_IH_both),
     xlab = "Time (days)", ylab = "Infected Humans",
     main = "Infected Humans (r = 0.037)")
mtext("B", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  lines(time_days, pack_r0037$series[[key]]$IH_timeseries[time_days],
        col = curve_colors[i], lty = line_types[i], lwd = 2)
}
# Inline table (colored rows)
legend_x <- 365 * 0.53
legend_y <- ymax_IH_both * 0.95
line_h   <- ymax_IH_both * 0.06
text(legend_x, legend_y, "dif   Harvest       Total        Eff", font = 2, cex = 0.95, pos = 4)
for (i in seq_len(nrow(pack_r0037$table))) {
  ry <- legend_y - i * line_h
  txt <- sprintf("%3d %12.1f %12.1f %9.3f",
                 pack_r0037$table$dif[i],
                 pack_r0037$table$harvest_infections[i],
                 pack_r0037$table$total_infections[i],
                 pack_r0037$table$efficacy[i])
  text(legend_x, ry, txt, col = curve_colors[i], cex = 0.9, pos = 4)
}

# Panel C: Accumulated Deaths, r = 0.02
plot(1, type = "n", xlim = c(1, 365), ylim = c(0, ymax_DH_both),
     xlab = "Time (days)", ylab = "Accumulated Deaths",
     main = "Accumulated Deaths (r = 0.02)")
mtext("C", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  lines(time_days, pack_r002$series[[key]]$DH_timeseries[time_days],
        col = curve_colors[i], lty = line_types[i], lwd = 2)
}

# Panel D: Accumulated Deaths, r = 0.037
plot(1, type = "n", xlim = c(1, 365), ylim = c(0, ymax_DH_both),
     xlab = "Time (days)", ylab = "Accumulated Deaths",
     main = "Accumulated Deaths (r = 0.037)")
mtext("D", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  lines(time_days, pack_r0037$series[[key]]$DH_timeseries[time_days],
        col = curve_colors[i], lty = line_types[i], lwd = 2)
}

mtext("Harvest-Focused Control (HFC): Recovery-Rate Comparison across Seasonal Shifts",
      outer = TRUE, cex = 1.2, font = 2)

dev.off()

cat("\nSaved: HFC_r_comparison_dif_-30_0_30_4panel-new.png\n")
