# ===============================
# Baseline-only seasonal-shift figure (2x2)
# Panels: SH, IH(+table), DH, NM, all with harvest shading
# ===============================

# ---- (optional) set working directory for outputs ----
# setwd("C:/path/to/output")

# ---- Global parameters ----
YEARS <- 1
harvest_start <- 213
harvest_duration <- 14
CD <- 30 # unused here (DFC/HFC not simulated), kept for parity

# ---- Simulation (Baseline only) ----
run_baseline_with_dif <- function(dif_value, 
                                  years = YEARS,
                                  b_H = (30/1000)/365,
                                  beta_H = 0.000015,
                                  mu_H = (8.6/1000)/365,
                                  r = 0.037,
                                  eta_M = 2,
                                  beta_M = 0.000030,
                                  mu_M = 1/13,
                                  delta_H = 0.000055) {
  
  duration <- 365 * years
  
  # time series
  SH <- IH <- RH <- DH <- new_infH <- numeric(duration + 1)
  SM <- IM <- NM <- numeric(duration + 1)
  seasonality <- numeric(duration + 1)
  
  # initials (match your baseline block)
  SH[1] <- 2000; IH[1] <- 1; RH[1] <- 0; DH[1] <- 0; new_infH[1] <- 0
  SM[1] <- 1000; IM[1] <- 0; NM[1] <- 1000
  
  # seasonality
  seasonality[1] <- 0.5 * cos(((2 * pi)/365) * (1 + 150 + dif_value)) + 0.5
  
  for (t in 1:duration) {
    # seasonality (shifted by dif)
    seasonality[t + 1] <- 0.5 * cos(((2 * pi)/365) * ((t + 1) + 150 + dif_value)) + 0.5
    
    # humans (baseline)
    inf_flow <- min(SH[t], beta_H * IM[t] * SH[t]) # new infections today
    rec_flow <- r * IH[t]
    
    SH[t + 1] <- max(0, SH[t] + (SH[t] + IH[t] + RH[t]) * b_H - inf_flow - mu_H * SH[t])
    IH[t + 1] <- max(0, IH[t] + inf_flow - rec_flow - mu_H * IH[t] - delta_H * IH[t])
    RH[t + 1] <- max(0, RH[t] + rec_flow - mu_H * RH[t])
    DH[t + 1] <- DH[t] + delta_H * IH[t]
    new_infH[t + 1] <- inf_flow
    
    # mosquitoes (baseline)
    # NOTE: births proportional to **current** SH+IH for stronger coupling.
    # If you want constant (initial) coupling, switch SH[t]+IH[t] -> SH[1]+IH[1]
    new_mosq <- (IH[t] + SH[t]) * eta_M * seasonality[t]
    inf_m_flow <- min(SM[t], beta_M * IH[t] * SM[t])
    death_term <- (mu_M * (1 - seasonality[t]))
    
    SM[t + 1] <- max(0, SM[t] + new_mosq - inf_m_flow - death_term * SM[t])
    IM[t + 1] <- max(0, IM[t] + inf_m_flow - death_term * IM[t])
    NM[t + 1] <- SM[t + 1] + IM[t + 1]
  }
  
  list(
    dif = dif_value,
    SH_timeseries = SH,
    IH_timeseries = IH,
    RH_timeseries = RH,
    DH_timeseries = DH,
    NM_timeseries = NM,
    new_infH_timeseries = new_infH,
    # "overall human infection" as IH person-days (prevalence, matches your DFC/HFC tables)
    total_IH_persondays = sum(IH, na.rm = TRUE),
    # Optionally, also expose total incidence (cases) if needed:
    total_incidence = sum(new_infH, na.rm = TRUE),
    seasonality = seasonality
  )
}

# ---- Harvest-window helpers ----
harvest_IH_persondays <- function(IH, hs = harvest_start, hd = harvest_duration) {
  s <- hs; e <- hs + hd - 1
  idx <- s:e
  idx <- idx[idx >= 1 & idx <= length(IH)]
  if (length(idx) == 0) return(NA_real_)
  sum(IH[idx], na.rm = TRUE)
}

# If you prefer cases instead of person-days:
# harvest_cases <- function(new_inf, hs = harvest_start, hd = harvest_duration) {
#   s <- hs; e <- hs + hd - 1
#   idx <- s:e
#   idx <- idx[idx >= 1 & idx <= length(new_inf)]
#   if (length(idx) == 0) return(NA_real_)
#   sum(new_inf[idx], na.rm = TRUE)
# }

# ---- Robust plotting helpers ----
safe_series <- function(x, idx) {
  if (is.null(x)) return(numeric(0))
  n <- length(x)
  if (n == 0) return(numeric(0))
  idx2 <- idx[idx >= 1 & idx <= n]
  y <- suppressWarnings(as.numeric(x[idx2]))
  y[is.finite(y)]
}

get_ylim <- function(results_list, field, idx) {
  y <- unlist(lapply(names(results_list), function(k) {
    res <- results_list[[k]]
    safe_series(res[[field]], idx)
  }))
  if (length(y) == 0) c(0, 1) else c(0, max(y))
}

shade_harvest <- function() {
  usr <- par("usr")
  rect(harvest_start, usr[3], harvest_start + harvest_duration - 1, usr[4], 
       col = rgb(0.3, 0.3, 0.3, 0.2), border = NA)
}

# ---- Sims to compare ----
dif_levels <- c(-30, 0, 30)
colors <- c("blue", "red", "darkgreen")
line_types <- c(1, 1, 1) # solid, dashed, dotted
time_days <- 1:(365 * YEARS)

baseline_results <- setNames(vector("list", length(dif_levels)), as.character(dif_levels))
for (dif_val in dif_levels) {
  cat("Running baseline for dif =", dif_val, "...\n")
  baseline_results[[as.character(dif_val)]] <- run_baseline_with_dif(dif_val)
}

# ---- Build inline legend table data (IH person-days) ----
table_data_base <- do.call(rbind, lapply(dif_levels, function(dv) {
  base <- baseline_results[[as.character(dv)]]
  data.frame(
    dif = dv,
    overall_infections = round(base$total_IH_persondays, 1), # IH person-days
    harvest_infections = round(harvest_IH_persondays(base$IH_timeseries), 1) # IH person-days in harvest
    # If using cases instead:
    # harvest_infections = round(harvest_cases(base$new_infH_timeseries), 1)
  )
}))

print(table_data_base)

# ---- Determine y-limits robustly ----
ylim_SH <- get_ylim(baseline_results, "SH_timeseries", time_days)
ylim_IH <- get_ylim(baseline_results, "IH_timeseries", time_days)
ylim_DH <- get_ylim(baseline_results, "DH_timeseries", time_days)
ylim_NM <- get_ylim(baseline_results, "NM_timeseries", time_days)

# ---- Plot: 2x2 with harvest shading and inline table on IH panel ----
png("baseline_seasonal_shift_comparison_2x2-10-22-25.png", 
    width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))

# Panel 1: Susceptible Humans (SH)
plot(1, type = "n", xlim = c(1, max(time_days)), ylim = ylim_SH,
     xlab = "Time (days)", ylab = "Susceptible Humans", 
     main = "Susceptible Humans")
mtext("A", side=3, adj=0, line=0.5, font=2, cex=1.2) # <— label top-left
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  y <- safe_series(baseline_results[[key]]$SH_timeseries, time_days)
  x <- seq_along(y)
  if (length(y)) lines(x, y, col = colors[i], lty = line_types[i], lwd = 2)
}

# Panel 2: Infected Humans (IH) + inline table
plot(1, type = "n", xlim = c(1, max(time_days)), ylim = ylim_IH,
     xlab = "Time (days)", ylab = "Infected Humans", 
     main = "Infected Humans")
mtext("B", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  y <- safe_series(baseline_results[[key]]$IH_timeseries, time_days)
  x <- seq_along(y)
  if (length(y)) lines(x, y, col = colors[i], lty = line_types[i], lwd = 2)
}

# Inline table (dif, Overall, Harvest)
legend_x <- max(time_days) * 0.55
legend_y <- 0.95 * max(1, ylim_IH[2])
line_h <- 0.06 * max(1, ylim_IH[2])

text(legend_x, legend_y, "dif Overall Harvest", font = 2, cex = 0.95, pos = 4)
for (i in seq_len(nrow(table_data_base))) {
  ry <- legend_y - i * line_h
  text(legend_x, ry, 
       sprintf("%3d %10.1f %11.1f", 
               table_data_base$dif[i], 
               table_data_base$overall_infections[i], 
               table_data_base$harvest_infections[i]),
       col = colors[i], cex = 0.9, pos = 4)
}

# Panel 3: Accumulated Deaths (DH)
plot(1, type = "n", xlim = c(1, max(time_days)), ylim = ylim_DH,
     xlab = "Time (days)", ylab = "Accumulated Deaths", 
     main = "Accumulated Deaths")
mtext("C", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  y <- safe_series(baseline_results[[key]]$DH_timeseries, time_days)
  x <- seq_along(y)
  if (length(y)) lines(x, y, col = colors[i], lty = line_types[i], lwd = 2)
}

# Panel 4: Mosquito Total Population (NM)
plot(1, type = "n", xlim = c(1, max(time_days)), ylim = ylim_NM,
     xlab = "Time (days)", ylab = "Total Mosquito Population", 
     main = "Mosquito Population")
mtext("D", side=3, adj=0, line=0.5, font=2, cex=1.2)
shade_harvest()
for (i in seq_along(dif_levels)) {
  key <- as.character(dif_levels[i])
  y <- safe_series(baseline_results[[key]]$NM_timeseries, time_days)
  x <- seq_along(y)
  if (length(y)) lines(x, y, col = colors[i], lty = line_types[i], lwd = 2)
}

mtext("Baseline - Seasonal Shift Comparison", outer = TRUE, cex = 1.2, font = 2)
dev.off()
