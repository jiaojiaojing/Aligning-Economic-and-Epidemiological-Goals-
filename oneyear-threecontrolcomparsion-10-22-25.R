# ---------------------------
# Vector-borne model plotting (1-year)
# Legend only in Panel B; A–D tags on panels
# ---------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# (Optional) set your working directory
# setwd("C:\\Users\\13522\\Downloads\\Crop ms")

# -------- global constants ----------
YEARS <- 1
harvest_start <- 213
harvest_duration <- 14
CD <- 30                  # control duration (days)
CA <- 30                  # yearly cap for DFC days
epsilon <- 0.5
gamma <- 0.05
alpha <- 100 * gamma

# Choose recovery rate and threshold used for BOTH simulation and filename
R_FOR_PLOTS <- 0.02
OMEGA_FOR_PLOTS <- 1

# -------- plotting palettes ----------
# All-black lines; distinguish scenarios by line type only
scenario_colors <- c(
  "Baseline"        = "black",
  "Disease Control" = "black",
  "Harvest Control" = "black"
)

# Linetype mapping (Baseline=solid, DFC=dashed, HFC=dotted)
line_types_common <- c(
  "Baseline"        = "solid",
  "Disease Control" = "dashed",
  "Harvest Control" = "dotted"
)

# For the p4 two-scenario plot (DFC vs HFC)
p4_colors <- c(
  "Disease Control" = "black",
  "Harvest Control" = "black"
)
p4_linetypes <- c(
  "Disease Control" = "dashed",
  "Harvest Control" = "dotted"
)

# ===============================
# Simulation for a given dif (parameterized r and Omega)
# ===============================
run_simulation_with_dif <- function(dif_value, r_val = R_FOR_PLOTS, Omega_thr = OMEGA_FOR_PLOTS) {
  # Parameters (tweak as needed)
  b_H     <- (30/1000)/365
  beta_H  <- 0.000015
  mu_H    <- (8.6/1000)/365
  r       <- r_val
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
    # DFC trigger uses Omega_thr
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
      IH_timeseries       = B_IH,
      DH_timeseries       = B_DH,
      NM_timeseries       = B_NM,
      total_infections    = sum(B_IH),     # IH person-days (prevalence integral)
      peak_infections     = max(B_IH),
      total_deaths        = B_DH[length(B_DH)],
      final_infected      = B_IH[length(B_IH)],
      new_infH_timeseries = B_new_infH
    ),
    disease_control = list(
      IH_timeseries       = D_IH,
      DH_timeseries       = D_DH,
      NM_timeseries       = D_NM,
      DFC_timeseries      = DFC,
      total_infections    = sum(D_IH),     # IH person-days
      peak_infections     = max(D_IH),
      total_deaths        = D_DH[length(D_DH)],
      final_infected      = D_IH[length(D_IH)],
      new_infH_timeseries = D_new_infH,
      control_days        = sum(DFC > 0)
    ),
    harvest_control = list(
      IH_timeseries       = H_IH,
      DH_timeseries       = H_DH,
      NM_timeseries       = H_NM,
      HFC_timeseries      = HFC,
      total_infections    = sum(H_IH),     # IH person-days
      peak_infections     = max(H_IH),
      total_deaths        = H_DH[length(H_DH)],
      final_infected      = H_IH[length(H_IH)],
      new_infH_timeseries = H_new_infH,
      control_days        = sum(HFC > 0)
    ),
    seasonality = seasonality
  )
}

# ===============================
# Helpers
# ===============================
# IH person-days in harvest window
harvest_IH_pd <- function(IH, hs = harvest_start, hd = harvest_duration) {
  s <- hs; e <- hs + hd - 1
  idx <- s:e
  idx <- idx[idx >= 1 & idx <= length(IH)]
  if (length(idx) == 0) return(NA_real_)
  sum(IH[idx], na.rm = TRUE)
}

# SAFE harvest infections (flexible years + clipping)
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

# ===============================
# Run analysis across dif values
# ===============================
test_dif_effects <- function(dif_values = c(-30, -15, -7, 0, 7, 15, 30),
                             r_here = R_FOR_PLOTS,
                             Omega_here = OMEGA_FOR_PLOTS) {
  results <- vector("list", length(dif_values))
  names(results) <- as.character(dif_values)
  for (dif_val in dif_values) {
    cat("Testing dif =", dif_val, "with r =", r_here, "and Omega =", Omega_here, "...\n")
    results[[as.character(dif_val)]] <- run_simulation_with_dif(dif_val, r_val = r_here, Omega_thr = Omega_here)
  }
  results
}

dif_values_all <- c(-30, -15, -7, 0, 7, 15, 30)
dif_results <- test_dif_effects(dif_values_all)

# ===============================
# Summary table (IH person-days totals, peaks, deaths)
# ===============================
summary_table <- bind_rows(lapply(names(dif_results), function(dif_val) {
  res <- dif_results[[dif_val]]
  tibble(
    dif = as.numeric(dif_val),
    scenario = c("Baseline", "Disease Control", "Harvest Control"),
    total_infections = c(res$baseline$total_infections,
                         res$disease_control$total_infections,
                         res$harvest_control$total_infections),
    peak_infections  = c(res$baseline$peak_infections,
                         res$disease_control$peak_infections,
                         res$harvest_control$peak_infections),
    total_deaths     = c(res$baseline$total_deaths,
                         res$disease_control$total_deaths,
                         res$harvest_control$total_deaths)
  )
}))

# ===============================
# Unified theme for all panels
# ===============================
big_theme <- theme_minimal(base_size = 18) +
  theme(
    axis.title  = element_text(size = 20, face = "bold"),
    axis.text   = element_text(size = 16),
    plot.title  = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    panel.grid  = element_blank(),
    axis.line   = element_line(linewidth = 0.6)
  )

# ===============================
# Plots p1–p4 (ggplot + patchwork) — legend only in Panel B, A–D tags
# ===============================

# Panel A (no legend): Total infections vs dif
p1 <- ggplot(summary_table,
             aes(x = dif, y = total_infections,
                 color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  labs(x = "", y = "Total Infections") +
  scale_color_manual(values = scenario_colors) +
  scale_linetype_manual(values = line_types_common) +
  big_theme +
  theme(legend.position = "none")

# Panel B (legend shown here only): Peak infections vs dif
# Panel B (legend centered inside, no box, thinner lines)
p2 <- ggplot(summary_table,
             aes(x = dif, y = peak_infections,
                 color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.5, show.legend = TRUE) +
  geom_point(size = 3, show.legend = FALSE) +
  labs(x = "", y = "Peak Infections") +
  scale_color_manual(values = scenario_colors, guide = "none") +
  scale_linetype_manual(values = line_types_common,
                        name = NULL,
                        labels = c("Baseline", "DFC", "HFC")) +
  guides(linetype = guide_legend(
    override.aes = list(color = "black", linewidth = 1.2),
    keywidth = unit(2.5, "cm")
  )) +
  big_theme +
  theme(
    legend.position = c(0.5, 0.2),
    legend.justification = c("center", "center"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.width = unit(2.5, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.text = element_text(size = 16)
  )

# Panel C (no legend): Harvest-window infection burden

harvest_table <- bind_rows(lapply(names(dif_results), function(dif_val) {
  res <- dif_results[[dif_val]]
  tibble(
    dif = as.numeric(dif_val),
    scenario = factor(
      c("Baseline", "Disease Control", "Harvest Control"),
      levels = c("Baseline", "Disease Control", "Harvest Control")
    ),
    harvest_IH = c(
      calculate_harvest_infections(res$baseline$IH_timeseries, years = YEARS,
                                   hs = harvest_start, hd = harvest_duration),
      calculate_harvest_infections(res$disease_control$IH_timeseries, years = YEARS,
                                   hs = harvest_start, hd = harvest_duration),
      calculate_harvest_infections(res$harvest_control$IH_timeseries, years = YEARS,
                                   hs = harvest_start, hd = harvest_duration)
    )
  )
}))

p3 <- ggplot(harvest_table,
             aes(x = dif, y = harvest_IH,
                 color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  labs(x = "", y = "Infection During Harvest") +
  scale_color_manual(values = scenario_colors) +
  scale_linetype_manual(values = line_types_common) +
  big_theme +
  theme(legend.position = "none")

# Panel D (no legend): Control efficacy ((baseline - control)/effort), only DFC and HFC
efficacy_df <- bind_rows(lapply(names(dif_results), function(d) {
  res <- dif_results[[d]]
  base_total_IH <- res$baseline$total_infections
  dfc_total_IH  <- res$disease_control$total_infections
  hfc_total_IH  <- res$harvest_control$total_infections
  dfc_effort    <- sum(res$disease_control$DFC_timeseries, na.rm = TRUE)
  hfc_effort    <- sum(res$harvest_control$HFC_timeseries, na.rm = TRUE)
  
  tibble(
    dif = as.numeric(d),
    scenario = c("Disease Control", "Harvest Control"),
    efficacy = c(
      ifelse(dfc_effort > 0, (base_total_IH - dfc_total_IH) / dfc_effort, NA_real_),
      ifelse(hfc_effort > 0, (base_total_IH - hfc_total_IH) / hfc_effort, NA_real_)
    )
  )
}))

p4 <- ggplot(efficacy_df,
             aes(x = dif, y = efficacy,
                 color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  labs(x = "Mosquito Seasonal Shift (days)",
       y = "Control Efficacy") +
  scale_color_manual(values = p4_colors, guide = "none") +
  scale_linetype_manual(values = p4_linetypes, guide = "none") +
  big_theme +
  theme(legend.position = "none")

# Combine, add A–D tags, and save
comb_shift <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")  # adds A, B, C, D

outfile <- sprintf("1year_disease_shift_r=%.3f_Omega=%g_%s.png",
                   R_FOR_PLOTS, OMEGA_FOR_PLOTS, format(Sys.Date(), "%Y-%m-%d"))
ggsave(outfile, comb_shift, width = 12, height = 13, dpi = 300)

cat("\nSaved plot:", outfile, "\n")

