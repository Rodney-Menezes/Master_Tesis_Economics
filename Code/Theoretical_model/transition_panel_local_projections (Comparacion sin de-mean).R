# =============================================================
# Comparative Local Projections: Simulated vs Empirical Panels
# =============================================================
#
# Este script estima las respuestas dinámicas de la inversión ante un shock
# monetario utilizando (i) los paneles simulados del modelo teórico y (ii) los
# datos firm-level empleados en el análisis empírico. Posteriormente, grafica
# ambas dinámicas en una única figura para facilitar la comparación.
#
# Para modificar las rutas de acceso se pueden definir las variables de
# entorno:
#   - TRANSITION_PANEL_DIR: carpeta con los CSV "mTransitionPanel_*.csv".
#   - TRANSITION_OUTPUT_DIR: carpeta donde se guardarán los resultados
#     intermedios del modelo teórico (opcional, sólo para depuración).
#   - EMPIRICAL_DATA_PATH: ruta al archivo .dta utilizado en el análisis
#     empírico.
#   - COMPARISON_OUTPUT_DIR: carpeta donde se guardará la figura comparativa.
#
# =============================================================

suppressPackageStartupMessages({
  library(readr)
  library(haven)
  library(zoo)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(fixest)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------
# Configuración
# ------------------------------
# Horizonte máximo de interés
sim_horizons <- 0:13
emp_horizons <- 0:12

# Parámetros del shock monetario en el modelo teórico
shock_length <- 12
shock_size   <- -0.0025
shock_decay  <- 0.9
shock_scale  <- -4

# Directorios y rutas de archivos
default_sim_dir   <- "C:/Users/joser/Documents"
transition_dir    <- Sys.getenv("TRANSITION_PANEL_DIR", unset = default_sim_dir)
transition_outdir <- Sys.getenv("TRANSITION_OUTPUT_DIR", unset = getwd())

default_emp_path <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.dta"
empirical_path   <- Sys.getenv("EMPIRICAL_DATA_PATH", unset = default_emp_path)
comparison_outdir <- Sys.getenv("COMPARISON_OUTPUT_DIR", unset = getwd())

if (!dir.exists(transition_dir)) {
  stop(
    "No se encontró la carpeta con los paneles de transición: ",
    transition_dir,
    "\nUtiliza Sys.setenv(TRANSITION_PANEL_DIR = 'ruta/a/paneles') antes de ejecutar."
  )
}

if (!dir.exists(comparison_outdir)) {
  dir.create(comparison_outdir, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------
# Funciones auxiliares (panel simulado)
# ------------------------------
parse_tpre <- function(file_name) {
  value <- stringr::str_extract(basename(file_name), "(?<=_)\\d+(?=\\.csv$)")
  ifelse(is.na(value), NA_integer_, as.integer(value))
}

build_shock_series <- function(max_period, t_pre) {
  shock <- rep(0, max_period)
  if (is.na(t_pre) || t_pre + 1 > max_period) {
    return(shock)
  }
  for (h in seq_len(shock_length)) {
    pos <- t_pre + h
    if (pos > max_period) {
      break
    }
    shock[pos] <- shock_size * (shock_decay^(h - 1))
  }
  shock * shock_scale
}

prepare_panel <- function(file_path) {
  raw <- readr::read_csv(
    file_path,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  names(raw) <- c(
    "firm_id", "quarter_id", "in_sample", "balanced_panel",
    "investment", "cash_flow", "market_value", "capital", "cash",
    "debt", "employment", "interest_rate", "unconstrained"
  )
  
  df <- raw %>%
    dplyr::filter(in_sample == 1) %>%
    dplyr::arrange(firm_id, quarter_id)
  
  df <- df %>%
    dplyr::mutate(
      log_capital = ifelse(capital > 0, log(capital), NA_real_),
      leverage = dplyr::if_else(capital != 0, debt / capital, NA_real_),
      cash_to_capital = dplyr::if_else(capital != 0, cash / capital, NA_real_)
    )
  
  df %>%
    dplyr::group_by(firm_id) %>%
    dplyr::mutate(
      lag_log_capital = dplyr::lag(log_capital),
      lag_leverage = dplyr::lag(leverage),
      lag_dd = dplyr::lag(cash_to_capital)
    ) %>%
    dplyr::ungroup()
}

run_local_projection <- function(df, interaction_col, horizon) {
  response_col <- paste0("response_h", horizon)
  
  df_h <- df %>%
    dplyr::group_by(firm_id) %>%
    dplyr::mutate(
      !!response_col := dplyr::lead(log_capital, horizon) - lag_log_capital
    ) %>%
    dplyr::ungroup()
  
  reg_data <- tibble::tibble(
    firm_id    = df_h$firm_id,
    quarter_id = df_h$quarter_id,
    response   = df_h[[response_col]],
    interaction = df_h[[interaction_col]]
  ) %>%
    dplyr::filter(!is.na(response), !is.na(interaction))
  
  if (nrow(reg_data) == 0) {
    return(tibble(
      horizon = horizon,
      coefficient = NA_real_,
      std_error = NA_real_,
      n_obs = 0L
    ))
  }
  
  has_within_quarter_variation <- reg_data %>%
    dplyr::group_by(quarter_id) %>%
    dplyr::summarise(n_unique = dplyr::n_distinct(interaction), .groups = "drop") %>%
    dplyr::filter(n_unique > 1) %>%
    nrow() > 0
  
  if (!has_within_quarter_variation) {
    return(tibble(
      horizon = horizon,
      coefficient = NA_real_,
      std_error = NA_real_,
      n_obs = 0L
    ))
  }
  
  model <- tryCatch(
    fixest::feols(
      response ~ interaction,
      data = reg_data,
      fixef = c("firm_id", "quarter_id"),
      cluster = ~firm_id
    ),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    return(tibble(
      horizon = horizon,
      coefficient = NA_real_,
      std_error = NA_real_,
      n_obs = as.integer(nrow(reg_data))
    ))
  }
  
  vcov_mat <- stats::vcov(model)
  coef_val <- stats::coef(model)["interaction"]
  se_val   <- sqrt(diag(vcov_mat))["interaction"]
  
  tibble(
    horizon = horizon,
    coefficient = unname(coef_val),
    std_error = unname(se_val),
    n_obs = stats::nobs(model)
  )
}

weighted_mean_safe <- function(x, w) {
  valid <- !is.na(x) & !is.na(w) & w > 0
  if (!any(valid)) {
    return(NA_real_)
  }
  sum(x[valid] * w[valid]) / sum(w[valid])
}

compute_simulated_dynamics <- function(panel_dir, horizons) {
  panel_files <- list.files(
    path = panel_dir,
    pattern = "^mTransitionPanel_\\d+\\.csv$",
    full.names = TRUE
  )
  
  horizon_grid <- sort(unique(c(0, horizons)))
  
  if (length(panel_files) == 0) {
    stop("No se encontraron archivos 'mTransitionPanel_*.csv' en ", panel_dir, ".")
  }
  
  message("Procesando ", length(panel_files), " panel(es) simulados en: ", panel_dir)
  
  results <- purrr::map_dfr(panel_files, function(path) {
    t_pre <- parse_tpre(path)
    df <- prepare_panel(path)
    
    max_q <- max(df$quarter_id, na.rm = TRUE)
    shock_series <- build_shock_series(max_q, t_pre)
    
    df <- df %>%
      dplyr::mutate(
        shock = shock_series[quarter_id],
        lev_shock = lag_leverage * shock,
        dd_shock = lag_dd * shock
      )
    
    message("  • ", basename(path), " (t_pre = ", t_pre, ")")
    
    lev_results <- purrr::map_dfr(horizon_grid, ~ run_local_projection(df, "lev_shock", .x)) %>%
      dplyr::mutate(measure = "leverage")
    
    dd_results <- purrr::map_dfr(horizon_grid, ~ run_local_projection(df, "dd_shock", .x)) %>%
      dplyr::mutate(measure = "default_distance")
    
    dplyr::bind_rows(lev_results, dd_results) %>%
      dplyr::mutate(
        panel_file = basename(path),
        t_pre = t_pre
      )
  })
  
  summary_results <- results %>%
    dplyr::group_by(measure, horizon) %>%
    dplyr::summarise(
      coefficient = weighted_mean_safe(coefficient, n_obs),
      std_error = sqrt(weighted_mean_safe(std_error^2, n_obs)),
      n_obs = sum(n_obs, na.rm = TRUE),
      n_panels = dplyr::n_distinct(panel_file),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      measure,
      horizon = horizon_grid,
      fill = list(
        coefficient = NA_real_,
        std_error = NA_real_,
        n_obs = 0,
        n_panels = 0
      )
    ) %>%
    dplyr::arrange(measure, horizon)
  
  list(raw = results, summary = summary_results)
}

# ------------------------------
# Funciones auxiliares (panel empírico)
# ------------------------------
winsorize <- function(x, p = 0.005) {
  lo <- stats::quantile(x, p, na.rm = TRUE)
  hi <- stats::quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

prep_fin_vars <- function(df, p = 0.005) {
  df %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      leverage_win = winsorize(leverage, p = p),
      dd_win       = winsorize(dd,       p = p)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(dateq, .by_group = TRUE) %>%
    dplyr::mutate(
      lev_dm   = leverage_win - mean(leverage_win, na.rm = TRUE),
      dd_dm    = dd_win       - mean(dd_win,       na.rm = TRUE),
      L1_lev_dm = dplyr::lag(lev_dm, 1),
      L1_dd_dm  = dplyr::lag(dd_dm,  1)
    ) %>%
    dplyr::ungroup()
}

prep_shock_var <- function(df, p = 0.005) {
  df %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      shock_win = winsorize(shock, p = p),
      shock_exp = -shock_win
    ) %>%
    dplyr::ungroup()
}

prep_ctrl_var <- function(df, var_in, prefix = var_in, p = 0.005) {
  df <- dplyr::ungroup(df)
  win_sym <- rlang::sym(paste0(prefix, "_win"))
  var_sym <- rlang::sym(var_in)
  
  df %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      !!win_sym := winsorize(!!var_sym, p = p)
    ) %>%
    dplyr::ungroup()
}

add_lagged_controls <- function(df, vars, lag = 1) {
  vars_present <- intersect(vars, names(df))
  if (length(vars_present) == 0) {
    return(df)
  }
  
  lag_prefix <- paste0("L", lag, "_")
  
  df %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(dateq, .by_group = TRUE) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(vars_present),
      ~ dplyr::lag(.x, lag),
      .names = paste0(lag_prefix, "{.col}")
    )) %>%
    dplyr::ungroup()
}

compute_empirical_dynamics <- function(data_path, horizons) {
  if (!file.exists(data_path)) {
    stop(
      "No se encontró el archivo de datos empíricos: ",
      data_path,
      "\nUtiliza Sys.setenv(EMPIRICAL_DATA_PATH = 'ruta/al/archivo.dta') antes de ejecutar."
    )
  }
  
  horizon_grid <- sort(unique(c(0, horizons)))
  
  df <- haven::read_dta(data_path) %>%
    dplyr::mutate(dateq = zoo::as.yearqtr(dateq)) %>%
    dplyr::arrange(name, dateq)
  
  df <- df %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(dateq, .by_group = TRUE) %>%
    dplyr::mutate(
      ratio_cap    = real_capital / dplyr::lag(real_capital),
      dlog_capital = dplyr::if_else(
        is.finite(ratio_cap) & ratio_cap > 0,
        log(ratio_cap),
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-ratio_cap) %>%
    dplyr::filter(!is.na(dlog_capital))
  
  df <- df %>%
    prep_shock_var(p = 0.005) %>%
    prep_fin_vars(p = 0.005)
  
  df <- df %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(dateq, .by_group = TRUE) %>%
    dplyr::mutate(
      rsales_g = log(saleq) - log(dplyr::lag(saleq))
    ) %>%
    dplyr::ungroup()
  
  df <- df %>%
    prep_ctrl_var(var_in = "rsales_g", prefix = "rsales_g", p = 0.005)
  
  if (!"size_raw" %in% names(df)) {
    df <- df %>%
      dplyr::group_by(name) %>%
      dplyr::arrange(dateq, .by_group = TRUE) %>%
      dplyr::mutate(size_raw = (log(atq) + log(saleq)) / 2) %>%
      dplyr::ungroup()
  }
  
  df <- df %>%
    prep_ctrl_var(var_in = "current_ratio", prefix = "current_ratio", p = 0.005) %>%
    add_lagged_controls(c("rsales_g_win", "size_raw", "current_ratio_win")) %>%
    dplyr::mutate(
      lev_shock = L1_lev_dm * shock_exp,
      d2d_shock = L1_dd_dm  * shock_exp
    )
  
  vars_cap <- paste0("cumF", horizon_grid, "_dlog_capital")
  df_dyn <- if (!all(vars_cap %in% names(df))) {
    df %>%
      dplyr::group_by(name) %>%
      dplyr::arrange(dateq, .by_group = TRUE) %>%
      dplyr::group_modify(~ {
        tmp <- .x
        for (h in horizon_grid) {
          tmp[[paste0("cumF", h, "_dlog_capital")]] <-
            rowSums(purrr::map_dfc(0:h, ~ dplyr::lead(tmp$dlog_capital, .x)), na.rm = TRUE)
        }
        tmp
      }) %>%
      dplyr::ungroup()
  } else {
    df
  }
  
  controls_firm  <- c("L1_rsales_g_win", "L1_size_raw", "L1_current_ratio_win")
  controls_macro <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn))
  controls_vec   <- c(controls_firm, controls_macro)
  all_controls   <- paste(controls_vec, collapse = " + ")
  
  lev_models <- purrr::map(horizon_grid, function(h) {
    fixest::feols(
      as.formula(paste0(
        "cumF", h, "_dlog_capital ~ lev_shock + ", all_controls,
        " | name + sec + dateq"
      )),
      data = df_dyn,
      cluster = ~ Country + dateq + name
    )
  })
  
  dd_models <- purrr::map(horizon_grid, function(h) {
    fixest::feols(
      as.formula(paste0(
        "cumF", h, "_dlog_capital ~ d2d_shock + ", all_controls,
        " | name + sec + dateq"
      )),
      data = df_dyn,
      cluster = ~ Country + dateq + name
    )
  })
  
  dyn_lev <- tibble(
    measure = "leverage",
    horizon = horizon_grid,
    coefficient = purrr::map_dbl(lev_models, ~ stats::coef(.x)["lev_shock"]),
    std_error  = purrr::map_dbl(lev_models, ~ sqrt(stats::vcov(.x)["lev_shock", "lev_shock"]))
  )
  
  dyn_dd <- tibble(
    measure = "default_distance",
    horizon = horizon_grid,
    coefficient = purrr::map_dbl(dd_models, ~ stats::coef(.x)["d2d_shock"]),
    std_error  = purrr::map_dbl(dd_models, ~ sqrt(stats::vcov(.x)["d2d_shock", "d2d_shock"]))
  )
  
  bind_rows(dyn_lev, dyn_dd) %>%
    tidyr::complete(
      measure,
      horizon = horizon_grid,
      fill = list(
        coefficient = NA_real_,
        std_error = NA_real_
      )
    ) %>%
    dplyr::arrange(measure, horizon)
}

# ------------------------------
# Ejecución principal
# ------------------------------

sim_results <- compute_simulated_dynamics(transition_dir, sim_horizons)
emp_results <- compute_empirical_dynamics(empirical_path, emp_horizons)

common_horizons <- sort(intersect(sim_results$summary$horizon, emp_results$horizon))
if (length(common_horizons) == 0) {
  stop("No hay horizontes comunes entre los resultados simulados y empíricos.")
}

sim_summary <- sim_results$summary %>%
  dplyr::filter(horizon %in% common_horizons) %>%
  dplyr::mutate(dataset = "Simulado (Modelo teórico)")

emp_summary <- emp_results %>%
  dplyr::filter(horizon %in% common_horizons) %>%
  dplyr::mutate(dataset = "Empírico (Datos firm-level)")

plot_data <- dplyr::bind_rows(sim_summary, emp_summary) %>%
  dplyr::mutate(
    ribbon_low  = coefficient - 1.645 * std_error,
    ribbon_high = coefficient + 1.645 * std_error
  )

if (nrow(plot_data) == 0) {
  stop("No hay observaciones válidas para graficar.")
}

palette_colors <- c(
  "Simulado (Modelo teórico)" = "firebrick",
  "Empírico (Datos firm-level)" = "steelblue"
)

plot_measure <- function(data, panel_title, horizons_axis) {
  ggplot(data, aes(x = horizon, y = coefficient, color = dataset, linetype = dataset)) +
    geom_ribbon(
      aes(ymin = ribbon_low, ymax = ribbon_high, fill = dataset),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = palette_colors) +
    scale_fill_manual(values = palette_colors) +
    scale_x_continuous(
      breaks = horizons_axis,
      limits = range(horizons_axis)
    ) +
    labs(
      title = panel_title,
      x = "Trimestres",
      y = "Efecto acumulado en la inversión",
      color = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "plain")
    )
}

plot_lev <- plot_measure(
  plot_data %>% dplyr::filter(measure == "leverage"),
  "Panel (a): Heterogeneidad por apalancamiento",
  common_horizons
)

plot_dd <- plot_measure(
  plot_data %>% dplyr::filter(measure == "default_distance"),
  "Panel (b): Heterogeneidad por distancia al default",
  common_horizons
)

comparison_plot <- (plot_lev | plot_dd) +
  patchwork::plot_annotation(
    title = paste0("Figura comparativa: Respuesta dinámica de la inversión ante un shock monetario\n",
                   "Simulaciones del modelo vs evidencia empírica"
    )
  )

print(comparison_plot)
print(plot_lev)
print(plot_dd)

output_path <- file.path(comparison_outdir, "comparacion_simulado_empirico_local_projections.png")

ggplot2::ggsave(
  filename = output_path,
  plot = comparison_plot,
  width = 11,
  height = 5,
  dpi = 300
)

message("Figura comparativa guardada en: ", output_path)

ggplot2::ggsave(
  filename = file.path(comparison_outdir, "comparacion_simulado_empirico_local_projections_panel_a.png"),
  plot = plot_lev,
  width = 5.5,
  height = 5,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(comparison_outdir, "comparacion_simulado_empirico_local_projections_panel_b.png"),
  plot = plot_dd,
  width = 5.5,
  height = 5,
  dpi = 300
)

# Guardar resultados agregados para referencia
summary_path <- file.path(comparison_outdir, "comparacion_simulado_empirico_local_projections.csv")
readr::write_csv(plot_data, summary_path)
message("Resultados agregados guardados en: ", summary_path)