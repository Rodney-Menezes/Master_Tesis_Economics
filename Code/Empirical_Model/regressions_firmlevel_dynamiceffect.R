#=====================================================
# Tesis: 'Efectos de la política monetaria sobre la dinámica de inversión en 
#         Economías Latinoamericanas mediante un modelo 
#         con empresas heterogéneas'
# Autor: Jose Rodney Menezes De la Cruz
# Sumbit: 2025
#=====================================================

# =======================================================
# Dynamic Local Projections - Modelos con Inversion Neta
# =======================================================

# 0) Cargar librerías necesarias
# --------------------------------------------------------
library(haven)
library(zoo)
library(dplyr)
library(purrr)
library(fixest)
library(readr)
library(broom)
library(ggplot2)
library(patchwork)
library(tidyr)
library(kableExtra)



# =========================
# Helpers (winsor/standard)
# =========================
winsorize <- function(x, p = 0.005) {
  lo <- stats::quantile(x, p, na.rm = TRUE)
  hi <- stats::quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

safe_scale <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  if (!is.finite(sigma) || sigma == 0) {
    return(ifelse(is.na(x), NA_real_, 0))
  }
  (x - mu) / sigma
}

# Winsoriza y estandariza leverage y dd a nivel país
prep_fin_vars <- function(df, p = 0.005) {
  df <- dplyr::ungroup(df)

  df %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      leverage_win = winsorize(leverage, p = p),
      dd_win       = winsorize(dd,       p = p),
      lev_std      = safe_scale(leverage_win),
      dd_std       = safe_scale(dd_win)
    ) %>%
    dplyr::ungroup()
}

# Winsoriza y estandariza el shock monetario a nivel país (z-score)
prep_shock_var <- function(df, p = 0.005) {
  df <- dplyr::ungroup(df)

  df %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      shock_win = winsorize(shock, p = p),
      shock_z   = safe_scale(shock_win)   # z-score por país
    ) %>%
    dplyr::ungroup()
}

# Winsoriza y estandariza controles firm-level a nivel país
prep_ctrl_var <- function(df, var_in, prefix, p = 0.005) {
  df <- dplyr::ungroup(df)

  win_sym <- rlang::sym(paste0(prefix, "_win"))
  std_sym <- rlang::sym(paste0(prefix, "_std"))
  var_sym <- rlang::sym(var_in)

  df %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      !!win_sym := winsorize(!!var_sym, p = p),
      !!std_sym := safe_scale(!!win_sym)
    ) %>%
    dplyr::ungroup()
}

# Crea rezagos (por defecto, un periodo) de los controles firm-level dentro de firma
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


# ===========================================
# Carga de datos y preprocesamiento consistente
# ===========================================
ruta_dta <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.dta"
df <- haven::read_dta(ruta_dta) %>%
  dplyr::mutate(dateq = zoo::as.yearqtr(dateq)) %>%
  dplyr::arrange(name, dateq)

# dlog_capital (Δ log K) y limpieza básica
df <- df %>%
  dplyr::group_by(name) %>%
  dplyr::arrange(dateq, .by_group = TRUE) %>%
  dplyr::mutate(
    ratio_cap    = real_capital / dplyr::lag(real_capital),
    dlog_capital = dplyr::if_else(is.finite(ratio_cap) & ratio_cap > 0,
                                  log(ratio_cap), NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-ratio_cap) %>%
  dplyr::filter(!is.na(dlog_capital))

# Shock: winsor + z-score global, luego cambio de signo
df <- df %>%
  prep_shock_var(p = 0.005) %>%
  dplyr::mutate(shock_std = -shock_z) %>%   # MP>0 ≈ recorte (expansivo)
  dplyr::select(-shock_win, -shock_z)    # limpiar intermedios; opcional: también -shock original

# Leverage y dd: winsor + estándar global (una sola vez)
df <- df %>%
  prep_fin_vars(p = 0.005)

# Niveles base por firma y país y dummies ex-ante (apalancamiento/distancia al default)
firm_fin_base <- df %>%
  dplyr::group_by(Country, name) %>%
  dplyr::summarise(
    lev_base = mean(leverage_win, na.rm = TRUE),
    dd_base  = mean(dd_win,       na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::group_by(Country) %>%
  dplyr::mutate(
    lev_p25 = stats::quantile(lev_base, 0.25, na.rm = TRUE),
    lev_p75 = stats::quantile(lev_base, 0.75, na.rm = TRUE),
    dd_p25  = stats::quantile(dd_base,  0.25, na.rm = TRUE),
    dd_p75  = stats::quantile(dd_base,  0.75, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    low_lev  = dplyr::if_else(lev_base <= lev_p25, 1L, 0L, missing = NA_integer_),
    high_lev = dplyr::if_else(lev_base >= lev_p75, 1L, 0L, missing = NA_integer_),
    near_def = dplyr::if_else(dd_base  <= dd_p25,  1L, 0L, missing = NA_integer_),
    far_def  = dplyr::if_else(dd_base  >= dd_p75,  1L, 0L, missing = NA_integer_)
  ) %>%
  dplyr::select(-lev_p25, -lev_p75, -dd_p25, -dd_p75)

df <- df %>%
  dplyr::left_join(firm_fin_base, by = c("Country", "name"))

# =================================================================
# Figure 1: Respuesta dinámica de la inversión al shock monetario
# =================================================================

# Parte 1: Sin controles cíclicos
# =================================================================

# 1) Controles firm-level y estandarización de ventas, tamaño y liquidez
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
    dplyr::arrange(dateq) %>%
    dplyr::mutate(
      size_raw = (log(atq) + log (saleq)) / 2
    ) %>%
    dplyr::ungroup()
}

df <- df %>%
  prep_ctrl_var(var_in = "current_ratio", prefix = "sh_current_a", p = 0.005)

df <- df %>%
  add_lagged_controls(c("rsales_g_std", "size_raw", "sh_current_a_std"))

# 2) Crear interacciones “raw” (winsorizadas) con el shock
df <- df %>%
mutate(
  lev_shock = lev_std * shock_std,
  d2d_shock = dd_std  * shock_std
)

# 3) Construir cumFh_dlog_capital para h = 0…12
# --------------------------------------------------------
vars_cap <- paste0("cumF", 0:12, "_dlog_capital")
if (!all(vars_cap %in% names(df))) {
  df_dyn_nocy <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    group_modify(~ {
      tmp <- .x
      for (h in 0:12) {
        tmp[[paste0("cumF", h, "_dlog_capital")]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  df_dyn_nocy <- df
}

# 4) Definir controles y cadena de regresores
controls_firm   <- c("L1_rsales_g_std", "L1_size_raw", "L1_sh_current_a_std")
controls_macro  <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn_nocy))
controls_vec_nocy <- c(controls_firm, controls_macro)
all_controls_nocy <- paste(controls_vec_nocy, collapse = " + ")

# 5) Estimar dinámicas h = 0…12 sin componente cíclico
res_lev_nocy <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ lev_shock + ", all_controls_nocy,
      " | name + sec + dateq"
    )),
    data    = df_dyn_nocy,
    cluster = ~ Country + dateq + name
  )
})
res_dd_nocy <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ d2d_shock + ", all_controls_nocy,
      " | name + sec + dateq"
    )),
    data    = df_dyn_nocy,
    cluster = ~ Country + dateq + name
  )
})

# 6) Extraer coeficientes y errores estándar
dyn_lev_nocy <- tibble(
  horizon = 0:12,
  beta    = map_dbl(res_lev_nocy, ~ coef(.x)["lev_shock"]),
  se      = map_dbl(res_lev_nocy, ~ sqrt(vcov(.x)["lev_shock","lev_shock"]))
)
dyn_dd_nocy <- tibble(
  horizon = 0:12,
  beta    = map_dbl(res_dd_nocy, ~ coef(.x)["d2d_shock"]),
  se      = map_dbl(res_dd_nocy, ~ sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
)

# 7) Graficar Figura 1
p_lev_nocy <- ggplot(dyn_lev_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              fill = "firebrick", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Panel (a): Heterogeneidad por apalancamiento",
    x        = "Trimestres",
    y        = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

p_dd_nocy <- ggplot(dyn_dd_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              fill = "steelblue", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Panel (b): Heterogeneidad por distancia al default",
    x        = "Trimestres",
    y        = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

figure1 <- (p_lev_nocy | p_dd_nocy) +
  plot_annotation(
    title = "Figura 1: Heterogeneidad financiera en la dinámica de la inversión ante un shock monetario expansivo"
  )

print(figure1)

# ============================================================
# Submuestras ex-ante: canales de apalancamiento y default
# ============================================================

estimate_lp_subset <- function(data_in) {
  purrr::map(0:12, function(h) {
    dep_var   <- paste0("cumF", h, "_dlog_capital")
    rhs_terms <- c("shock_std", controls_vec_nocy)
    rhs_terms <- rhs_terms[nzchar(rhs_terms)]
    rhs_str   <- paste(rhs_terms, collapse = " + ")
    fml_str   <- paste(dep_var, "~", rhs_str, "| name + sec + dateq")

    feols(
      stats::as.formula(fml_str),
      data    = data_in,
      cluster = ~ Country
    )
  })
}

extract_shock_path <- function(models) {
  tibble(
    horizon = 0:12,
    beta = purrr::map_dbl(models, function(m) {
      coefs <- coef(m)
      if ("shock_std" %in% names(coefs)) {
        coefs[["shock_std"]]
      } else {
        NA_real_
      }
    }),
    se = purrr::map_dbl(models, function(m) {
      vc <- vcov(m)
      if (!is.null(dim(vc)) &&
          "shock_std" %in% rownames(vc) &&
          "shock_std" %in% colnames(vc)) {
        sqrt(vc["shock_std", "shock_std"])
      } else {
        NA_real_
      }
    })
  )
}

df_low_lev  <- df_dyn_nocy %>% dplyr::filter(!is.na(low_lev),  low_lev == 1L)
df_high_lev <- df_dyn_nocy %>% dplyr::filter(!is.na(high_lev), high_lev == 1L)
df_near_def <- df_dyn_nocy %>% dplyr::filter(!is.na(near_def), near_def == 1L)
df_far_def  <- df_dyn_nocy %>% dplyr::filter(!is.na(far_def),  far_def == 1L)

res_low_lev  <- estimate_lp_subset(df_low_lev)
res_high_lev <- estimate_lp_subset(df_high_lev)
res_near_def <- estimate_lp_subset(df_near_def)
res_far_def  <- estimate_lp_subset(df_far_def)

dyn_low_lev <- extract_shock_path(res_low_lev) %>%
  dplyr::mutate(grupo = "P25: Bajo apalancamiento")
dyn_high_lev <- extract_shock_path(res_high_lev) %>%
  dplyr::mutate(grupo = "P75: Alto apalancamiento")
dyn_near_def <- extract_shock_path(res_near_def) %>%
  dplyr::mutate(grupo = "P25: Cercanas al default")
dyn_far_def <- extract_shock_path(res_far_def) %>%
  dplyr::mutate(grupo = "P75: Lejos del default")

lev_palette <- c(
  "P25: Bajo apalancamiento" = "#1f78b4",
  "P75: Alto apalancamiento" = "#e31a1c"
)
def_palette <- c(
  "P25: Cercanas al default" = "#33a02c",
  "P75: Lejos del default"   = "#6a3d9a"
)

dyn_lev_channel <- dplyr::bind_rows(dyn_low_lev, dyn_high_lev) %>%
  dplyr::mutate(grupo = factor(grupo, levels = names(lev_palette)))
dyn_def_channel <- dplyr::bind_rows(dyn_near_def, dyn_far_def) %>%
  dplyr::mutate(grupo = factor(grupo, levels = names(def_palette)))

p_leverage_channel <- ggplot(dyn_lev_channel,
                              aes(x = horizon, y = beta, color = grupo, fill = grupo)) +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              alpha = 0.15,
              color = NA) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = lev_palette) +
  scale_fill_manual(values = lev_palette) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Canal de apalancamiento",
    x     = "Trimestres",
    y     = "Efecto acumulado de inversión",
    color = "Submuestra",
    fill  = "Submuestra"
  ) +
  theme_minimal()

p_default_channel <- ggplot(dyn_def_channel,
                             aes(x = horizon, y = beta, color = grupo, fill = grupo)) +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              alpha = 0.15,
              color = NA) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = def_palette) +
  scale_fill_manual(values = def_palette) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Canal distancia al default",
    x     = "Trimestres",
    y     = "Efecto acumulado de inversión",
    color = "Submuestra",
    fill  = "Submuestra"
  ) +
  theme_minimal()

print(p_leverage_channel)
print(p_default_channel)


# =================================================================
# Figure 2: Respuesta Promedio de Inversión neta al Shock Monetario
# =================================================================


# 3) Construir variables dinámicas cumFh_dlog_capital para h = 0…12
# --------------------------------------------------------
vars_cap <- paste0("cumF", 0:12, "_dlog_capital")
df_dyn <- if (!all(vars_cap %in% names(df))) {
  df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    group_modify(~ {
      tmp <- .x
      for (h in 0:12) {
        tmp[[paste0("cumF", h, "_dlog_capital")]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  df
}

# 4) Definir controles firm-level y macro contemporáneos
# --------------------------------------------------------
controls_firm  <- c("L1_rsales_g_std", "L1_size_raw", "L1_sh_current_a_std")
controls_macro <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn))
all_controls   <- paste(c(controls_firm, controls_macro), collapse = " + ")

# 5) Estimar respuesta promedio al shock (horizonte 0–12)
# --------------------------------------------------------
res_avg <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var,
    " ~ shock_std + lev_shock + d2d_shock + ", all_controls,
    " | name + sec + dateq"
  ))
  feols(fml, data = df_dyn, cluster = ~ Country + dateq + name)
})

# 6) Extraer coeficientes y errores estándar para 'shock_std'
# -------------------------------------------------------
avg_coefs <- tibble::tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_avg, ~ coef(.x)["shock_std"]),
  se_shock   = map_dbl(res_avg, ~ sqrt(vcov(.x)["shock_std","shock_std"]))
)

# 7) Graficar la respuesta promedio al shock
# -------------------------------------------
p_avg <- ggplot(avg_coefs, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "darkgreen") +
  geom_point(size = 2, color = "darkgreen") +
  geom_ribbon(aes(
    ymin = beta_shock - 1.645 * se_shock,
    ymax = beta_shock + 1.645 * se_shock
  ), alpha = 0.2, fill = "darkgreen") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 2: Respuesta Dinámica de la Inversión ante un Shock Monetario",
    x     = "Trimestres",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

print(p_avg)




# ===================================================================================
# Figura 3: Dinámica heterognea de Inversión Neta Controlando por Inversión Rezagada
# ===================================================================================


# 1) Preprocesar y crear dlog_capital y Ldl_capital
# ---------------------------------------------------

if (!"Ldl_capital" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(Ldl_capital = lag(dlog_capital, 1)) %>%
    ungroup()
}

# 2) Winsorizar y estandarizar leverage y dd a nivel global
# ----------------------------------------------------------
df <- prep_fin_vars(df)

# 3) Interacciones con shock y con dlog_gdp (si existe)
# ------------------------------------------------------
if ("dlog_gdp" %in% names(df)) {
  df <- df %>%
    mutate(
      lev_shock_gdp = lev_std * dlog_gdp,
      d2d_shock_gdp = dd_std * dlog_gdp
    )
}

df <- df %>%
  mutate(
    lev_shock = lev_std * shock_std,
    d2d_shock = dd_std  * shock_std
  )

# 4) Construir cumFh_dlog_capital para h = 0…12
# --------------------------------------
vars_cap11 <- paste0("cumF", 0:12, "_dlog_capital")
df_dyn11 <- if (!all(vars_cap11 %in% names(df))) {
  df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    group_modify(~ {
      tmp <- .x
      for (h in 0:12) {
        tmp[[paste0("cumF", h, "_dlog_capital")]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  df
}

# 5) Definir controles firm-level y macro
# --------------------------------------
controls_firm  <- c("L1_rsales_g_std", "L1_size_raw", "L1_sh_current_a_std")
controls_macro <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn11))
base_controls  <- c(controls_firm, controls_macro)

# 6) Estimar dinámicas controlando Ldl_capital
# --------------------------------------
has_lev_gdp <- "lev_shock_gdp" %in% names(df_dyn11)
has_dd_gdp  <- "d2d_shock_gdp" %in% names(df_dyn11)

res_lev_lag <- map(0:12, function(h) {
  dep_var   <- paste0("cumF", h, "_dlog_capital")
  rhs_terms <- c(
    if (has_lev_gdp) "lev_shock_gdp",
    "lev_shock",
    "Ldl_capital",
    base_controls
  )
  fml_str   <- paste(dep_var, "~", paste(rhs_terms, collapse = " + "), "| name + sec + dateq")
  feols(as.formula(fml_str), data = df_dyn11, cluster = ~ Country + dateq + name)
})

res_dd_lag <- map(0:12, function(h) {
  dep_var   <- paste0("cumF", h, "_dlog_capital")
  rhs_terms <- c(
    if (has_dd_gdp) "d2d_shock_gdp",
    "d2d_shock",
    "Ldl_capital",
    base_controls
  )
  fml_str   <- paste(dep_var, "~", paste(rhs_terms, collapse = " + "), "| name + sec + dateq")
  feols(as.formula(fml_str), data = df_dyn11, cluster = ~ Country + dateq + name)
})

# 7) Extraer coeficientes y errores para plot
# --------------------------------------
coef_lev_lag <- tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_lev_lag, ~ coef(.x)["lev_shock"]),
  se_shock   = map_dbl(res_lev_lag, ~ sqrt(vcov(.x)["lev_shock","lev_shock"]))
)

coef_dd_lag <- tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_dd_lag, ~ coef(.x)["d2d_shock"]),
  se_shock   = map_dbl(res_dd_lag, ~ sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
)

# 8) Graficar respuestas dinámicas
# --------------------------------------
p11a <- ggplot(coef_lev_lag, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_ribbon(aes(ymin = beta_shock - 1.645 * se_shock,
                  ymax = beta_shock + 1.645 * se_shock),
              fill = "firebrick", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Panel (a): Heterogeneidad por apalancamiento",
    x     = "Trimestres",
    y     = "Efecto acumulado de la inversion residual"
  ) +
  theme_minimal()

p11b <- ggplot(coef_dd_lag, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta_shock - 1.645 * se_shock,
                  ymax = beta_shock + 1.645 * se_shock),
              fill = "steelblue", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Panel (b): Heterogeneidad por distancia al default",
    x     = "Trimestres",
    y     = "Efecto acumulado de la inversion residual"
  ) +
  theme_minimal()

# 9) Montaje final 2×1
# --------------------------------------
(p11a / p11b) +
  plot_annotation(
    title = "Figura 3: Heterogeneidad financiera en la dinámica de la inversión residual ante un shock monetario expansivo"
  )




#===============================================================
# Figure 4: Respuesta Promedio de Inversión al Shock Monetario
# ==============================================================

# 3) Construir Ldl_capital (lag de dlog_capital)
# --------------------------------------------------------
if (!"Ldl_capital" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(Ldl_capital = lag(dlog_capital, 1)) %>%
    ungroup()
}

# 4) Construir variables dinámicas cumFh_dlog_capital para h = 0…12
# --------------------------------------------------------
vars_cap <- paste0("cumF", 0:12, "_dlog_capital")
df_dyn <- if (!all(vars_cap %in% names(df))) {
  df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    group_modify(~ {
      tmp <- .x
      for (h in 0:12) {
        tmp[[paste0("cumF", h, "_dlog_capital")]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  df
}

# 5) Definir controles firm-level y macro contemporáneos
# --------------------------------------------------------
controls_firm  <- c("L1_rsales_g_std", "L1_size_raw", "L1_sh_current_a_std")
controls_macro <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn))
all_controls   <- paste(c(controls_firm, controls_macro), collapse = " + ")

# 6) Estimar respuesta promedio al shock (horizonte 0–12) con Ldl_capital
# --------------------------------------------------------
res_avg <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var,
    " ~ shock_std + lev_shock + d2d_shock + Ldl_capital + ", all_controls,
    " | name + sec + dateq"
  ))
  feols(fml, data = df_dyn, cluster = ~ Country + dateq + name)
})

# 7) Extraer coeficientes y errores estándar para 'shock_std'
# -------------------------------------------------------
avg_coefs <- tibble::tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_avg, ~ coef(.x)["shock_std"]),
  se_shock   = map_dbl(res_avg, ~ sqrt(vcov(.x)["shock_std","shock_std"]))
)

# 8) Graficar la respuesta promedio al shock
# -------------------------------------------
p_avg <- ggplot(avg_coefs, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "purple") +
  geom_point(size = 2, color = "purple") +
  geom_ribbon(aes(
    ymin = beta_shock - 1.645 * se_shock,
    ymax = beta_shock + 1.645 * se_shock
  ), alpha = 0.2, fill = "purple") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 4: Respuesta Dinámica de la Inversión residual ante un Shock Monetario",
    x     = "Trimestres",
    y     = "Efecto acumulado de inversión residual"
  ) +
  theme_minimal()

print(p_avg)


# ======================================================================
# Figura 5: Dinámica de Respuestas al Shock Monetario por Tamaño 
# ======================================================================

df_cntl <- df_cntl %>%
  dplyr::group_by(name) %>%
  dplyr::arrange(dateq, .by_group = TRUE) %>%
  dplyr::mutate(
    # crecimiento de ventas
    rsales_g = log(saleq) - log(dplyr::lag(saleq)),
    # tamaño bruto sin estandarizar
    size_raw = (log(atq) + log(saleq)) / 2
  ) %>%
  dplyr::ungroup() %>%
  prep_ctrl_var(var_in = "rsales_g", prefix = "rsales_g", p = 0.005) %>%
  prep_ctrl_var(var_in = "current_ratio", prefix = "sh_current_a", p = 0.005) %>%
  add_lagged_controls(c("rsales_g_std", "size_raw", "sh_current_a_std")) %>%
  dplyr::mutate(
    size_win = winsorize(size_raw)
  ) %>%
  dplyr::select(-size_raw) %>%
  prep_fin_vars()

# 2) Crear interacciones size × GDP y size × shock
# --------------------------------------------------------
if ("dlog_gdp" %in% names(df_cntl)) {
  df_cntl <- df_cntl %>%
    mutate(
      size_gdp   = if (!"size_gdp"   %in% names(.)) size_win * dlog_gdp else size_gdp
    )
} else {
  warning("No existe 'dlog_gdp'; se omite creación de 'size_gdp'.")
}

if ("shock_std" %in% names(df_cntl)) {
  df_cntl <- df_cntl %>%
    mutate(
      size_shock = if (!"size_shock" %in% names(.)) size_win * shock_std else size_shock
    )
} else {
  stop("Falta la variable 'shock_std' en el data frame.")
}

# 3) Construir dinámicas acumuladas cumFh_dlog_capital para h = 0…12
# ------------------------------------------------------------------
vars_cum <- paste0("cumF", 0:12, "_dlog_capital")
df_dyn12 <- if (!all(vars_cum %in% names(df_cntl))) {
  df_cntl %>%
    group_by(name) %>%
    arrange(dateq, .by_group = TRUE) %>%
    group_modify(~{
      tmp <- .x
      for (h in 0:12) {
        tmp[[ vars_cum[h+1] ]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  message("Variables dinámicas cumFh_dlog_capital ya existen; omitiendo.")
  df_cntl
}

# 4) Definir controles firm-level y macro presentes
# --------------------------------------------------
controls_firm  <- c("L1_rsales_g_std", "L1_sh_current_a_std")
controls_macro <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")
present_macro  <- intersect(controls_macro, names(df_dyn12))
if (length(present_macro) < length(controls_macro)) {
  warning("Faltan controles macro: ",
          paste(setdiff(controls_macro, present_macro), collapse = ", "),
          ". Se omitirán.")
}
all_controls <- paste(c(controls_firm, present_macro), collapse = " + ")


# 5) Estimar efecto heterogéneo size × shock para h = 0…12
# ---------------------------------------------------------
res_size <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      vars_cum[h+1], " ~ size_shock + ", all_controls,
      " | name + sec + dateq"
    )),
    data    = df_dyn12,
    cluster = ~ Country + dateq + name
  )
})

# 6) Extraer coeficiente y error estándar de size_shock
# ------------------------------------------------------
coef_size <- tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_size, ~ coef(.x)["size_shock"]),
  se_shock   = map_dbl(res_size, ~ sqrt(vcov(.x)["size_shock", "size_shock"]))
)

# 7) Graficar la dinámica de size_shock
# --------------------------------------
ggplot(coef_size, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "orange") +
  geom_point(size = 2, color = "orange") +
  geom_ribbon(aes(
    ymin = beta_shock - 1.645 * se_shock,
    ymax = beta_shock + 1.645 * se_shock
  ), alpha = 0.2, fill = "orange") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 5: Heterogeneidad de la inversion por tamaño de empresa",
    x     = "Trimestres",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()




# ==============================================================
# Figura 6: Dinámica Conjunta de Posición Financiera y Tamaño 
# ==============================================================

# Partimos de df original
df_cntl13 <- df

# 1) Calcular tamaño raw y winsorización
# ---------------------------------------
df_cntl13 <- df_cntl13 %>%
  dplyr::mutate(
    # tamaño proxy bruto sin estandarizar
    size_raw = (log(atq) + log (saleq)) / 2,
    size_win = winsorize(size_raw)
  ) %>%
  dplyr::select(-size_raw) %>%
  prep_fin_vars()

# 2) rear interacciones con shock y dlog_gdp
# -------------------------------------------
df_cntl13 <- df_cntl13 %>%
  mutate(
    size_shock    = size_win * shock_std,
    size_gdp      = size_win * dlog_gdp,
    lev_shock     = lev_std * shock_std,
    lev_shock_gdp = lev_std * dlog_gdp,
    d2d_shock     = dd_std  * shock_std,
    d2d_shock_gdp = dd_std  * dlog_gdp
  )

# 3) Construir dlog_capital y dinámicas cumFh_dlog_capital
# ---------------------------------------------------------
vars13 <- paste0("cumF", 0:12, "_dlog_capital")
df_dyn13 <- if (!all(vars13 %in% names(df_cntl13))) {
  df_cntl13 %>%
    group_by(name) %>%
    arrange(dateq) %>%
    group_modify(~ {
      tmp <- .x
      for (h in 0:12) {
        tmp[[vars13[h+1]]] <- rowSums(
          map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)),
          na.rm = TRUE
        )
      }
      tmp
    }) %>%
    ungroup()
} else df_cntl13

# 4) Definir controles firm-level y macro
# ----------------------------------------
controls_firm13  <- c("L1_rsales_g_std", "L1_sh_current_a_std")
controls_macro13 <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn13))
all_controls13   <- paste(c(controls_firm13, controls_macro13), collapse = " + ")

# 4a) Estimar dinámica Size + Leverage
res_lev_size <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      vars13[h+1],
      " ~ size_shock + lev_shock + lev_shock_gdp + ", all_controls13,
      " | name + sec + dateq"
    )),
    data    = df_dyn13,
    cluster = ~ Country + dateq + name
  )
})

# 4b) Estimar dinámica Size + DD
res_dd_size <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      vars13[h+1],
      " ~ size_shock + d2d_shock + d2d_shock_gdp + ", all_controls13,
      " | name + sec + dateq"
    )),
    data    = df_dyn13,
    cluster = ~ Country + dateq + name
  )
})

# 5) Extraer coeficientes y errores
# ----------------------------------
lev_size_coefs <- tibble(
  horizon   = 0:12,
  beta_size = map_dbl(res_lev_size, ~ coef(.x)["size_shock"]),
  se_size   = map_dbl(res_lev_size, ~ sqrt(vcov(.x)["size_shock","size_shock"])),
  beta_lev  = map_dbl(res_lev_size, ~ coef(.x)["lev_shock"]),
  se_lev    = map_dbl(res_lev_size, ~ sqrt(vcov(.x)["lev_shock","lev_shock"]))
)

dd_size_coefs <- tibble(
  horizon   = 0:12,
  beta_size = map_dbl(res_dd_size, ~ coef(.x)["size_shock"]),
  se_size   = map_dbl(res_dd_size, ~ sqrt(vcov(.x)["size_shock","size_shock"])),
  beta_dd   = map_dbl(res_dd_size, ~ coef(.x)["d2d_shock"]),
  se_dd     = map_dbl(res_dd_size, ~ sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
)

# 5a) Gráfico Size × Shock
p13a <- ggplot(lev_size_coefs, aes(x = horizon, y = beta_size)) +
  geom_line(size = 1, color = "darkgreen") +
  geom_ribbon(aes(ymin = beta_size - 1.645 * se_size,
                  ymax = beta_size + 1.645 * se_size),
              fill = "darkgreen", alpha = 0.2) +
    labs(title = "Panel (a): Heterogeneidad por tamaño", x = "Trimestres", y = "Efecto acumulado de inversión") +
  theme_minimal()

# 5b) Gráfico Leverage × Shock
p13b <- ggplot(lev_size_coefs, aes(x = horizon, y = beta_lev)) +
  geom_line(size = 1, color = "firebrick") +
  geom_ribbon(aes(ymin = beta_lev - 1.645 * se_lev,
                  ymax = beta_lev + 1.645 * se_lev),
              fill = "firebrick", alpha = 0.2) +
    labs(title = "Panel (b): Heterogeneidad por apalancamiento", x = "Trimestres", y = "Efecto acumulado de inversión") +
  theme_minimal()

# 5c) Gráfico DD × Shock
p13c <- ggplot(dd_size_coefs, aes(x = horizon, y = beta_dd)) +
  geom_line(size = 1, color = "steelblue") +
  geom_ribbon(aes(ymin = beta_dd - 1.645 * se_dd,
                  ymax = beta_dd + 1.645 * se_dd),
              fill = "steelblue", alpha = 0.2) +
    labs(title = "Panel (c): Heterogeneidad por distancia al default", x = "Trimestres", y = "Efecto acumulado de inversión") +
  theme_minimal()

# 6) Montaje final 3×1
# ----------------------
(p13a / p13b / p13c) +
  plot_annotation(
    title = "Figura 6: Heterogeneidad conjunta de la posición Financiera y tamaño en la dinamica de la inversion"
  )
