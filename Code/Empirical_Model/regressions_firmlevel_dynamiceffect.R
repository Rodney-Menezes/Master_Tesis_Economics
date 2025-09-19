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

# Winsoriza y calcula (i) de-mean y (ii) estándar dentro de firma
prep_fin_vars <- function(df, p = 0.005) {
  df %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(
      leverage_win = winsorize(leverage, p = p),
      dd_win       = winsorize(dd,       p = p),
      # De-mean dentro de firma
      lev_std       = leverage_win - mean(leverage_win, na.rm = TRUE),
      dd_dm        = dd_win       - mean(dd_win,       na.rm = TRUE),
      # Estandarización dentro de firma (comparabilidad entre gráficos)
      lev_std      = (leverage_win - mean(leverage_win, na.rm = TRUE)) /
        sd(leverage_win,  na.rm = TRUE),
      dd_std       = (dd_win       - mean(dd_win,       na.rm = TRUE)) /
        sd(dd_win,        na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
}

# Winsoriza y estandariza el shock monetario por país (z-score)
prep_shock_var <- function(df, p = 0.005) {
  safe_zscore <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    if (!is.finite(sigma) || sigma == 0) return(ifelse(is.na(x), NA_real_, 0))
    (x - mu) / sigma
  }
  df %>%
    dplyr::group_by(Country) %>%
    dplyr::mutate(
      shock_win = winsorize(shock, p = p),
      shock_z   = safe_zscore(shock_win)   # z-score por país
    ) %>%
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

# Shock: winsor + z-score por país, luego cambio de signo y renombrar a 'shock1'
df <- df %>%
  prep_shock_var(p = 0.005) %>%
  dplyr::mutate(shock1 = -shock_z) %>%   # MP>0 ≈ recorte (expansivo)
  dplyr::select(-shock_win, -shock_z)    # limpiar intermedios; opcional: también -shock original

# Leverage y dd: winsor + de-mean + estándar dentro de firma (una sola vez)
df <- df %>%
  prep_fin_vars(p = 0.005)


# =================================================================
# Figure 1: Respuesta dinámica de la inversión al shock monetario
# =================================================================

# Parte 1: Sin controles cíclicos
# =================================================================

# 1) Controles firm-level y estandarización de ventas, tamaño y liquidez
if (!"rsales_g_std" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      rsales_g_std = {
        x <- log(saleq) - log(lag(saleq))
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      }
    ) %>%
    ungroup()
}

if (!"size_index" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      size_index = {
        x <- log(atq)
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      }
    ) %>%
    ungroup()
}

if (!"sh_current_a_std" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      sh_current_a_std = {
        x <- current_ratio
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      }
    ) %>%
    ungroup()
}

# 2) Crear interacciones “raw” (winsorizadas) con el shock
df <- df %>%
mutate(
  lev_shock = lev_std * shock,
  d2d_shock = dd_dm  * shock
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
controls_firm   <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro  <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn_nocy))
all_controls_nocy <- paste(c(controls_firm, controls_macro), collapse = " + ")

# 5) Estimar dinámicas h = 0…12 sin componente cíclico
res_lev_nocy <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ lev_shock + ", all_controls_nocy,
      " | name + sec + dateq"
    )),
    data    = df_dyn_nocy,
    cluster = ~ Country
  )
})
res_dd_nocy <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ d2d_shock + ", all_controls_nocy,
      " | name + sec + dateq"
    )),
    data    = df_dyn_nocy,
    cluster = ~ Country
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

# 7) Graficar sin control cíclico
p_lev_nocy <- ggplot(dyn_lev_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              fill = "firebrick", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
     title    = "Panel (a): Heterogeneidad por apalancamiento",
    subtitle = "Sin control cíclico",
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
    title    = "Panel (a): Heterogeneidad por distancia al default",
    subtitle = "Sin control cíclico",
    x        = "Trimestres",
    y        = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

# Parte 2: Con controles cíclicos
# ==================================

# 1) Crear sensibilidad cíclica por firma
sens_cyc <- df %>%
  group_by(name) %>%
  filter(is.finite(dlog_capital), is.finite(dlog_gdp), n() >= 6) %>%
  do(tidy(lm(dlog_capital ~ dlog_gdp, data = .))) %>%
  ungroup() %>%
  filter(term == "dlog_gdp") %>%
  select(name, sens_ciclica = estimate)

# 2) Unir y crear choque cíclico
df_dyn_cyc <- df_dyn_nocy %>%
  left_join(sens_cyc, by = "name") %>%
  mutate(
    sens_ciclica   = replace_na(sens_ciclica, 0),
    sens_cyc_shock = sens_ciclica * shock
  )

# 3) Cadena de controles cíclicos
all_controls_cyc <- paste(c(controls_firm, controls_macro, "sens_cyc_shock"), collapse = " + ")

# 4) Estimar dinámicas con componente cíclico
res_lev_cyc <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ lev_shock + ", all_controls_cyc,
      " | name + sec + dateq"
    )),
    data    = df_dyn_cyc,
    cluster = ~ Country
  )
})
res_dd_cyc <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ d2d_shock + ", all_controls_cyc,
      " | name + sec + dateq"
    )),
    data    = df_dyn_cyc,
    cluster = ~ Country
  )
})

# 5) Extraer coeficientes y errores estándar con control cíclico
dyn_lev_cyc <- tibble(
  horizon = 0:12,
  beta    = map_dbl(res_lev_cyc, ~ coef(.x)["lev_shock"]),
  se      = map_dbl(res_lev_cyc, ~ sqrt(vcov(.x)["lev_shock","lev_shock"]))
)
dyn_dd_cyc <- tibble(
  horizon = 0:12,
  beta    = map_dbl(res_dd_cyc, ~ coef(.x)["d2d_shock"]),
  se      = map_dbl(res_dd_cyc, ~ sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
)

# 6) Graficar con control cíclico
p_lev_cyc <- ggplot(dyn_lev_cyc, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              fill = "firebrick", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Panel (b): Heterogeneidad por apalancamiento",
    subtitle = "Con control cíclico",
    x        = "Trimestres",
    y        = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

p_dd_cyc <- ggplot(dyn_dd_cyc, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta - 1.645 * se,
                  ymax = beta + 1.645 * se),
              fill = "steelblue", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Panel (b): Heterogeneidad por distancia al default",
    subtitle = "Con control cíclico",
    x        = "Trimestres",
    y        = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

# 7) Montaje final en dos figuras (Figura 1a y Figura 1b)
figure1a <- (p_lev_nocy / p_lev_cyc) +
  plot_annotation(
    title = "Figura 1a: Heterogeneidad por apalancamiento en la dinámica de la inversión ante un shock monetario expansivo"
  )

figure1b <- (p_dd_nocy / p_dd_cyc) +
  plot_annotation(
    title = "Figura 1b: Heterogeneidad por distancia al default en la dinámica de la inversión ante un shock monetario expansivo"
  )

figure1c <- (p_lev_nocy | p_dd_nocy) +
  plot_annotation(
    title = "Figura 1c: Heterogeneidad financiera en la dinámica de la inversión ante un shock monetario expansivo"
  )

print(figure1a)
print(figure1b)
print(figure1c)


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
controls_firm  <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn))
all_controls   <- paste(c(controls_firm, controls_macro), collapse = " + ")

# 5) Estimar respuesta promedio al shock (horizonte 0–12)
# --------------------------------------------------------
res_avg <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var,
    " ~ shock + lev_shock + d2d_shock + ", all_controls,
    " | name + sec + dateq"
  ))
  feols(fml, data = df_dyn, cluster = ~ Country)
})

# 6) Extraer coeficientes y errores estándar para 'shock'
# -------------------------------------------------------
avg_coefs <- tibble::tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_avg, ~ coef(.x)["shock"]),
  se_shock   = map_dbl(res_avg, ~ sqrt(vcov(.x)["shock","shock"]))
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
# --------------------------------------

if (!"Ldl_capital" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(Ldl_capital = lag(dlog_capital, 1)) %>%
    ungroup()
}

# 2) Winsorizar y diazotizar leverage y dd dentro de firma
# ------------------------------------------------------
df <- prep_fin_vars(df)

# 3) Interacciones con shock y con dlog_gdp (si existe)
# ------------------------------------------------------
if ("dlog_gdp" %in% names(df)) {
  df <- df %>%
    mutate(
      lev_shock_gdp = lev_std * dlog_gdp,
      d2d_shock_gdp = dd_dm * dlog_gdp
    )
}

df <- df %>%
  mutate(
    lev_shock = lev_std * shock,
    d2d_shock = dd_dm  * shock
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
controls_firm  <- c("rsales_g_std", "size_index", "sh_current_a_std")
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
  feols(as.formula(fml_str), data = df_dyn11, cluster = ~ Country)
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
  feols(as.formula(fml_str), data = df_dyn11, cluster = ~ Country)
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
controls_firm  <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro <- intersect(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), names(df_dyn))
all_controls   <- paste(c(controls_firm, controls_macro), collapse = " + ")

# 6) Estimar respuesta promedio al shock (horizonte 0–12) con Ldl_capital
# --------------------------------------------------------
res_avg <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var,
    " ~ shock + lev_shock + d2d_shock + Ldl_capital + ", all_controls,
    " | name + sec + dateq"
  ))
  feols(fml, data = df_dyn, cluster = ~ Country)
})

# 7) Extraer coeficientes y errores estándar para 'shock'
# -------------------------------------------------------
avg_coefs <- tibble::tibble(
  horizon    = 0:12,
  beta_shock = map_dbl(res_avg, ~ coef(.x)["shock"]),
  se_shock   = map_dbl(res_avg, ~ sqrt(vcov(.x)["shock","shock"]))
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

# Partimos de df original
df_cntl <- df

# 1) Calcular ventas, liquidez y winsorización de tamaño y apalancamiento
# ------------------------------------------------------------------------

df_cntl <- df_cntl %>%
group_by(name) %>%
arrange(dateq, .by_group = TRUE) %>%
mutate(
    #  crecimiento de ventas
    rsales_g = if (!"rsales_g" %in% names(.)) log(saleq) - log(lag(saleq)) else rsales_g,
    rsales_g_std = if (!"rsales_g_std" %in% names(.))
      (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE)
    else rsales_g_std,
    
    # winsorización de tamaño bruto
    size_raw = (log(atq) + saleq) / 2,
    size_win = winsorize(size_raw),
    
    # liquidez corriente estandarizada
    sh_current_a_std = if (!"sh_current_a_std" %in% names(.))
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE)
    else sh_current_a_std
) %>%
  ungroup() %>%
  select(-size_raw) %>%
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

if ("shock" %in% names(df_cntl)) {
  df_cntl <- df_cntl %>%
    mutate(
      size_shock = if (!"size_shock" %in% names(.)) size_win * shock else size_shock
    )
} else {
  stop("Falta la variable 'shock' en el data frame.")
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
controls_firm  <- c("rsales_g_std", "sh_current_a_std")
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
    cluster = ~ Country
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
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # tamaño proxy bruto sin estandarizar
    size_raw = (log(atq) + saleq) / 2,
    size_win = winsorize(size_raw)
  ) %>%
  ungroup() %>%
  select(-size_raw) %>%
  prep_fin_vars()

# 2) rear interacciones con shock y dlog_gdp
# -------------------------------------------
df_cntl13 <- df_cntl13 %>%
  mutate(
    size_shock    = size_win * shock,
    size_gdp      = size_win * dlog_gdp,
    lev_shock     = lev_std * shock,
    lev_shock_gdp = lev_std * dlog_gdp,
    d2d_shock     = dd_dm  * shock,
    d2d_shock_gdp = dd_dm  * dlog_gdp
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
controls_firm13  <- c("rsales_g_std", "sh_current_a_std")
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
    cluster = ~ Country
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
    cluster = ~ Country
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
