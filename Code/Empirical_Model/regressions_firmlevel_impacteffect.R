#=====================================================
# Tesis: 'Efectos de la política monetaria sobre la dinámica de inversión en 
#         Economías Latinoamericanas mediante un modelo 
#         con empresas heterogéneas'
# Autor: Jose Rodney Menezes De la Cruz
# Sumbit: 2025
#=====================================================

# =====================================
# Modelos de Panel de Efectos Fijos
# =====================================

# ================================================
# 0) Instalar paquetes necesarios (si no están instalados)
# ================================================
required_packages <- c(
  "modelsummary", "fixest", "haven", "dplyr",
  "tidyr", "zoo", "tibble", "writexl"
)
new_packages <- setdiff(required_packages, installed.packages()[, "Package"])
if (length(new_packages) > 0) {
  install.packages(new_packages)
}

# ================================================
# 1) Cargar librerías
# ================================================
library(haven)        # read_dta()
library(zoo)          # as.yearqtr()
library(fixest)       # feols()
library(modelsummary) # presentación de tablas
library(tibble)       # tibble
library(writexl)      # para exportar
library(tidyr)        # drop_na()
library(dplyr)        # manipulación de datos

# ================================================
# 2) Definir funciones auxiliares
# ================================================
winsorize <- function(x, p = 0.005) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

run_model <- function(df, vars) {
  tmp <- df %>%
    filter(!is.na(dlog_capital)) %>%
    select(name, Country, dlog_capital, all_of(setdiff(vars, "dlog_capital"))) %>%
    drop_na(all_of(setdiff(vars, "dlog_capital")))
  
  if (nrow(tmp) < 5 || length(unique(tmp$dlog_capital)) <= 1) {
    message("-> Insuficientes datos o dlog_capital constante. Modelo omitido.")
    return(NULL)
  }
  
  formula_str <- paste0(
    "dlog_capital ~ ",
    paste(setdiff(vars, "dlog_capital"), collapse = " + "),
    " | name + Country"
  )
  feols(as.formula(formula_str), data = tmp, cluster = c("name", "Country"))
}

# ================================================
# 3) Leer datos
# ================================================
ruta_dta <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.dta"
df_raw <- read_dta(ruta_dta)

# ================================================
# 4) Preprocesamiento: dlog_capital
# ================================================
df <- df_raw %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    ratio_cap    = real_capital / lag(real_capital),
    dlog_capital = suppressWarnings(ifelse(ratio_cap > 0, log(ratio_cap), NA_real_)),
    dlog_capital = ifelse(is.finite(dlog_capital), dlog_capital, NA_real_),
    shock        = -shock
  ) %>%
  select(-ratio_cap) %>%
  ungroup() %>%
  filter(!is.na(dlog_capital))

# ================================================
# 5) Variables macro y rezagos (por Country)
# ================================================
df <- df %>%
  arrange(Country, dateq) %>%
  group_by(Country) %>%
  mutate(
    dlog_gdp = log(gdp) - log(lag(gdp)),
    dlog_cpi = log(cpi) - log(lag(cpi)),
    unemp    = unemp,
    embigl   = embigl,
    across(
      c(dlog_gdp, dlog_cpi, unemp, embigl),
      list(L1 = ~lag(.x,1), L2 = ~lag(.x,2), L3 = ~lag(.x,3), L4 = ~lag(.x,4)),
      .names = "{fn}_{.col}"
    )
  ) %>%
  ungroup()

# ================================================
# 6) Respuestas heterogéneas a nivel firma
# ================================================
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    rsales_g     = log(saleq) - log(lag(saleq)),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    log_atq      = log(atq),
    log_atq_std  = (log_atq - mean(log_atq, na.rm = TRUE)) / sd(log_atq, na.rm = TRUE),
    saleq_std    = (saleq - mean(saleq, na.rm = TRUE)) / sd(saleq, na.rm = TRUE),
    size_index   = ((log_atq_std + saleq_std)/2 - mean((log_atq_std + saleq_std)/2, na.rm = TRUE)) /
      sd((log_atq_std + saleq_std)/2, na.rm = TRUE),
    sh_current_a_std = (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE)
  ) %>%
  select(-log_atq, -log_atq_std, -saleq_std) %>%
  ungroup()

# ================================================
# 7) Winsorizar y crear versiones raw/demeaned
# ================================================
df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# ================================================
# 8) Crear interacciones raw/demeaned con "shock" y controls
# ================================================
df <- df %>%
  mutate(
    lev_raw_shock = lev_dev_raw * shock,
    dd_raw_shock  = dd_dev_raw  * shock,
    lev_std_gdp   = lev_dev_raw * dlog_gdp,  # si deseas raw con GDP
    dd_std_gdp    = dd_dev_raw  * dlog_gdp
  )

# ================================================
# 9) Definir variables para los modelos “raw”
# ================================================
vars_m1_raw <- c("dlog_capital", "lev_raw_shock", "dd_raw_shock")
vars_m2_raw <- c(vars_m1_raw, "lev_std_gdp", "dd_std_gdp")
vars_m3_raw <- c(vars_m2_raw, "rsales_g_std", "size_index", "sh_current_a_std")
vars_m4_raw <- c(
  vars_m3_raw,
  "L1_dlog_gdp", "L1_dlog_cpi", "L1_unemp", "L1_embigl"
)

# ================================================
# 10) Estimar sólo los modelos raw
# ================================================
m1_raw <- run_model(df, vars_m1_raw)
m2_raw <- run_model(df, vars_m2_raw)
m3_raw <- run_model(df, vars_m3_raw)
m4_raw <- run_model(df, vars_m4_raw)

# ================================================
# 11) Presentar resultados sin la etiqueta “raw”
# ================================================
modelsummary(
  list(
    "M1" = m1_raw,
    "M2" = m2_raw,
    "M3" = m3_raw,
    "M4" = m4_raw
  ),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_raw_shock" = "Leverage × Shock",
    "dd_raw_shock"  = "DD × Shock",
    "lev_std_gdp"   = "Leverage × ΔGDP",
    "dd_std_gdp"    = "DD × ΔGDP"
  )
)


# --------------------------------
# Empirical results: Model vs Data 
# --------------------------------

df <- df %>%
  # Asegúrate de tener winsorizadas leverage_win y dd_win
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  mutate(
    # 1) Versiones raw/demeaned
    lev_dev           = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev            = dd_win       - mean(dd_win,       na.rm = TRUE),
    lev_wins_dem_wide = lev_dev      * shock,
    dd_wins_dem_wide  = dd_dev       * shock,
    
    # 2) Versiones estandarizadas
    lev_std           = (leverage_win - mean(leverage_win, na.rm = TRUE)) / sd(leverage_win, na.rm = TRUE),
    dd_std            = (dd_win       - mean(dd_win,       na.rm = TRUE)) / sd(dd_win,       na.rm = TRUE),
    lev_std_shock     = lev_std      * shock,
    dd_std_shock      = dd_std       * shock
  ) %>%
  ungroup()


# 1) Controles a nivel firma (mismos que antes)
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")

# 2) Modelo (1) Std: leverage × shock + DD × shock + controles
m7_1 <- feols(
  dlog_capital ~
    lev_std_gdp +
    lev_std_shock +
    dd_std_shock +
    rsales_g_std + size_index + sh_current_a_std
  | name + Country,
  data    = df,
  cluster = ~name + Country
)

# 3) Modelo (2) Raw: leverage(raw) × shock + DD(raw) × shock + controles
m7_2 <- feols(
  dlog_capital ~
    lev_std_gdp +
    lev_wins_dem_wide +
    dd_wins_dem_wide +
    rsales_g_std + size_index + sh_current_a_std
  | name + Country,
  data    = df,
  cluster = ~name + Country
)

# 4) Presentar Table 7 ampliada con DD
modelsummary(
  list(`(1) Std` = m7_1, `(2) Raw` = m7_2),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_shock"     = "Leverage × Shock",
    "lev_wins_dem_wide" = "Leverage (dem) × Shock",
    "dd_std_shock"      = "DD × Shock",
    "dd_wins_dem_wide"  = "DD (dem) × Shock"
  ),
  gof_map  = c("r2" = "R²"),
  gof_omit = "^(?!r2$)[A-Za-z0-9_.]+",
  title    = "Model vs Data"
)



# --------------------------------
# Table 9: Heterogeneous responses, sin demeaning de posiciones financieras
# --------------------------------

# 0) Asegurarnos de tener winsorizadas las posiciones financieras
df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  )

# 1) Crear dummy “aboveA_dummy” según dd_win
df <- df %>%
  group_by(name) %>%
  mutate(
    A_threshold  = median(dd_win, na.rm = TRUE),
    aboveA_dummy = as.numeric(dd_win >= A_threshold)
  ) %>%
  ungroup()

# 2) Calcular versiones raw/demeaned de leverage y dd
df <- df %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 3) Crear interacciones raw/demeaned con Shock y ΔGDP
df <- df %>%
  mutate(
    lev_raw_wide = lev_dev_raw * shock,
    lev_raw_gdp  = lev_dev_raw * dlog_gdp,
    dd_raw_wide  = dd_dev_raw  * shock,
    dd_raw_gdp   = dd_dev_raw  * dlog_gdp,
    aboveA_wide  = aboveA_dummy * shock,
    aboveA_gdp   = aboveA_dummy * dlog_gdp
  )

# 4) Controles a nivel firma
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")

# 5) Controles agregados (rezagos macro)
controls_aggregate <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 6) Definir especificaciones m1–m7 usando solo raw
vars9_raw <- list(
  m1 = c("dlog_capital", "lev_raw_gdp", "lev_raw_wide"),
  m2 = c("dlog_capital", "lev_raw_gdp", "lev_raw_wide", controls_firm),
  m3 = c("dlog_capital", "aboveA_gdp",   "aboveA_wide",  controls_firm),
  m4 = c("dlog_capital", "dd_raw_gdp",  "dd_raw_wide",  controls_firm),
  m5 = c("dlog_capital",
         "lev_raw_gdp", "aboveA_gdp",
         "lev_raw_wide","aboveA_wide",
         controls_firm),
  m6 = c("dlog_capital",
         "lev_raw_gdp", "dd_raw_gdp",
         "lev_raw_wide","dd_raw_wide",
         controls_firm),
  m7 = c("dlog_capital", "shock",
         "lev_raw_gdp",  "dd_raw_gdp",
         "lev_raw_wide", "dd_raw_wide",
         "aboveA_wide","aboveA_gdp",
         controls_firm, controls_aggregate)
)

# 7) Estimar m1–m7 con lapply
models9_raw <- lapply(names(vars9_raw), function(m) {
  vars_i <- vars9_raw[[m]]
  df_i   <- df %>%
    select(name, Country, all_of(vars_i)) %>%
    drop_na()
  fmla   <- as.formula(
    paste("dlog_capital ~",
          paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
          "| name + Country")
  )
  feols(fmla, data = df_i, cluster = ~name + Country)
})
names(models9_raw) <- names(vars9_raw)

# 8) Table 9 raw con R², AIC y BIC
modelsummary(
  models9_raw,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_raw_wide" = "Leverage (wins, dem) × Shock",
    "dd_raw_wide"  = "DD (wins, dem) × Shock",
    "aboveA_wide"  = "Above-A × Shock",
    "shock"        = "Shock (puro)"
  ),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = NULL,
  title    = "Heterogeneous Responses (raw, no demeaning)"
)



# --------------------------------
# Table 9 (con demeaning): Heterogeneous responses, con demeaning de posiciones financieras
# --------------------------------

# 0) Asegúrate de tener winsorizadas las posiciones financieras
df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  )

# 1) Dummy “aboveA_dummy” según dd_win
df <- df %>%
  group_by(name) %>%
  mutate(
    A_threshold    = median(dd_win, na.rm = TRUE),
    aboveA_dummy   = as.numeric(dd_win >= A_threshold)
  ) %>%
  ungroup()

# 2) Demeaning de leverage y dd (raw / demeaned)
df <- df %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 3) Crear interacciones raw/demeaned con Shock y ΔGDP
df <- df %>%
  mutate(
    # leverage raw × shock / GDP
    lev_raw_wide = lev_dev_raw * shock,
    lev_raw_gdp  = lev_dev_raw * dlog_gdp,
    # dd raw × shock / GDP
    dd_raw_wide  = dd_dev_raw  * shock,
    dd_raw_gdp   = dd_dev_raw  * dlog_gdp,
    # above‐A dummy × shock / GDP
    aboveA_wide  = aboveA_dummy * shock,
    aboveA_gdp   = aboveA_dummy * dlog_gdp
  )

# 4) Controles a nivel firma
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")

# 5) Controles agregados (rezagos macro)
controls_aggregate <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 6) Definición de especificaciones m1–m7 (usando raw)
vars9_raw <- list(
  m1 = c("dlog_capital", "lev_raw_gdp", "lev_raw_wide"),
  m2 = c("dlog_capital", "lev_raw_gdp", "lev_raw_wide", controls_firm),
  m3 = c("dlog_capital", "aboveA_gdp",   "aboveA_wide",  controls_firm),
  m4 = c("dlog_capital", "dd_raw_gdp",  "dd_raw_wide",  controls_firm),
  m5 = c("dlog_capital",
         "lev_raw_gdp", "aboveA_gdp",
         "lev_raw_wide","aboveA_wide",
         controls_firm),
  m6 = c("dlog_capital",
         "lev_raw_gdp", "dd_raw_gdp",
         "lev_raw_wide","dd_raw_wide",
         controls_firm),
  m7 = c("dlog_capital", "shock",
         "lev_raw_gdp", "dd_raw_gdp",
         "lev_raw_wide","dd_raw_wide",
         "aboveA_wide","aboveA_gdp",
         controls_firm, controls_aggregate)
)

# 7) Estimar m1–m7 con lapply
models9_raw <- lapply(names(vars9_raw), function(m) {
  vars_i <- vars9_raw[[m]]
  df_i   <- df %>% select(name, Country, all_of(vars_i)) %>% drop_na()
  fmla   <- as.formula(
    paste("dlog_capital ~",
          paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
          "| name + Country")
  )
  feols(fmla, data = df_i, cluster = ~name + Country)
})
names(models9_raw) <- names(vars9_raw)

# 8) Table 9 raw con R², AIC y BIC
modelsummary(
  models9_raw,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_raw_wide" = "Leverage (wins, dem) × Shock",
    "dd_raw_wide"  = "DD (wins, dem) × Shock",
    "aboveA_wide"  = "Above-A × Shock",
    "shock"        = "Shock (puro)"
  )
)




# =======================================
# Table 10: Resultados principales, sin sensibilidades cíclicas
#            (solo shock + FE firma + FE país + EMBIGL)
#  — usando únicamente interacciones raw/demeaned —
# =======================================

# 1) Asegurar controls_firm
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
  message("Aviso: 'controls_firm' no existía y se creó aquí.")
}

# 2) Asegurar controls_agg
if (!exists("controls_agg")) {
  if (exists("controls_aggregate")) {
    controls_agg <- controls_aggregate
    message("Aviso: 'controls_agg' asignado desde 'controls_aggregate'.")
  } else {
    controls_agg <- c(
      paste0("L", 1:4, "_dlog_gdp"),
      paste0("L", 1:4, "_dlog_cpi"),
      paste0("L", 1:4, "_unemp"),
      paste0("L", 1:4, "_embigl")
    )
    message("Aviso: 'controls_agg' no existía y se creó aquí.")
  }
}

# 3) Crear interacciones raw/demeaned con shock
#    (lev_dev_raw y dd_dev_raw deben existir de pasos previos)
df <- df %>%
  mutate(
    lev_raw_shock = lev_dev_raw * shock,
    dd_raw_shock  = dd_dev_raw  * shock
  )

# 4) Especificaciones para m1–m4 usando raw
vars10_raw <- list(
  m1 = c("dlog_capital", "lev_raw_shock", controls_firm),
  m2 = c("dlog_capital", "dd_raw_shock",  controls_firm),
  m3 = c("dlog_capital", "lev_raw_shock", "dd_raw_shock", controls_firm),
  m4 = c("dlog_capital",
         "shock",
         "lev_raw_shock", "dd_raw_shock",
         controls_firm, controls_agg)
)

# 5) Estimar m1–m4 usando lapply
models10_raw <- lapply(names(vars10_raw), function(m) {
  vars_i <- vars10_raw[[m]]
  df_i   <- df %>%
    select(name, Country, all_of(vars_i)) %>%
    drop_na()
  f_i <- as.formula(
    paste("dlog_capital ~",
          paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
          "| name + Country")
  )
  feols(f_i, data = df_i, cluster = ~name + Country)
})
names(models10_raw) <- names(vars10_raw)

# 6) Mostrar Table 10 con R², AIC y BIC (sin etiquetar 'raw')
modelsummary(
  models10_raw,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_raw_shock" = "Leverage × Shock",
    "dd_raw_shock"  = "DD × Shock",
    "shock"         = "Shock (puro)"
  ),
)





# ----------------------------------------
# Table 11: regresiones por separado para shocks positivos y negativos
#            (usando únicamente interacciones raw/demeaned)
# ----------------------------------------

# 0) Asegurarse de tener winsorizadas las posiciones financieras y dev_raw
df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 1) Crear variables de shock asimétrico
df <- df %>%
  mutate(
    shock_p = if_else(shock > 0, shock, 0),
    shock_n = if_else(shock < 0, shock, 0)
  )

# 2) Normalizar shock_p y shock_n (z‐scores globales)
df <- df %>%
  mutate(
    shock_p_z = (shock_p - mean(shock_p, na.rm = TRUE)) / sd(shock_p, na.rm = TRUE),
    shock_n_z = (shock_n - mean(shock_n, na.rm = TRUE)) / sd(shock_n, na.rm = TRUE)
  )

# 3) Crear interacciones raw/demeaned con los shocks normalizados
df <- df %>%
  mutate(
    lev_raw_shock_p = lev_dev_raw * shock_p_z,
    dd_raw_shock_p  = dd_dev_raw  * shock_p_z,
    lev_raw_shock_n = lev_dev_raw * shock_n_z,
    dd_raw_shock_n  = dd_dev_raw  * shock_n_z
  )

# 4) Definir controles de firma y controles agregados
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_agg  <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 5) Especificación única con ambos shocks raw
m11_raw <- feols(
  dlog_capital ~ 
    lev_raw_shock_p + dd_raw_shock_p +    # interacción con shocks positivos
    lev_raw_shock_n + dd_raw_shock_n +    # interacción con shocks negativos
    rsales_g_std + size_index + sh_current_a_std  # controles de firma
  | name + Country,                       # efectos fijos
  data    = df,
  cluster = ~ name + Country              # clustering
)

# 6) Mostrar resultados sin “raw” en las etiquetas
modelsummary(
  m11_raw,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    lev_raw_shock_p = "Leverage × Shock (+)",
    dd_raw_shock_p  = "DD       × Shock (+)",
    lev_raw_shock_n = "Leverage × Shock (–)",
    dd_raw_shock_n  = "DD       × Shock (–)"
  )
)








# ================================================
# PARA ANEXOS:
# ================================================


# ================================================
# Table 20: Alternative Time Aggregation
#            (solo “shock_sum” + FE firma + FE país + EMBIGL)
#            — usando únicamente interacciones raw & winsorize —
# ================================================

# 0) Winsorizar y crear versiones raw/demeaned de leverage y dd
df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 1) Crear rezagos de shock y sumar
if (!"L0_shock" %in% names(df)) {
  df <- df %>%
    arrange(name, dateq) %>%
    group_by(name) %>%
    mutate(
      L0_shock  = shock,
      L1_shock  = lag(shock, 1),
      L2_shock  = lag(shock, 2),
      L3_shock  = lag(shock, 3),
      L4_shock  = lag(shock, 4),
      shock_sum = L0_shock + L1_shock + L2_shock + L3_shock + L4_shock
    ) %>%
    ungroup()
  message("Aviso: rezagos L0–L4 de shock y shock_sum creados aquí.")
}

# 2) Normalizar shock_sum (z‐score global)
if (!"shock_sum_z" %in% names(df)) {
  df <- df %>%
    mutate(
      shock_sum_z = (shock_sum - mean(shock_sum, na.rm = TRUE)) /
        sd(shock_sum, na.rm = TRUE)
    )
  message("Aviso: 'shock_sum_z' creado aquí.")
}

# 3) Crear interacciones raw/demeaned con shock_sum_z y ciclo raw
df <- df %>%
  mutate(
    lev_raw_shock_sum = lev_dev_raw * shock_sum_z,
    dd_raw_shock_sum  = dd_dev_raw  * shock_sum_z,
    lev_raw_gdp       = lev_dev_raw * dlog_gdp,
    dd_raw_gdp        = dd_dev_raw  * dlog_gdp
  )

# 4) Asegurar controles a nivel firma y agregados
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
  message("Aviso: 'controls_firm' creada aquí.")
}
if (!exists("controls_agg")) {
  controls_agg <- c(
    paste0("L", 1:4, "_dlog_gdp"),
    paste0("L", 1:4, "_dlog_cpi"),
    paste0("L", 1:4, "_unemp"),
    paste0("L", 1:4, "_embigl")
  )
  message("Aviso: 'controls_agg' creada aquí.")
}

# 5) Estimar modelos m20_1–m20_4 usando raw
# M1: canal apalancamiento + ciclo y suma
f20_1 <- as.formula(paste0(
  "dlog_capital ~ lev_raw_gdp + lev_raw_shock_sum + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
m20_1 <- feols(f20_1, data = df %>% drop_na(all.vars(f20_1)[-1]),
               cluster = ~name + Country)

# M2: canal solvencia + ciclo y suma
f20_2 <- as.formula(paste0(
  "dlog_capital ~ dd_raw_gdp + dd_raw_shock_sum + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
m20_2 <- feols(f20_2, data = df %>% drop_na(all.vars(f20_2)[-1]),
               cluster = ~name + Country)

# M3: ambos canales + ciclo y suma
f20_3 <- as.formula(paste0(
  "dlog_capital ~ lev_raw_gdp + dd_raw_gdp + ",
  "lev_raw_shock_sum + dd_raw_shock_sum + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
m20_3 <- feols(f20_3, data = df %>% drop_na(all.vars(f20_3)[-1]),
               cluster = ~name + Country)

# M4: shock_sum puro + canales + controles agregados
rhs20_4 <- c(
  "lev_raw_gdp", "shock_sum_z", "lev_raw_shock_sum", "dd_raw_shock_sum",
  controls_firm, "embigl", controls_agg
)
f20_4 <- as.formula(
  paste("dlog_capital ~", paste(rhs20_4, collapse = " + "), "| name + Country")
)
m20_4 <- feols(f20_4, data = df %>% drop_na(all_of(rhs20_4)),
               cluster = ~name + Country)

# 6) Mostrar Table 20 con etiquetas limpias
modelsummary(
  list(M1 = m20_1, M2 = m20_2, M3 = m20_3, M4 = m20_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    lev_raw_shock_sum = "Leverage × shock (sum)",
    dd_raw_shock_sum  = "DD       × shock (sum)",
    shock_sum_z       = "Sum of shocks (z-score)",
    lev_raw_gdp       = "Leverage × ΔGDP",
    dd_raw_gdp        = "DD       × ΔGDP"
  ),
  keep     = c("lev_raw_shock_sum", "dd_raw_shock_sum", "shock_sum_z"),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 20: Alternative Time Aggregation"
)



# ================================================
# Table 21: Interaction with Other Firm‐Level Covariates (con winsorización & raw/demeaned para leverage y DD)
# ================================================

# 0) Winsorizar (1%) los shocks y las variables financieras
winsorize <- function(x, p = 0.01) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

df <- df %>%
  # winsorize shock, leverage y dd
  mutate(
    shock_w       = winsorize(shock),
    leverage_win  = winsorize(leverage),
    dd_win        = winsorize(dd)
  )

# 1) Calcular crecimiento futuro de ventas (4 trimestres adelante) y winsorizar + estandarizar
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    rsales_Fg4_raw   = log(dplyr::lead(saleq, 4)) - log(dplyr::lead(saleq, 3)),
    rsales_Fg4_w     = winsorize(rsales_Fg4_raw)
  ) %>%
  ungroup() %>%
  mutate(
    rsales_Fg4_std_w = (rsales_Fg4_w - mean(rsales_Fg4_w, na.rm = TRUE)) /
      sd(rsales_Fg4_w, na.rm = TRUE)
  )

# 2) Winsorizar y re‐estandarizar las demás covariables (global)
df <- df %>%
  mutate(
    rsales_g_std_w      = winsorize(rsales_g_std),
    size_index_w        = winsorize(size_index),
    sh_current_a_std_w  = winsorize(sh_current_a_std)
  ) %>%
  mutate(
    rsales_g_std_win      = (rsales_g_std_w      - mean(rsales_g_std_w,      na.rm = TRUE)) /
      sd(rsales_g_std_w,      na.rm = TRUE),
    size_index_win        = (size_index_w        - mean(size_index_w,        na.rm = TRUE)) /
      sd(size_index_w,        na.rm = TRUE),
    sh_current_a_std_win  = (sh_current_a_std_w  - mean(sh_current_a_std_w,   na.rm = TRUE)) /
      sd(sh_current_a_std_w,   na.rm = TRUE)
  )

# 3) Calcular raw deviations para leverage y dd
df <- df %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 4) Crear interacciones con el choque winsorizado
df <- df %>%
  mutate(
    # leverage y DD raw × shock
    lev_raw_shock          = lev_dev_raw            * shock_w,
    dd_raw_shock           = dd_dev_raw             * shock_w,
    # otras covariables estandarizadas × shock
    rsales_g_std_shock     = rsales_g_std_win       * shock_w,
    rsales_Fg4_std_shock   = rsales_Fg4_std_w       * shock_w,
    size_index_shock       = size_index_win         * shock_w,
    sh_current_a_std_shock = sh_current_a_std_win   * shock_w
  )

# 5) Definir controles
controls_firm <- c("rsales_g_std_win", "size_index_win", "sh_current_a_std_win")
controls_agg  <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 6) Especificaciones m1–m8 usando raw interactions para leverage y DD
vars21_raw <- list(
  m1 = c("dlog_capital",
         "lev_raw_shock",         "rsales_g_std_shock",
         controls_firm,           "embigl"),
  m2 = c("dlog_capital",
         "dd_raw_shock",          "rsales_g_std_shock",
         controls_firm,           "embigl"),
  m3 = c("dlog_capital",
         "lev_raw_shock",         "rsales_Fg4_std_shock",
         "rsales_Fg4_std_w",      controls_firm, "embigl"),
  m4 = c("dlog_capital",
         "dd_raw_shock",          "rsales_Fg4_std_shock",
         "rsales_Fg4_std_w",      controls_firm, "embigl"),
  m5 = c("dlog_capital",
         "lev_raw_shock",         "size_index_shock",
         controls_firm,           "embigl"),
  m6 = c("dlog_capital",
         "dd_raw_shock",          "size_index_shock",
         controls_firm,           "embigl"),
  m7 = c("dlog_capital",
         "lev_raw_shock",         "sh_current_a_std_shock",
         "sh_current_a_std_w",    controls_firm, "embigl"),
  m8 = c("dlog_capital",
         "dd_raw_shock",          "sh_current_a_std_shock",
         "sh_current_a_std_w",    controls_firm, "embigl")
)

# 7) Estimar m1–m8 con feols()
models21_raw <- lapply(vars21_raw, function(vars) {
  df_sub <- df %>%
    select(name, Country, all_of(vars)) %>%
    drop_na()
  fmla <- as.formula(
    paste("dlog_capital ~",
          paste(setdiff(vars, "dlog_capital"), collapse = " + "),
          "| name + Country")
  )
  feols(fmla, data = df_sub, cluster = ~name + Country)
})

# 8) Mostrar Table 21 con etiquetas limpias
modelsummary(
  models21_raw,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_raw_shock", "dd_raw_shock",
               "rsales_g_std_shock", "rsales_Fg4_std_shock",
               "size_index_shock", "sh_current_a_std_shock"),
  order    = c("lev_raw_shock", "dd_raw_shock",
               "rsales_g_std_shock", "rsales_Fg4_std_shock",
               "size_index_shock", "sh_current_a_std_shock"),
  coef_map = c(
    lev_raw_shock          = "Leverage × shock",
    dd_raw_shock           = "DD       × shock",
    rsales_g_std_shock     = "Sales growth × shock",
    rsales_Fg4_std_shock   = "Future sales growth × shock",
    size_index_shock       = "Size     × shock",
    sh_current_a_std_shock = "Liquidity× shock"
  ),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  title    = "Table 21: Interaction with Other Firm‐Level Covariates"
)




# ===============================================================
# Table 22: Interacción con Otras Medidas de Posición Financiera
#            (solo “shock” + FE firma + FE país + EMBIGL)
#            — usando interacciones raw/demeaned para leverage y DD —
# ===============================================================

# 0) Winsorizar shock y posiciones financieras
winsorize <- function(x, p = 0.01) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

df <- df %>%
  mutate(
    shock_w      = winsorize(shock),
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  )

# 1) Calcular raw deviations de leverage y dd para interacción
df <- df %>%
  group_by(name) %>%
  mutate(
    lev_dev_raw = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev_raw  = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 2) Estandarizar otras medidas de posición financiera (si faltan)
df <- df %>%
  mutate(
    cashop_win = winsorize(cashop),
    cashop_std = (cashop_win - mean(cashop_win, na.rm = TRUE)) / sd(cashop_win, na.rm = TRUE),
    div_win    = winsorize(dividends),
    div_std    = (div_win    - mean(div_win,    na.rm = TRUE)) / sd(div_win,    na.rm = TRUE),
    liq_win    = winsorize(sh_current_a_std),
    liq_std    = (liq_win    - mean(liq_win,    na.rm = TRUE)) / sd(liq_win,    na.rm = TRUE),
    size_win   = winsorize(size_index),
    size_std   = (size_win   - mean(size_win,   na.rm = TRUE)) / sd(size_win,   na.rm = TRUE)
  )

# 3) Crear interacciones con shock_w
df <- df %>%
  mutate(
    lev_raw_shock       = lev_dev_raw * shock_w,
    dd_raw_shock        = dd_dev_raw  * shock_w,
    size_std_shock      = size_std     * shock_w,
    cashop_std_shock    = cashop_std   * shock_w,
    div_std_shock       = div_std      * shock_w,
    liq_std_shock       = liq_std      * shock_w
  )

# 4) Asegurar controles de firma y agregados
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
}
if (!exists("controls_agg")) {
  controls_agg <- if (exists("controls_aggregate")) controls_aggregate else c(
    paste0("L",1:4,"_dlog_gdp"),
    paste0("L",1:4,"_dlog_cpi"),
    paste0("L",1:4,"_unemp"),
    paste0("L",1:4,"_embigl")
  )
}

# 5) Definir variables para m1–m8
vars22_raw <- list(
  m1 = c("dlog_capital", "lev_raw_shock",    "size_std_shock",    controls_firm, "embigl"),
  m2 = c("dlog_capital", "dd_raw_shock",     "size_std_shock",    controls_firm, "embigl"),
  m3 = c("dlog_capital", "lev_raw_shock",    "cashop_std_shock",  "cashop_std",  controls_firm, "embigl"),
  m4 = c("dlog_capital", "dd_raw_shock",     "cashop_std_shock",  "cashop_std",  controls_firm, "embigl"),
  m5 = c("dlog_capital", "lev_raw_shock",    "div_std_shock",     "div_std",     controls_firm, "embigl"),
  m6 = c("dlog_capital", "dd_raw_shock",     "div_std_shock",     "div_std",     controls_firm, "embigl"),
  m7 = c("dlog_capital", "lev_raw_shock",    "liq_std_shock",     "liq_std",     controls_firm, "embigl"),
  m8 = c("dlog_capital", "dd_raw_shock",     "liq_std_shock",     "liq_std",     controls_firm, "embigl")
)

# 6) Estimar modelos m1–m8
models22_raw <- lapply(vars22_raw, function(v) {
  df_sub <- df %>% select(name, Country, all_of(v)) %>% drop_na()
  f      <- as.formula(paste(
    "dlog_capital ~",
    paste(setdiff(v, "dlog_capital"), collapse = " + "),
    "| name + Country"
  ))
  feols(f, data = df_sub, cluster = ~name + Country)
})

# 7) Mostrar Table 22 (solo interacciones, R², AIC, BIC)
modelsummary(
  models22_raw,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_raw_shock", "dd_raw_shock",
               "size_std_shock","cashop_std_shock",
               "div_std_shock", "liq_std_shock"),
  order    = c("lev_raw_shock", "dd_raw_shock",
               "size_std_shock","cashop_std_shock",
               "div_std_shock", "liq_std_shock"),
  coef_map = c(
    "lev_raw_shock"    = "Leverage × shock",
    "dd_raw_shock"     = "DD       × shock",
    "size_std_shock"   = "Size     × shock",
    "cashop_std_shock" = "Cash flow × shock",
    "div_std_shock"    = "Dividend  × shock",
    "liq_std_shock"    = "Liquidity × shock"
  ),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 22: Interaction with Other Measures of Financial Positions"
)
