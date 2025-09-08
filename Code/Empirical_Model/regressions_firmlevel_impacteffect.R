
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
library(dplyr)        # manipulación de datos
library(zoo)          # as.yearqtr()
library(fixest)      # feols()
library(modelsummary) # presentación de tablas
library(tidyr)        # drop_na()
library(tibble)       # tibble
library(writexl)      # para exportar
#library(plm)  # Me da error, importalo despues de crear dlog_capital

# ================================================
# 2) Definir funciones auxiliares
# ================================================
# Winsorizar valores en percentiles p y (1-p)
winsorize <- function(x, p = 0.005) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

# Función para estimar modelo de efectos fijos
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
    ratio_cap     = real_capital / lag(real_capital),
    dlog_capital  = suppressWarnings(ifelse(ratio_cap > 0, log(ratio_cap), NA_real_)),
    shock = -shock,
  ) %>%
  select(-ratio_cap) %>%
  ungroup()


# ================================================
# 5) Variables macro y rezagos (por Country)
# ================================================
df <- df %>%
  arrange(Country, dateq) %>%
  group_by(Country) %>%
  mutate(
    dlog_gdp   = log(gdp) - log(lag(gdp)),
    dlog_cpi   = log(cpi) - log(lag(cpi)),
    unemp      = unemp,
    embigl     = embigl,
    across(c(dlog_gdp, dlog_cpi, unemp, embigl), list(
      L1 = ~lag(.x, 1),
      L2 = ~lag(.x, 2),
      L3 = ~lag(.x, 3),
      L4 = ~lag(.x, 4)
    ), .names = "{fn}_{.col}")
  ) %>%
  ungroup()

# ================================================
# 6) Heterogeneous Responses: variables a nivel firma
# ================================================
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    # 6.1) Crecimiento de ventas trimestral estandarizado
    rsales_g     = log(saleq) - log(lag(saleq)),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 6.2) Proxy compuesto de tamaño (Opción B, más robusto)
    log_atq      = log(atq),
    log_atq_std  = (log_atq - mean(log_atq, na.rm = TRUE)) / sd(log_atq, na.rm = TRUE),
    saleq_std    = (saleq - mean(saleq, na.rm = TRUE)) / sd(saleq, na.rm = TRUE),
    size_index   = (log_atq_std + saleq_std) / 2,
    size_index   = (size_index - mean(size_index, na.rm = TRUE)) / sd(size_index, na.rm = TRUE),
    
    # 6.3) Liquidez corriente estandarizada
    sh_current_a_std = (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE)
  ) %>%
  select(-log_atq, -log_atq_std, -saleq_std) %>%
  ungroup()

# ================================================
# 7) Winsorizar y estandarizar variables financieras
# ================================================
df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  mutate(
    lev_std = (leverage_win - mean(leverage_win, na.rm = TRUE)) / sd(leverage_win, na.rm = TRUE),
    dd_std  = (dd_win       - mean(dd_win, na.rm = TRUE))       / sd(dd_win, na.rm = TRUE)
  ) %>%
  ungroup()

# ================================================
# 8) Crear interacciones con "shock" y dlog_gdp
# ================================================
df <- df %>%
  mutate(
    lev_std_shock = lev_std * shock,
    dd_std_shock  = dd_std  * shock,
    lev_std_gdp   = lev_std * dlog_gdp,
    dd_std_gdp    = dd_std  * dlog_gdp
  )

# ================================================
# 9) Definir variables para los modelos M1–M4
# ================================================
vars_m1 <- c("dlog_capital", "lev_std_shock", "dd_std_shock")
vars_m2 <- c(vars_m1, "lev_std_gdp", "dd_std_gdp")
vars_m3 <- c(vars_m2, "rsales_g_std", "size_index", "sh_current_a_std")
vars_m4 <- c(
  vars_m3,
  paste0("L", 1:1, "_dlog_gdp"),
  paste0("L", 1:1, "_dlog_cpi"),
  paste0("L", 1:1, "_unemp"),
  paste0("L", 1:1, "_embigl")
)

# ================================================
# 10) Estimar modelos
# ================================================
m1 <- run_model(df, vars_m1)
m2 <- run_model(df, vars_m2)
m3 <- run_model(df, vars_m3)
m4 <- run_model(df, vars_m4)

# ================================================
# 11) Presentar resultados
# ================================================
modelsummary(
  list(M1 = m1, M2 = m2, M3 = m3, M4 = m4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_shock" = "Leverage × Shock",
    "dd_std_shock"  = "DD × Shock",
    "lev_std_gdp"   = "Leverage × ΔGDP",
    "dd_std_gdp"    = "DD × ΔGDP"
  )
)





# --------------------------------
# Empirical results: Model vs Data 
# --------------------------------

# 0) Crear interacciones “std” y “raw” para Leverage y DD
df <- df %>%
  group_by(name) %>%
  mutate(
    # Desviaciones “raw” dentro de firma
    lev_dev           = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev            = dd_win       - mean(dd_win,       na.rm = TRUE),
    # Interacciones “raw” con el shock
    lev_wins_dem_wide = lev_dev      * shock,
    dd_wins_dem_wide  = dd_dev       * shock
  ) %>%
  ungroup() %>%
  mutate(
    # Interacciones “estándar” con el shock
    lev_std_shock     = lev_std      * shock,
    dd_std_shock      = dd_std       * shock
  )

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

# 1) Dummy “aboveA_dummy” según dd_win
df <- df %>%
  group_by(name) %>%
  mutate(
    A_threshold  = median(dd_win, na.rm = TRUE),
    aboveA_dummy = as.numeric(dd_win >= A_threshold)
  ) %>%
  ungroup()

# 2) Interacciones “no-demeaning” ('nodem,std') con shock y ΔGDP
df <- df %>%
  mutate(
    lev_wins_nodem_std_wide = (leverage_win / sd(leverage_win, na.rm = TRUE)) * shock,
    lev_wins_nodem_std_gdp  = (leverage_win / sd(leverage_win, na.rm = TRUE)) * dlog_gdp,
    d2d_wins_nodem_std_wide = (dd_win      / sd(dd_win,      na.rm = TRUE)) * shock,
    d2d_wins_nodem_std_gdp  = (dd_win      / sd(dd_win,      na.rm = TRUE)) * dlog_gdp,
    aboveA_dummy_wide       = aboveA_dummy * shock,
    aboveA_dummy_gdp        = aboveA_dummy * dlog_gdp
  )

# 3) Controles a nivel firma (definidos previamente)
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")

# 4) Controles agregados (rezagos macro)
controls_aggregate <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 5) Definición de especificaciones m1–m7
vars9 <- list(
  m1 = c("dlog_capital", "lev_wins_nodem_std_gdp", "lev_wins_nodem_std_wide"),
  m2 = c("dlog_capital", "lev_wins_nodem_std_gdp", "lev_wins_nodem_std_wide",
         controls_firm),
  m3 = c("dlog_capital", "aboveA_dummy_gdp", "aboveA_dummy_wide",
         controls_firm),
  m4 = c("dlog_capital", "d2d_wins_nodem_std_gdp", "d2d_wins_nodem_std_wide",
         controls_firm),
  m5 = c("dlog_capital", "lev_wins_nodem_std_gdp", "aboveA_dummy_gdp",
         "lev_wins_nodem_std_wide", "aboveA_dummy_wide",
         controls_firm),
  m6 = c("dlog_capital", "lev_wins_nodem_std_gdp", "d2d_wins_nodem_std_gdp",
         "lev_wins_nodem_std_wide", "d2d_wins_nodem_std_wide",
         controls_firm),
  m7 = c("dlog_capital", "shock",
         "lev_wins_nodem_std_gdp", "d2d_wins_nodem_std_gdp",
         "lev_wins_nodem_std_wide", "d2d_wins_nodem_std_wide",
         "aboveA_dummy_wide", "aboveA_dummy_gdp",
         controls_firm, controls_aggregate)
)

# 6) Estimar m1–m7 con lapply
models9 <- lapply(names(vars9), function(m) {
  vars_i <- vars9[[m]]
  df_i   <- df %>%
    select(name, Country, all_of(vars_i)) %>%
    drop_na()
  fmla_i <- as.formula(paste(
    "dlog_capital ~",
    paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
    "| name + Country"
  ))
  feols(fmla_i, data = df_i, cluster = ~name + Country)
})
names(models9) <- names(vars9)

# 7) Table 9 con R², AIC y BIC
modelsummary(
  models9,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_wins_nodem_std_wide" = "Leverage (wins,nodem,std) × Shock",
    "d2d_wins_nodem_std_wide" = "DD (wins,nodem,std) × Shock",
    "aboveA_dummy_wide"       = "Above-A × Shock",
    "shock"                   = "Shock (puro)"
  ),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = NULL,
  title    = "Heterogeneous Responses (No demeaning financial position)"
)




# --------------------------------
# Table 9 (con demeaning): Heterogeneous responses, con demeaning de posiciones financieras
# --------------------------------

# 1) Dummy “aboveA_dummy” según dd_win (igual que antes)
df <- df %>%
  group_by(name) %>%
  mutate(
    A_threshold  = median(dd_win, na.rm = TRUE),
    aboveA_dummy = as.numeric(dd_win >= A_threshold)
  ) %>%
  ungroup()

# 2) Demeaning de leverage y dd (restar media por firma)
df <- df %>%
  group_by(name) %>%
  mutate(
    lev_dev      = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev       = dd_win       - mean(dd_win,       na.rm = TRUE)
  ) %>%
  ungroup()

# 3) Crear interacciones “dem” (demeaned & std) con Shock y ΔGDP
df <- df %>%
  mutate(
    lev_dev_std_wide = (lev_dev / sd(leverage_win, na.rm = TRUE)) * shock,
    lev_dev_std_gdp  = (lev_dev / sd(leverage_win, na.rm = TRUE)) * dlog_gdp,
    dd_dev_std_wide  = (dd_dev  / sd(dd_win,      na.rm = TRUE)) * shock,
    dd_dev_std_gdp   = (dd_dev  / sd(dd_win,      na.rm = TRUE)) * dlog_gdp,
    aboveA_dummy_wide = aboveA_dummy * shock,
    aboveA_dummy_gdp  = aboveA_dummy * dlog_gdp
  )

# 4) Controles a nivel firma (definidos previamente)
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")

# 5) Controles agregados (rezagos macro)
controls_aggregate <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 6) Definición de especificaciones m1–m7 (usando las interacciones “dem”)
vars9_dem <- list(
  m1 = c("dlog_capital", "lev_dev_std_gdp", "lev_dev_std_wide"),
  m2 = c("dlog_capital", "lev_dev_std_gdp", "lev_dev_std_wide", controls_firm),
  m3 = c("dlog_capital", "aboveA_dummy_gdp", "aboveA_dummy_wide", controls_firm),
  m4 = c("dlog_capital", "dd_dev_std_gdp", "dd_dev_std_wide", controls_firm),
  m5 = c("dlog_capital",
         "lev_dev_std_gdp", "aboveA_dummy_gdp",
         "lev_dev_std_wide", "aboveA_dummy_wide",
         controls_firm),
  m6 = c("dlog_capital",
         "lev_dev_std_gdp", "dd_dev_std_gdp",
         "lev_dev_std_wide", "dd_dev_std_wide",
         controls_firm),
  m7 = c("dlog_capital", "shock",
         "lev_dev_std_gdp", "dd_dev_std_gdp",
         "lev_dev_std_wide", "dd_dev_std_wide",
         "aboveA_dummy_wide", "aboveA_dummy_gdp",
         controls_firm, controls_aggregate)
)

# 7) Estimar m1–m7 con lapply
models9_dem <- lapply(names(vars9_dem), function(m) {
  vars_i <- vars9_dem[[m]]
  df_i   <- df %>%
    select(name, Country, all_of(vars_i)) %>%
    drop_na()
  fmla_i <- as.formula(
    paste("dlog_capital ~",
          paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
          "| name + Country")
  )
  feols(fmla_i, data = df_i, cluster = ~name + Country)
})
names(models9_dem) <- names(vars9_dem)

# 8) Table 9 dem con R², AIC y BIC
modelsummary(
  models9_dem,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_dev_std_wide" = "Leverage (wins, dem,std) × Shock",
    "dd_dev_std_wide"  = "DD (wins, dem,std) × Shock",
    "aboveA_dummy_wide"= "Above-A × Shock",
    "shock"            = "Shock (puro)"
  )
)






# =======================================
# Table 10: Resultados principales, sin sensibilidades cíclicas
#            (solo shock + FE firma + FE país + EMBIGL)
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

# 3) Interacciones estándar “demeaning & std” con shock
#    (lev_std y dd_std ya existen)
df <- df %>%
  mutate(
    lev_std_shock = lev_std * shock,
    dd_std_shock  = dd_std  * shock
  )

# 4) Especificaciones para m1–m4
vars10 <- list(
  m1 = c("dlog_capital",
         "lev_std_shock",
         controls_firm),
  m2 = c("dlog_capital",
         "dd_std_shock",
         controls_firm),
  m3 = c("dlog_capital",
         "lev_std_shock", "dd_std_shock",
         controls_firm),
  m4 = c("dlog_capital",
         "shock",
         "lev_std_shock", "dd_std_shock",
         controls_firm,
         controls_agg)
)

# 5) Estimar m1–m4 usando lapply
models10 <- lapply(names(vars10), function(m) {
  vars_i <- vars10[[m]]
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
names(models10) <- names(vars10)

# 6) Mostrar Table 10 con R², AIC y BIC
modelsummary(
  models10,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_shock" = "Leverage × Shock",
    "dd_std_shock"  = "DD × Shock",
    "shock"         = "Shock (puro)"
  )
)




# ----------------------------------------
# Table 11: regresiones por separado para shocks positivos y negativos
# ----------------------------------------

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

# 3) Crear interacciones “demeaning & std” con los shocks normalizados
df <- df %>%
  mutate(
    lev_std_shock_p = lev_std * shock_p_z,
    dd_std_shock_p  = dd_std  * shock_p_z,
    lev_std_shock_n = lev_std * shock_n_z,
    dd_std_shock_n  = dd_std  * shock_n_z
  )

# 4) Definir controles de firma y controles agregados
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_agg  <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 5) Estimar dos regresiones por separado

# 5.1) Shock expansivo
m_shock_p <- feols(
  dlog_capital ~ 
    lev_std_shock_p + dd_std_shock_p +
    rsales_g_std + size_index + sh_current_a_std
  | name + Country,
  data    = df,
  cluster = ~name + Country
)

# 5.2) Shock contractivo
m_shock_n <- feols(
  dlog_capital ~ 
    lev_std_shock_n + dd_std_shock_n +
    rsales_g_std + size_index + sh_current_a_std
  | name + Country,
  data    = df,
  cluster = ~name + Country
)

# 6) Presentar resultados lado a lado
modelsummary(
  list(
    "Shock Expansivo"   = m_shock_p,
    "Shock Contractivo" = m_shock_n
  ),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    lev_std_shock_p = "Leverage × Shock (+)",
    dd_std_shock_p  = "DD × Shock (+)",
    lev_std_shock_n = "Leverage × Shock (–)",
    dd_std_shock_n  = "DD × Shock (–)"
  ),
  gof_map  = c("r2" = "R²"),
  title    = "Respuestas heterogéneas por tipo de shock"
)


# Una sola regresion para shocks positivos y negativos

# Modelo único con ambos shocks
m_both <- feols(
  dlog_capital ~ 
    lev_std_shock_p + dd_std_shock_p +    # interacción con shocks positivos
    lev_std_shock_n + dd_std_shock_n +    # interacción con shocks negativos
    rsales_g_std + size_index + sh_current_a_std  # controles de firma
  | name + Country,                       # efectos fijos
  data    = df,
  cluster = ~ name + Country              # clustering
)

# Mostrar resultados
modelsummary(
  m_both,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    lev_std_shock_p = "Leverage × Shock (+)",
    dd_std_shock_p  = "DD × Shock (+)",
    lev_std_shock_n = "Leverage × Shock (–)",
    dd_std_shock_n  = "DD × Shock (–)"
  )
)



# Ver si es necesario o no. (Me da resultados un poco extraños).

# ----------------------------------------
# Table 15 (revisado v2): Arellano–Bond con rezagos macro y shock incluido
# ----------------------------------------

# 1) Definir controles
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_agg2  <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp"),
  paste0("L", 1:4, "_embigl")
)

# 2) Asegurar que existan las interacciones
if (!"lev_std_shock" %in% names(df)) {
  df <- df %>% mutate(lev_std_shock = lev_std * shock)
}
if (!"dd_std_shock" %in% names(df)) {
  df <- df %>% mutate(dd_std_shock  = dd_std  * shock)
}

# 3) Limpiar NAs e infinitos en todas las variables necesarias
vars_needed <- c(
  "dlog_capital", "lev_std_shock", "dd_std_shock",
  "shock",                # <<-- incluir shock
  controls_firm,
  controls_agg2,
  "name", "dateq"
)
df_gmm <- df %>%
  select(all_of(vars_needed)) %>%
  filter_if(is.numeric, all_vars(is.finite(.))) %>%
  drop_na()

# 4) Panel data
df_p <- pdata.frame(df_gmm, index = c("name", "dateq"))

# 5) Instrumentos: lag 2 de dlog_capital
inst_lag <- 2

# 6) Construir fórmulas
base_fml <- function(rhs_vars, inst) {
  paste0(
    "dlog_capital ~ lag(dlog_capital, 1) + ",
    paste(rhs_vars, collapse = " + "),
    " | lag(dlog_capital, ", inst, ")"
  ) %>% as.formula()
}

# 6.1) m1: solo lev_std_shock
f1 <- base_fml(c("lev_std_shock", controls_firm), inst_lag)

# 6.2) m2: solo dd_std_shock
f2 <- base_fml(c("dd_std_shock", controls_firm), inst_lag)

# 6.3) m3: ambas interacciones
f3 <- base_fml(c("lev_std_shock", "dd_std_shock", controls_firm), inst_lag)

# 6.4) m4: + shock puro + controles de ciclo (2 rezagos)
rhs4 <- c("shock", "lev_std_shock", "dd_std_shock", controls_firm, controls_agg2)
f4    <- base_fml(rhs4, inst_lag)

# 7) Estimar con pgmm (difference GMM, two‐step)
m1_gmm <- pgmm(f1, data = df_p,
               effect = "individual", model = "twosteps",
               transformation = "d", collapse = FALSE)

m2_gmm <- pgmm(f2, data = df_p,
               effect = "individual", model = "twosteps",
               transformation = "d", collapse = FALSE)

m3_gmm <- pgmm(f3, data = df_p,
               effect = "individual", model = "twosteps",
               transformation = "d", collapse = FALSE)

m4_gmm <- pgmm(f4, data = df_p,
               effect = "individual", model = "twosteps",
               transformation = "d", collapse = FALSE)

# 8) Mostrar resultados
modelsummary(
  list(
    "m1 (Lev×Shock)"             = m1_gmm,
    "m2 (DD×Shock)"              = m2_gmm,
    "m3 (Ambas interac.)"        = m3_gmm,
    "m4 (+Shock + 2 rezagos macro)" = m4_gmm
  ),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lag(dlog_capital, 1)" = "Lagged Δ log k",
    "lev_std_shock"        = "Leverage × Shock",
    "dd_std_shock"         = "DD × Shock",
    "shock"                = "Shock (puro)"
  ),
  gof_map  = c("r.squared" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  title    = "Table 15: Dynamic Panel (Arellano–Bond, 2 rezagos macro)"
)




# ==============================================
# Table 18 (modificado): Sensibilidades cíclicas con demeaning de variables cíclicas,
# sin EMBIGL y en dos conjuntos de modelos (m1–m3 sin macro, m4–m6 con macro l1–l2)
# ==============================================

# 1) Demean y estandarizar variables cíclicas por Country
df <- df %>%
  group_by(Country) %>%
  mutate(
    dlog_gdp_z = (dlog_gdp - mean(dlog_gdp, na.rm = TRUE)) / sd(dlog_gdp, na.rm = TRUE),
    dlog_cpi_z = (dlog_cpi - mean(dlog_cpi, na.rm = TRUE)) / sd(dlog_cpi, na.rm = TRUE),
    unemp_z    = (unemp    - mean(unemp,    na.rm = TRUE)) / sd(unemp,    na.rm = TRUE)
  ) %>%
  ungroup()

# 2) Crear interacciones “demeaning & std” con shock y con las variables cíclicas estandarizadas
df <- df %>%
  mutate(
    lev_std_shock = lev_std * shock,
    dd_std_shock  = dd_std  * shock,
    lev_std_gdp   = lev_std * dlog_gdp_z,
    lev_std_cpi   = lev_std * dlog_cpi_z,
    lev_std_ur    = lev_std * unemp_z,
    dd_std_gdp    = dd_std  * dlog_gdp_z,
    dd_std_cpi    = dd_std  * dlog_cpi_z,
    dd_std_ur     = dd_std  * unemp_z
  )

# 3) Definir controles
controls_firm      <- c("rsales_g_std", "size_index", "sh_current_a_std")
# Macro–controles reducidos: solo lags 1 y 2 de cada variable cíclica
controls_agg_small <- c(
  "L1_dlog_gdp", "L2_dlog_gdp",
  "L1_dlog_cpi", "L2_dlog_cpi",
  "L1_unemp",    "L2_unemp"
)

# 4) Conjunto A: sin controles macro
vars18a <- list(
  m1 = c("dlog_capital",
         "lev_std_shock", "lev_std_gdp", "lev_std_cpi", "lev_std_ur",
         controls_firm),
  m2 = c("dlog_capital",
         "dd_std_shock",  "dd_std_gdp",  "dd_std_cpi",  "dd_std_ur",
         controls_firm),
  m3 = c("dlog_capital",
         "lev_std_shock", "dd_std_shock",
         "lev_std_gdp",   "dd_std_gdp",
         "lev_std_cpi",   "dd_std_cpi",
         "lev_std_ur",    "dd_std_ur",
         controls_firm)
)

# 5) Conjunto B: con controles macro l1–l2
vars18b <- list(
  m4 = c("dlog_capital",
         "lev_std_shock", "lev_std_gdp",
         controls_firm, controls_agg_small),
  m5 = c("dlog_capital",
         "dd_std_shock",  "dd_std_gdp",
         controls_firm, controls_agg_small),
  m6 = c("dlog_capital",
         "lev_std_shock", "dd_std_shock",
         controls_firm, controls_agg_small)
)

# 6) Función auxiliar para estimar una lista de modelos
run_models <- function(df, vars_list) {
  models <- lapply(names(vars_list), function(m) {
    vars_i <- vars_list[[m]]
    df_i   <- df %>%
      select(name, Country, all_of(vars_i)) %>%
      drop_na()
    fmla <- as.formula(
      paste("dlog_capital ~",
            paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
            "| name + Country")
    )
    feols(fmla, data = df_i, cluster = ~name + Country)
  })
  names(models) <- names(vars_list)
  models
}

models18a <- run_models(df, vars18a)
models18b <- run_models(df, vars18b)

# 7) Presentar ambos conjuntos de modelos
modelsummary(
  c(models18a, models18b),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_shock" = "Leverage × Shock",
    "dd_std_shock"  = "DD × Shock",
    "lev_std_gdp"   = "Leverage × ΔGDP (z)",
    "lev_std_cpi"   = "Leverage × Δcpi (z)",
    "lev_std_ur"    = "Leverage × Unemp (z)",
    "dd_std_gdp"    = "DD × ΔGDP (z)",
    "dd_std_cpi"    = "DD × Δcpi (z)",
    "dd_std_ur"     = "DD × Unemp (z)"
  ),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = NULL,
  title    = "Table 18 Modificado: Sensibilidades cíclicas con demeaning y controles reducidos"
)





# ================================================
# Table 19 (alineado con Table 20): dos modelos por separado sin winsorización del choque
# ================================================

# 1) Crear versiones positiva y negativa del choque original
df <- df %>%
  mutate(
    shock_pos = pmax(shock,  0),
    shock_neg = pmin(shock,  0)
  )

# 2) Crear interacciones “demeaning & std” con apalancamiento y DD
df <- df %>%
  mutate(
    lev_std_shock_pos = lev_std * shock_pos,
    dd_std_shock_pos  = dd_std  * shock_pos,
    lev_std_shock_neg = lev_std * shock_neg,
    dd_std_shock_neg  = dd_std  * shock_neg
  )

# 3) Definir controles a nivel firma
controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")

# 4) Estimar dos regresiones por separado

# 4.1) Shock expansivo
m19_pos <- feols(
  dlog_capital ~ 
    lev_std_shock_pos + dd_std_shock_pos +
    rsales_g_std + size_index + sh_current_a_std
  | name + Country,
  data    = df,
  cluster = ~name + Country
)

# 4.2) Shock contractivo
m19_neg <- feols(
  dlog_capital ~ 
    lev_std_shock_neg + dd_std_shock_neg +
    rsales_g_std + size_index + sh_current_a_std
  | name + Country,
  data    = df,
  cluster = ~name + Country
)

# 5) Presentar resultados lado a lado
modelsummary(
  list(
    "Shock Expansivo"   = m19_pos,
    "Shock Contractivo" = m19_neg
  ),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    lev_std_shock_pos = "Leverage × Shock (+)",
    dd_std_shock_pos  = "DD × Shock (+)",
    lev_std_shock_neg = "Leverage × Shock (–)",
    dd_std_shock_neg  = "DD × Shock (–)"
  ),
  gof_map  = c("r2" = "R²"),
  title    = "Table 19: Respuestas heterogéneas a shocks"
)





# ================================================
# Table 20: Alternative Time Aggregation
#            (solo “shock” + FE firma + FE país + EMBIGL)
# ================================================

# 1) Crear rezagos de shock y shock_sum si no existen
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

# 2) Crear interacciones de suma si faltan
df <- df %>%
  mutate(
    lev_std_sum = if (!"lev_std_sum" %in% names(df)) { message("Aviso: 'lev_std_sum' creado aquí."); lev_std * shock_sum } else lev_std_sum,
    dd_std_sum  = if (!"dd_std_sum"  %in% names(df)) { message("Aviso: 'dd_std_sum' creado aquí.");  dd_std  * shock_sum } else dd_std_sum
  )

# 3) Asegurar controles a nivel firma y agregados
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
  message("Aviso: 'controls_firm' creada aquí.")
}
if (!exists("controls_agg")) {
  controls_agg <- if (exists("controls_aggregate")) controls_aggregate else c(
    paste0("L", 1:4, "_dlog_gdp"),
    paste0("L", 1:4, "_dlog_cpi"),
    paste0("L", 1:4, "_unemp"),
    paste0("L", 1:4, "_embigl")
  )
  message("Aviso: 'controls_agg' creada aquí.")
}

# 4) Asegurar interacciones trimestrales de ciclo y choque si faltan
df <- df %>%
  mutate(
    lev_std_gdp   = if (!"lev_std_gdp"   %in% names(df)) { message("Aviso: 'lev_std_gdp' creado aquí.");   lev_std * dlog_gdp } else lev_std_gdp,
    dd_std_gdp    = if (!"dd_std_gdp"    %in% names(df)) { message("Aviso: 'dd_std_gdp' creado aquí.");    dd_std  * dlog_gdp } else dd_std_gdp,
    lev_std_shock = if (!"lev_std_shock" %in% names(df)) { message("Aviso: 'lev_std_shock' creado aquí."); lev_std * shock  } else lev_std_shock,
    dd_std_shock  = if (!"dd_std_shock"  %in% names(df)) { message("Aviso: 'dd_std_shock' creado aquí.");  dd_std  * shock  } else dd_std_shock
  )

# 5) Estimar modelos m20_1–m20_4  
#    Construyo fórmulas y datos filtrados para cada caso:

# M1: leverage only  
f20_1 <- as.formula(paste0(
  "dlog_capital ~ lev_std_gdp + lev_std_sum + lev_std_shock + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
df20_1 <- df %>% drop_na(all.vars(f20_1)[-1])
m20_1 <- feols(f20_1, data = df20_1, cluster = ~name + Country)

# M2: DD only  
f20_2 <- as.formula(paste0(
  "dlog_capital ~ dd_std_gdp + dd_std_sum + dd_std_shock + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
df20_2 <- df %>% drop_na(all.vars(f20_2)[-1])
m20_2 <- feols(f20_2, data = df20_2, cluster = ~name + Country)

# M3: Leverage + DD  
f20_3 <- as.formula(paste0(
  "dlog_capital ~ lev_std_gdp + dd_std_gdp + lev_std_sum + lev_std_shock + dd_std_sum + dd_std_shock + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
df20_3 <- df %>% drop_na(all.vars(f20_3)[-1])
m20_3 <- feols(f20_3, data = df20_3, cluster = ~name + Country)

# M4: shock_sum + controles agregados  
rhs20_4 <- c(
  "lev_std_gdp", "shock_sum", "lev_std_sum", "lev_std_shock",
  "dd_std_sum", "dd_std_shock", controls_firm, "embigl", controls_agg
)
f20_4 <- as.formula(
  paste("dlog_capital ~", paste(rhs20_4, collapse = " + "), "| name + Country")
)
df20_4 <- df %>% drop_na(all_of(rhs20_4))
m20_4 <- feols(f20_4, data = df20_4, cluster = ~name + Country)

# 6) Mostrar Table 20  
modelsummary(
  list(M1 = m20_1, M2 = m20_2, M3 = m20_3, M4 = m20_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_sum"   = "Leverage × shock (sum)",
    "dd_std_sum"    = "DD × shock (sum)",
    "shock_sum"     = "Sum of shocks",
    "lev_std_gdp"   = "Leverage × ΔGDP",
    "dd_std_gdp"    = "DD × ΔGDP",
    "lev_std_shock" = "Leverage × shock",
    "dd_std_shock"  = "DD × shock"
  ),
  keep     = c("lev_std_sum", "dd_std_sum", "shock_sum"),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 20: Alternative Time Aggregation"
)






## Version 20, similar al 11:

# ================================================
# Table 20: Alternative Time Aggregation
#            (solo “shock_sum” + FE firma + FE país + EMBIGL)
#            Adaptado a partir del código de choques asimétricos
# ================================================

# 1) Crear rezagos de shock y sumar
if (!"L0_shock" %in% names(df)) {
  df <- df %>%
    arrange(name, dateq) %>%
    group_by(name) %>%
    mutate(
      L0_shock   = shock,
      L1_shock   = lag(shock, 1),
      L2_shock   = lag(shock, 2),
      L3_shock   = lag(shock, 3),
      L4_shock   = lag(shock, 4),
      shock_sum  = L0_shock + L1_shock + L2_shock + L3_shock + L4_shock
    ) %>%
    ungroup()
  message("Aviso: rezagos L0–L4 de shock y shock_sum creados aquí.")
}

# 2) Normalizar shock_sum (z-score global) para hacerlo comparable
if (!"shock_sum_z" %in% names(df)) {
  df <- df %>%
    mutate(
      shock_sum_z = (shock_sum - mean(shock_sum, na.rm = TRUE)) / sd(shock_sum, na.rm = TRUE)
    )
  message("Aviso: 'shock_sum_z' creado aquí.")
}

# 3) Crear interacciones “demeaning & std” con shock_sum_z
df <- df %>%
  mutate(
    lev_std_shock_sum = lev_std * shock_sum_z,
    dd_std_shock_sum  = dd_std  * shock_sum_z
  )

# 4) Asegurar controles a nivel firma y agregados (si no existen)
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

# 5) Asegurar interacción de ciclo (apalancamiento × ΔGDP) si faltan
df <- df %>%
  mutate(
    lev_std_gdp = if (!"lev_std_gdp"   %in% names(df)) { 
      message("Aviso: 'lev_std_gdp' creado aquí."); 
      lev_std * dlog_gdp 
    } else lev_std_gdp,
    dd_std_gdp  = if (!"dd_std_gdp"    %in% names(df)) { 
      message("Aviso: 'dd_std_gdp' creado aquí.");  
      dd_std  * dlog_gdp 
    } else dd_std_gdp
  )

# 6) Estimar modelos m20_1–m20_4

# M1: canal apalancamiento + ciclo y suma
f20_1 <- as.formula(paste0(
  "dlog_capital ~ lev_std_gdp + lev_std_shock_sum + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
df20_1 <- df %>% drop_na(all.vars(f20_1)[-1])
m20_1 <- feols(f20_1, data = df20_1, cluster = ~name + Country)

# M2: canal solvencia + ciclo y suma
f20_2 <- as.formula(paste0(
  "dlog_capital ~ dd_std_gdp + dd_std_shock_sum + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
df20_2 <- df %>% drop_na(all.vars(f20_2)[-1])
m20_2 <- feols(f20_2, data = df20_2, cluster = ~name + Country)

# M3: ambos canales + ciclo y suma
f20_3 <- as.formula(paste0(
  "dlog_capital ~ lev_std_gdp + dd_std_gdp + ",
  "lev_std_shock_sum + dd_std_shock_sum + ",
  paste(controls_firm, collapse = " + "), " + embigl | name + Country"
))
df20_3 <- df %>% drop_na(all.vars(f20_3)[-1])
m20_3 <- feols(f20_3, data = df20_3, cluster = ~name + Country)

# M4: shock_sum puro + canales + controles agregados
rhs20_4 <- c(
  "lev_std_gdp", "shock_sum_z", "lev_std_shock_sum", "dd_std_shock_sum",
  controls_firm, "embigl", controls_agg
)
f20_4 <- as.formula(
  paste("dlog_capital ~", paste(rhs20_4, collapse = " + "), "| name + Country")
)
df20_4 <- df %>% drop_na(all_of(rhs20_4))
m20_4 <- feols(f20_4, data = df20_4, cluster = ~name + Country)

# 7) Mostrar Table 20
modelsummary(
  list(M1 = m20_1, M2 = m20_2, M3 = m20_3, M4 = m20_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    lev_std_shock_sum = "Leverage × shock (sum)",
    dd_std_shock_sum  = "DD × shock (sum)",
    shock_sum_z       = "Sum of shocks (z-score)",
    lev_std_gdp       = "Leverage × ΔGDP",
    dd_std_gdp        = "DD × ΔGDP"
  ),
  keep     = c("lev_std_shock_sum", "dd_std_shock_sum", "shock_sum_z"),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 20: Alternative Time Aggregation"
)







# ================================================
# Table 21: Interaction with Other Firm‐Level Covariates (con winsorización)
# ================================================

# 0) Definir función winsorize
winsorize <- function(x, p = 0.01) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

# 1) Calcular crecimiento futuro de ventas (4 trimestres adelante) y winsorizarlo
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    rsales_Fg4_raw = log(dplyr::lead(saleq, 4)) - log(dplyr::lead(saleq, 3)),
    rsales_Fg4_w   = winsorize(rsales_Fg4_raw)
  ) %>%
  ungroup() %>%
  mutate(
    rsales_Fg4_std_w = (rsales_Fg4_w - mean(rsales_Fg4_w, na.rm = TRUE)) /
      sd(rsales_Fg4_w,    na.rm = TRUE)
  )

# 2) Winsorizar y re‐estandarizar las demás covariables
df <- df %>%
  mutate(
    shock_w            = winsorize(shock),
    lev_std_w          = winsorize(lev_std),
    dd_std_w           = winsorize(dd_std),
    rsales_g_std_w     = winsorize(rsales_g_std),
    size_index_w       = winsorize(size_index),
    sh_current_a_std_w = winsorize(sh_current_a_std)
  ) %>%
  mutate(
    lev_std_win         = (lev_std_w          - mean(lev_std_w,          na.rm=TRUE)) /
      sd(lev_std_w,          na.rm=TRUE),
    dd_std_win          = (dd_std_w           - mean(dd_std_w,           na.rm=TRUE)) /
      sd(dd_std_w,           na.rm=TRUE),
    rsales_g_std_win    = (rsales_g_std_w     - mean(rsales_g_std_w,     na.rm=TRUE)) /
      sd(rsales_g_std_w,     na.rm=TRUE),
    size_index_win      = (size_index_w       - mean(size_index_w,       na.rm=TRUE)) /
      sd(size_index_w,       na.rm=TRUE),
    sh_current_a_std_win= (sh_current_a_std_w - mean(sh_current_a_std_w, na.rm=TRUE)) /
      sd(sh_current_a_std_w, na.rm=TRUE)
  )

# 3) Crear interacciones con el choque winsorizado
df <- df %>%
  mutate(
    lev_std_shock          = lev_std_win          * shock_w,
    dd_std_shock           = dd_std_win           * shock_w,
    rsales_g_std_shock     = rsales_g_std_win     * shock_w,
    rsales_Fg4_std_shock   = rsales_Fg4_std_w     * shock_w,
    size_index_shock       = size_index_win       * shock_w,
    sh_current_a_std_shock = sh_current_a_std_win * shock_w
  )

# 4) Definir controles
controls_firm <- c("rsales_g_std_win", "size_index_win", "sh_current_a_std_win")
controls_agg  <- c(
  paste0("L",1:4,"_dlog_gdp"),
  paste0("L",1:4,"_dlog_cpi"),
  paste0("L",1:4,"_unemp"),
  paste0("L",1:4,"_embigl")
)

# 5) Especificaciones m1–m8
vars21 <- list(
  m1 = c("dlog_capital",
         "lev_std_shock", "rsales_g_std_shock",
         controls_firm, "embigl"),
  m2 = c("dlog_capital",
         "dd_std_shock",  "rsales_g_std_shock",
         controls_firm, "embigl"),
  m3 = c("dlog_capital",
         "lev_std_shock", "rsales_Fg4_std_shock",
         "rsales_Fg4_std_w", controls_firm, "embigl"),
  m4 = c("dlog_capital",
         "dd_std_shock",  "rsales_Fg4_std_shock",
         "rsales_Fg4_std_w", controls_firm, "embigl"),
  m5 = c("dlog_capital",
         "lev_std_shock", "size_index_shock",
         controls_firm, "embigl"),
  m6 = c("dlog_capital",
         "dd_std_shock",  "size_index_shock",
         controls_firm, "embigl"),
  m7 = c("dlog_capital",
         "lev_std_shock", "sh_current_a_std_shock",
         "sh_current_a_std_w", controls_firm, "embigl"),
  m8 = c("dlog_capital",
         "dd_std_shock",  "sh_current_a_std_shock",
         "sh_current_a_std_w", controls_firm, "embigl")
)

# 6) Estimar m1–m8 con feols()
models21 <- lapply(vars21, function(vars) {
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

# 7) Mostrar Table 21
modelsummary(
  models21,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_std_shock", "dd_std_shock",
               "rsales_g_std_shock", "rsales_Fg4_std_shock",
               "size_index_shock", "sh_current_a_std_shock"),
  order    = c("lev_std_shock", "dd_std_shock",
               "rsales_g_std_shock", "rsales_Fg4_std_shock",
               "size_index_shock", "sh_current_a_std_shock"),
  coef_map = c(
    "lev_std_shock"           = "Leverage × shock",
    "dd_std_shock"            = "DD × shock",
    "rsales_g_std_shock"      = "Sales growth × shock",
    "rsales_Fg4_std_shock"    = "Future sales growth × shock",
    "size_index_shock"        = "Size × shock",
    "sh_current_a_std_shock"  = "Liquidity × shock"
  ),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  title    = "Table 21: Interaction with Other Firm‐Level Covariates (winsorizadas)"
)



# ===============================================================
# Table 22: Interacción con Otras Medidas de Posición Financiera
#            (solo “shock” + FE firma + FE país + EMBIGL)
# ===============================================================

# 1) Estandarizar y crear nuevos proxies solo si faltan
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    lev_std = if (!"lev_std" %in% names(df)) {
      message("Aviso: 'lev_std' creada aquí.")
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE)
    } else lev_std,
    
    dd_std = if (!"dd_std" %in% names(df)) {
      message("Aviso: 'dd_std' creada aquí.")
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE)
    } else dd_std,
    
    cashop_std = if (!"cashop_std" %in% names(df)) {
      message("Aviso: 'cashop_std' creada aquí.")
      (cashop - mean(cashop, na.rm = TRUE)) / sd(cashop, na.rm = TRUE)
    } else cashop_std,
    
    div_std = if (!"div_std" %in% names(df)) {
      message("Aviso: 'div_std' creada aquí.")
      (dividends - mean(dividends, na.rm = TRUE)) / sd(dividends, na.rm = TRUE)
    } else div_std,
    
    liq_std = if (!"liq_std" %in% names(df)) {
      message("Aviso: 'liq_std' creada aquí.")
      sh_current_a_std
    } else liq_std
  ) %>%
  ungroup()

# 2) Crear interacciones con shock y ΔGDP solo si faltan
df <- df %>%
  mutate(
    lev_std_shock  = if (!"lev_std_shock" %in% names(df))  { message("Aviso: 'lev_std_shock' creada aquí.");  lev_std * shock }        else lev_std_shock,
    dd_std_shock   = if (!"dd_std_shock"  %in% names(df))  { message("Aviso: 'dd_std_shock' creada aquí.");   dd_std  * shock }        else dd_std_shock,
    size_idx_shock = if (!"size_idx_shock" %in% names(df)){ message("Aviso: 'size_idx_shock' creada aquí."); size_index * shock }       else size_idx_shock,
    cashop_std_shock = if (!"cashop_std_shock" %in% names(df)) { message("Aviso: 'cashop_std_shock' creada aquí."); cashop_std * shock } else cashop_std_shock,
    div_std_shock    = if (!"div_std_shock"    %in% names(df)) { message("Aviso: 'div_std_shock' creada aquí.");    div_std * shock }    else div_std_shock,
    liq_std_shock    = if (!"liq_std_shock"    %in% names(df)) { message("Aviso: 'liq_std_shock' creada aquí.");    liq_std * shock }    else liq_std_shock
  )

# 3) Asegurar controles de firma y agregados
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
  message("Aviso: 'controls_firm' creada aquí.")
}
if (!exists("controls_agg")) {
  controls_agg <- if (exists("controls_aggregate")) controls_aggregate else c(
    paste0("L",1:4,"_dlog_gdp"),
    paste0("L",1:4,"_dlog_cpi"),
    paste0("L",1:4,"_unemp"),
    paste0("L",1:4,"_embigl")
  )
  message("Aviso: 'controls_agg' creada aquí.")
}

# 4) Definir variables para m1–m8
vars22 <- list(
  m1 = c("dlog_capital","lev_std_gdp","lev_std_shock","size_idx_shock","size_index",controls_firm,"embigl"),
  m2 = c("dlog_capital","dd_std_gdp","dd_std_shock","size_idx_shock","size_index",controls_firm,"embigl"),
  m3 = c("dlog_capital","lev_std_gdp","lev_std_shock","cashop_std_shock","cashop_std",controls_firm,"embigl"),
  m4 = c("dlog_capital","dd_std_gdp","dd_std_shock","cashop_std_shock","cashop_std",controls_firm,"embigl"),
  m5 = c("dlog_capital","lev_std_gdp","lev_std_shock","div_std_shock","div_std",controls_firm,"embigl"),
  m6 = c("dlog_capital","dd_std_gdp","dd_std_shock","div_std_shock","div_std",controls_firm,"embigl"),
  m7 = c("dlog_capital","lev_std_gdp","lev_std_shock","liq_std_shock","liq_std",controls_firm,"embigl"),
  m8 = c("dlog_capital","dd_std_gdp","dd_std_shock","liq_std_shock","liq_std",controls_firm,"embigl")
)

# 5) Estimar modelos m1–m8
models22 <- lapply(vars22, function(v) {
  df_sub <- df %>% select(name, Country, all_of(v)) %>% drop_na()
  f      <- as.formula(paste("dlog_capital ~",
                             paste(setdiff(v,"dlog_capital"), collapse=" + "),
                             "| name + Country"))
  feols(f, data = df_sub, cluster = ~name + Country)
})

# 6) Mostrar Table 22 (solo interacciones, R², AIC, BIC)
modelsummary(
  models22,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_std_shock","dd_std_shock",
               "size_idx_shock","cashop_std_shock",
               "div_std_shock","liq_std_shock"),
  order    = c("lev_std_shock","dd_std_shock",
               "size_idx_shock","cashop_std_shock",
               "div_std_shock","liq_std_shock"),
  coef_map = c(
    "lev_std_shock"      = "Leverage × shock",
    "dd_std_shock"       = "DD × shock",
    "size_idx_shock"     = "Size_index × shock",
    "cashop_std_shock"   = "Cash flow × shock",
    "div_std_shock"      = "Dividend × shock",
    "liq_std_shock"      = "Liquidity × shock"
  ),
  gof_map  = c("r2" = "R²","AIC" = "AIC","BIC" = "BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 22: Interaction with Other Measures of Financial Positions"
)





# -------------------------------------------------------------------------
# Table 23: IV usando posiciones financieras rezagadas
#            (solo “shock” + FE firma + país + EMBIGL)
#            con rezagos 1, 2 y 3 trimestres y usando size_index
# -------------------------------------------------------------------------

# 0) Crear interacciones endógenas y sus rezagos solo si aún no existen
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    # Interacciones trimestrales
    lev_std_gdp   = if (!"lev_std_gdp"   %in% names(df)) { message("Aviso: 'lev_std_gdp' creada aquí.");   lev_std * dlog_gdp   } else lev_std_gdp,
    lev_std_shock = if (!"lev_std_shock" %in% names(df)) { message("Aviso: 'lev_std_shock' creada aquí."); lev_std * shock     } else lev_std_shock,
    dd_std_gdp    = if (!"dd_std_gdp"    %in% names(df)) { message("Aviso: 'dd_std_gdp' creada aquí.");    dd_std  * dlog_gdp } else dd_std_gdp,
    dd_std_shock  = if (!"dd_std_shock"  %in% names(df)) { message("Aviso: 'dd_std_shock' creada aquí.");  dd_std  * shock     } else dd_std_shock
  ) %>%
  mutate(
    # Rezagos de leverage
    levL1_std_gdp   = if (!"levL1_std_gdp"   %in% names(df)) { message("Aviso: 'levL1_std_gdp' creada aquí.");   lag(lev_std_gdp,   1) } else levL1_std_gdp,
    levL1_std_shock = if (!"levL1_std_shock" %in% names(df)) { message("Aviso: 'levL1_std_shock' creada aquí."); lag(lev_std_shock, 1) } else levL1_std_shock,
    levL2_std_gdp   = if (!"levL2_std_gdp"   %in% names(df)) { message("Aviso: 'levL2_std_gdp' creada aquí.");   lag(lev_std_gdp,   2) } else levL2_std_gdp,
    levL2_std_shock = if (!"levL2_std_shock" %in% names(df)) { message("Aviso: 'levL2_std_shock' creada aquí."); lag(lev_std_shock, 2) } else levL2_std_shock,
    levL3_std_gdp   = if (!"levL3_std_gdp"   %in% names(df)) { message("Aviso: 'levL3_std_gdp' creada aquí.");   lag(lev_std_gdp,   3) } else levL3_std_gdp,
    levL3_std_shock = if (!"levL3_std_shock" %in% names(df)) { message("Aviso: 'levL3_std_shock' creada aquí."); lag(lev_std_shock, 3) } else levL3_std_shock
  ) %>%
  mutate(
    # Rezagos de DD
    ddL1_std_gdp   = if (!"ddL1_std_gdp"   %in% names(df)) { message("Aviso: 'ddL1_std_gdp' creada aquí.");   lag(dd_std_gdp,   1) } else ddL1_std_gdp,
    ddL1_std_shock = if (!"ddL1_std_shock" %in% names(df)) { message("Aviso: 'ddL1_std_shock' creada aquí."); lag(dd_std_shock, 1) } else ddL1_std_shock,
    ddL2_std_gdp   = if (!"ddL2_std_gdp"   %in% names(df)) { message("Aviso: 'ddL2_std_gdp' creada aquí.");   lag(dd_std_gdp,   2) } else ddL2_std_gdp,
    ddL2_std_shock = if (!"ddL2_std_shock" %in% names(df)) { message("Aviso: 'ddL2_std_shock' creada aquí."); lag(dd_std_shock, 2) } else ddL2_std_shock,
    ddL3_std_gdp   = if (!"ddL3_std_gdp"   %in% names(df)) { message("Aviso: 'ddL3_std_gdp' creada aquí.");   lag(dd_std_gdp,   3) } else ddL3_std_gdp,
    ddL3_std_shock = if (!"ddL3_std_shock" %in% names(df)) { message("Aviso: 'ddL3_std_shock' creada aquí."); lag(dd_std_shock, 3) } else ddL3_std_shock
  ) %>%
  ungroup()

# 1) Asegurar controles a nivel firma
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
  message("Aviso: 'controls_firm' creada aquí.")
}

# 2) Definir y estimar Modelos IV

## Modelo 1: Leverage ~ lag1
df23_1 <- df %>% drop_na(dlog_capital, rsales_g_std, size_index, sh_current_a_std, embigl,
                         lev_std_gdp, lev_std_shock,
                         levL1_std_gdp, levL1_std_shock)
m23_1 <- feols(
  dlog_capital ~ rsales_g_std + size_index + sh_current_a_std + embigl |
    name + Country |
    (lev_std_gdp + lev_std_shock) ~ (levL1_std_gdp + levL1_std_shock),
  data    = df23_1,
  cluster = ~name + Country
)

## Modelo 2: Leverage ~ lag2
df23_2 <- df %>% drop_na(dlog_capital, rsales_g_std, size_index, sh_current_a_std, embigl,
                         lev_std_gdp, lev_std_shock,
                         levL2_std_gdp, levL2_std_shock)
m23_2 <- feols(
  dlog_capital ~ rsales_g_std + size_index + sh_current_a_std + embigl |
    name + Country |
    (lev_std_gdp + lev_std_shock) ~ (levL2_std_gdp + levL2_std_shock),
  data    = df23_2,
  cluster = ~name + Country
)

## Modelo 3: Leverage ~ lag3
df23_3 <- df %>% drop_na(dlog_capital, rsales_g_std, size_index, sh_current_a_std, embigl,
                         lev_std_gdp, lev_std_shock,
                         levL3_std_gdp, levL3_std_shock)
m23_3 <- feols(
  dlog_capital ~ rsales_g_std + size_index + sh_current_a_std + embigl |
    name + Country |
    (lev_std_gdp + lev_std_shock) ~ (levL3_std_gdp + levL3_std_shock),
  data    = df23_3,
  cluster = ~name + Country
)

## Modelo 4–6: DD instrumentado igual para lag1, lag2, lag3
for (i in 1:3) {
  lag_gdp <- paste0("ddL", i, "_std_gdp")
  lag_shk <- paste0("ddL", i, "_std_shock")
  df23   <- df %>% drop_na(dlog_capital, rsales_g_std, size_index, sh_current_a_std, embigl,
                           dd_std_gdp, dd_std_shock,
                           !!sym(lag_gdp), !!sym(lag_shk))
  assign(paste0("m23_", i+3),
         feols(
           as.formula(
             paste0(
               "dlog_capital ~ rsales_g_std + size_index + sh_current_a_std + embigl | name + Country | ",
               "(dd_std_gdp + dd_std_shock) ~ (", lag_gdp, " + ", lag_shk, ")"
             )
           ),
           data    = df23,
           cluster = ~name + Country
         )
  )
}

# 3) Mostrar Table 23 con R², AIC y BIC (solo coeficientes de shock)
modelsummary(
  list(
    "Lev 1q" = m23_1,
    "Lev 2q" = m23_2,
    "Lev 3q" = m23_3,
    "DD 1q"  = m23_4,
    "DD 2q"  = m23_5,
    "DD 3q"  = m23_6
  ),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_std_shock", "dd_std_shock"),
  order    = c("lev_std_shock", "dd_std_shock"),
  gof_map  = c("r2" = "R²", "AIC" = "AIC", "BIC" = "BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 23: IV usando posiciones financieras rezagadas"
)



# --------------------------------------------------------
# Table 24: Decomposition of Leverage
#            (solo “shock” + FE firma + país + EMBIGL)
# --------------------------------------------------------

# 1) Crear y estandarizar medidas de posición financiera solo si faltan
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    lev_std = if (!"lev_std" %in% names(df)) {
      message("Aviso: 'lev_std' creada aquí.")
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE)
    } else lev_std,
    
    levnet     = if (!"levnet" %in% names(df)) {
      message("Aviso: 'levnet' creada aquí.")
      (tdebt - cheq1) / atq
    } else levnet,
    levnet_std = if (!"levnet_std" %in% names(df)) {
      message("Aviso: 'levnet_std' creada aquí.")
      (levnet - mean(levnet, na.rm = TRUE)) / sd(levnet, na.rm = TRUE)
    } else levnet_std,
    
    shstdt     = if (!"shstdt" %in% names(df)) {
      message("Aviso: 'shstdt' creada aquí.")
      dlcq / tdebt
    } else shstdt,
    shstdt_std = if (!"shstdt_std" %in% names(df)) {
      message("Aviso: 'shstdt_std' creada aquí.")
      (shstdt - mean(shstdt, na.rm = TRUE)) / sd(shstdt, na.rm = TRUE)
    } else shstdt_std,
    
    shltdt     = if (!"shltdt" %in% names(df)) {
      message("Aviso: 'shltdt' creada aquí.")
      dlttq / tdebt
    } else shltdt,
    shltdt_std = if (!"shltdt_std" %in% names(df)) {
      message("Aviso: 'shltdt_std' creada aquí.")
      (shltdt - mean(shltdt, na.rm = TRUE)) / sd(shltdt, na.rm = TRUE)
    } else shltdt_std,
    
    sh_ol      = if (!"sh_ol" %in% names(df)) {
      message("Aviso: 'sh_ol' creada aquí.")
      (atq - tdebt) / atq
    } else sh_ol,
    sh_ol_std  = if (!"sh_ol_std" %in% names(df)) {
      message("Aviso: 'sh_ol_std' creada aquí.")
      (sh_ol - mean(sh_ol, na.rm = TRUE)) / sd(sh_ol, na.rm = TRUE)
    } else sh_ol_std,
    
    sh_l_std = if (!"sh_l_std" %in% names(df)) {
      message("Aviso: 'sh_l_std' creada aquí.")
      lev_std
    } else sh_l_std
  ) %>%
  ungroup()

# 2) Crear interacciones estándar con ΔGDP y shock solo si faltan
df <- df %>%
  mutate(
    lev_std_gdp      = if (!"lev_std_gdp" %in% names(df))      { message("Aviso: 'lev_std_gdp' creada aquí.");      lev_std    * dlog_gdp }    else lev_std_gdp,
    lev_std_shock    = if (!"lev_std_shock" %in% names(df))    { message("Aviso: 'lev_std_shock' creada aquí.");    lev_std    * shock   }    else lev_std_shock,
    
    levnet_std_gdp   = if (!"levnet_std_gdp" %in% names(df))   { message("Aviso: 'levnet_std_gdp' creada aquí.");   levnet_std * dlog_gdp }    else levnet_std_gdp,
    levnet_std_shock = if (!"levnet_std_shock" %in% names(df)) { message("Aviso: 'levnet_std_shock' creada aquí."); levnet_std * shock   }    else levnet_std_shock,
    
    shstdt_std_gdp   = if (!"shstdt_std_gdp" %in% names(df))   { message("Aviso: 'shstdt_std_gdp' creada aquí.");   shstdt_std * dlog_gdp }  else shstdt_std_gdp,
    shstdt_std_shock = if (!"shstdt_std_shock" %in% names(df)) { message("Aviso: 'shstdt_std_shock' creada aquí."); shstdt_std * shock   }  else shstdt_std_shock,
    
    shltdt_std_gdp   = if (!"shltdt_std_gdp" %in% names(df))   { message("Aviso: 'shltdt_std_gdp' creada aquí.");   shltdt_std * dlog_gdp }  else shltdt_std_gdp,
    shltdt_std_shock = if (!"shltdt_std_shock" %in% names(df)) { message("Aviso: 'shltdt_std_shock' creada aquí."); shltdt_std * shock   }  else shltdt_std_shock,
    
    sh_ol_std_gdp    = if (!"sh_ol_std_gdp" %in% names(df))    { message("Aviso: 'sh_ol_std_gdp' creada aquí.");    sh_ol_std  * dlog_gdp }    else sh_ol_std_gdp,
    sh_ol_std_shock  = if (!"sh_ol_std_shock" %in% names(df))  { message("Aviso: 'sh_ol_std_shock' creada aquí.");  sh_ol_std  * shock   }    else sh_ol_std_shock,
    
    sh_l_std_gdp     = if (!"sh_l_std_gdp" %in% names(df))     { message("Aviso: 'sh_l_std_gdp' creada aquí.");     sh_l_std   * dlog_gdp }    else sh_l_std_gdp,
    sh_l_std_shock   = if (!"sh_l_std_shock" %in% names(df))   { message("Aviso: 'sh_l_std_shock' creada aquí.");   sh_l_std   * shock   }    else sh_l_std_shock
  )

# 3) Asegurar controles de firma y agregados
if (!exists("controls_firm")) {
  controls_firm <- c("rsales_g_std", "size_index", "sh_current_a_std")
  message("Aviso: 'controls_firm' creada aquí.")
}
if (!exists("controls_agg")) {
  controls_agg <- if (exists("controls_aggregate")) controls_aggregate else c(
    paste0("L",1:4,"_dlog_gdp"),
    paste0("L",1:4,"_dlog_cpi"),
    paste0("L",1:4,"_unemp"),
    paste0("L",1:4,"_embigl")
  )
  message("Aviso: 'controls_agg' creada aquí.")
}

# 4) Listas de variables para m1–m7
vars24 <- list(
  m1 = c("dlog_capital","lev_std_gdp","lev_std_shock",controls_firm,"embigl"),
  m2 = c("dlog_capital","levnet_std_gdp","levnet_std_shock",controls_firm,"embigl"),
  m3 = c("dlog_capital","shstdt_std_gdp","shstdt_std_shock",controls_firm,"embigl"),
  m4 = c("dlog_capital","shltdt_std_gdp","shltdt_std_shock",controls_firm,"embigl"),
  m5 = c("dlog_capital","shstdt_std_gdp","shstdt_std_shock","shltdt_std_shock",controls_firm,"embigl"),
  m6 = c("dlog_capital","sh_ol_std_gdp","sh_ol_std_shock",controls_firm,"embigl"),
  m7 = c("dlog_capital","sh_l_std_gdp","sh_l_std_shock",controls_firm,"embigl")
)

# 5) Estimar modelos m1–m7
models24 <- lapply(vars24, function(vs) {
  df_sub <- df %>% select(name, Country, all_of(vs)) %>% drop_na()
  fml    <- as.formula(paste(
    "dlog_capital ~", paste(setdiff(vs,"dlog_capital"), collapse=" + "), "| name + Country"
  ))
  feols(fml, data = df_sub, cluster = ~name + Country)
})

# 6) Mostrar Table 24 (solo coeficientes de “shock”)
modelsummary(
  models24,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c(
    "lev_std_shock","levnet_std_shock","shstdt_std_shock",
    "shltdt_std_shock","sh_ol_std_shock","sh_l_std_shock"
  ),
  order    = c(
    "lev_std_shock","levnet_std_shock","shstdt_std_shock",
    "shltdt_std_shock","sh_ol_std_shock","sh_l_std_shock"
  ),
  coef_map = c(
    "lev_std_shock"       = "Leverage × shock",
    "levnet_std_shock"    = "Net leverage × shock",
    "shstdt_std_shock"    = "ST debt share × shock",
    "shltdt_std_shock"    = "LT debt share × shock",
    "sh_ol_std_shock"     = "Other liabilities share × shock",
    "sh_l_std_shock"      = "Total liabilities share × shock"
  ),
  gof_map  = c("r2" = "R²","AIC"="AIC","BIC"="BIC"),
  gof_omit = "^(?!r2$|AIC$|BIC$)[A-Za-z0-9_.]+",
  title    = "Table 24: Decomposition of Leverage"
)

