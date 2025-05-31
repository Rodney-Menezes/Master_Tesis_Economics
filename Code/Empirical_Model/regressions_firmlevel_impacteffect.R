# Instalar paquetes necesarios (si aún no lo están)
install.packages(c(
  "modelsummary",
  "fixest",
  "haven",
  "dplyr",
  "tidyr",
  "zoo"
))

# 0) Cargar librerías
library(haven)        # read_dta()
library(dplyr)        # manipulación de datos
library(zoo)          # as.yearqtr()
library(fixest)       # feols()
library(modelsummary) # presentación de tablas
library(tidyr)        # drop_na()

# 1) Leer datos y convertir fecha trimestral
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(dateq = as.yearqtr(dateq))

# 2) Construir variables macro y sus lags 1–4
df <- df %>%
  arrange(name, dateq) %>%
  mutate(
    dlog_gdp = log(gdp) - log(lag(gdp)),
    dlog_cpi = log(cpi) - log(lag(cpi)),
    unemp    = unemp
  ) %>%
  group_by(name) %>%
  mutate(
    L1_dlog_gdp = lag(dlog_gdp, 1),
    L2_dlog_gdp = lag(dlog_gdp, 2),
    L3_dlog_gdp = lag(dlog_gdp, 3),
    L4_dlog_gdp = lag(dlog_gdp, 4),
    L1_dlog_cpi = lag(dlog_cpi, 1),
    L2_dlog_cpi = lag(dlog_cpi, 2),
    L3_dlog_cpi = lag(dlog_cpi, 3),
    L4_dlog_cpi = lag(dlog_cpi, 4),
    L1_unemp    = lag(unemp,    1),
    L2_unemp    = lag(unemp,    2),
    L3_unemp    = lag(unemp,    3),
    L4_unemp    = lag(unemp,    4)
  ) %>%
  ungroup()

# 3) Controles a nivel firma (estandarizados dentro de cada firma)
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    rsales_g     = log(saleq) - log(lag(saleq)),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm=TRUE)) / sd(rsales_g, na.rm=TRUE),
    size_std     = (log(atq)   - mean(log(atq), na.rm=TRUE))   / sd(log(atq), na.rm=TRUE),
    sh_current_a_std = (current_ratio - mean(current_ratio, na.rm=TRUE)) / sd(current_ratio, na.rm=TRUE)
  ) %>%
  ungroup()

# 4) Winsorizar y estandarizar las variables financieras
winsorize <- function(x, p=0.005) {
  lo <- quantile(x, p, na.rm=TRUE)
  hi <- quantile(x, 1-p, na.rm=TRUE)
  pmin(pmax(x, lo), hi)
}

df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  mutate(
    lev_dev = leverage_win - mean(leverage_win, na.rm=TRUE),
    dd_dev  = dd_win       - mean(dd_win,       na.rm=TRUE),
    lev_std = lev_dev / sd(lev_dev, na.rm=TRUE),
    dd_std  = dd_dev  / sd(dd_dev,  na.rm=TRUE)
  ) %>%
  ungroup()

# 5) Crear las interacciones de Table 3
df <- df %>%
  mutate(
    # Con dlog_gdp
    lev_std_gdp = lev_std * dlog_gdp,
    dd_std_gdp  = dd_std  * dlog_gdp,
    # Con el shock monetario ortogonalizado (mpso)
    lev_std_wide = lev_std * mpso,
    dd_std_wide  = dd_std  * mpso,
    # Con el choque bruto (mps)
    lev_std_mps  = lev_std * mps,
    dd_std_mps   = dd_std  * mps
  )

# ---------------------------------------------------
# 6) Estimar los cinco modelos de Table 3
#    (efectos fijos: name + Country_x, sin embigl)
# ---------------------------------------------------

# Definir variables de cada modelo
vars_m1 <- c("dlog_capital",
             "lev_std_gdp","lev_std_wide","lev_std_mps")

vars_m2 <- c(vars_m1,
             "rsales_g_std","size_std","sh_current_a_std")

vars_m3 <- c("dlog_capital",
             "dd_std_gdp","dd_std_wide","dd_std_mps",
             "rsales_g_std","size_std","sh_current_a_std")

vars_m4 <- c("dlog_capital",
             "lev_std_gdp","dd_std_gdp",
             "lev_std_wide","lev_std_mps",
             "dd_std_wide","dd_std_mps",
             "rsales_g_std","size_std","sh_current_a_std")

vars_m5 <- c(vars_m4,
             "mpso",
             # lags macro
             paste0("L",1:4,"_dlog_gdp"),
             paste0("L",1:4,"_dlog_cpi"),
             paste0("L",1:4,"_unemp")
)

# Función para preparar dataset y estimar
run_model <- function(vars) {
  df_sub <- df %>% select(name, Country_x, all_of(vars)) %>% drop_na()
  feols(
    as.formula(paste0("dlog_capital ~ ", paste(setdiff(vars, "dlog_capital"), collapse = " + "),
                      " | name + Country_x")),
    data    = df_sub,
    cluster = ~name + Country_x
  )
}

# Estimar modelos M1–M5
m1 <- run_model(vars_m1)
m2 <- run_model(vars_m2)
m3 <- run_model(vars_m3)
m4 <- run_model(vars_m4)
m5 <- run_model(vars_m5)

# 7) Mostrar en consola la Tabla 3
modelsummary(
  list(M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"  = "Leverage × shock (mpso)",
    "dd_std_wide"   = "DD × shock (mpso)",
    "mpso"          = "Monetary shock (wide)",
    "lev_std_gdp"   = "Leverage × GDP growth",
    "dd_std_gdp"    = "DD × GDP growth"
  ),
  gof_omit = ".*"
)



# --- 0) Crear las interacciones "raw" (winsorized & demean) sin estandarizar ---
# (Usamos lev_dev de tu flujo anterior y mpso/mps ya definidos)
df <- df %>%
  mutate(
    lev_wins_dem_wide = lev_dev * mpso,  # “raw” × shock ortogonalizado
    lev_wins_dem      = lev_dev * mps    # “raw” × choque bruto
  )

# Variables de control a nivel firma
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")


# --------------------------------
# Table 7: Empirical results, model vs data
# --------------------------------

## Modelo 1: interacciones estándar (std & demean) + controles
vars_t7_m1 <- c(
  "dlog_capital",
  "lev_std_gdp",    # lev_wins_dem_std_gdp
  "lev_std_wide",   # lev_wins_dem_std_wide
  "lev_std_mps",    # lev_wins_dem_std
  controls_firm
)

df_t7_m1 <- df %>%
  select(name, Country_x, all_of(vars_t7_m1)) %>%
  drop_na()

m7_1 <- feols(
  dlog_capital ~
    lev_std_gdp +
    lev_std_wide +
    lev_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df_t7_m1,
  cluster = ~name + Country_x
)

## Modelo 2: interacciones “raw” (winsorized & demean) + controles
vars_t7_m2 <- c(
  "dlog_capital",
  "lev_std_gdp",       # seguimos controlando GDP growth estandarizado
  "lev_wins_dem_wide",
  "lev_wins_dem",
  controls_firm
)
df_t7_m2 <- df %>%
  select(name, Country_x, all_of(vars_t7_m2)) %>%
  drop_na()

m7_2 <- feols(
  dlog_capital ~
    lev_std_gdp +
    lev_wins_dem_wide +
    lev_wins_dem +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df_t7_m2,
  cluster = ~name + Country_x
)

# Mostrar Table 7
modelsummary(
  list(`(1) Std` = m7_1, `(2) Raw` = m7_2),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"       = "Leverage × shock (wide)",
    "lev_std_mps"        = "Leverage × raw shock",
    "lev_wins_dem_wide"  = "Leverage (dem) × shock (wide)",
    "lev_wins_dem"       = "Leverage (dem) × raw shock"
  ),
  gof_omit = ".*",
  title    = "Table 7: Model vs. Data"
)


# --------------------------------
# Table 9: Heterogeneous responses, no demeaning (con FE firma + país)
# --------------------------------

# 8) Table 9: Heterogeneous responses, no demeaning
df <- df %>%
  group_by(name) %>%
  mutate(
    A_threshold      = median(dd_win, na.rm = TRUE),
    aboveA_dummy     = as.numeric(dd_win >= A_threshold)
  ) %>%
  ungroup() %>%
  mutate(
    lev_nodem_std_gdp  = (leverage_win / sd(leverage_win, na.rm = TRUE)) * dlog_gdp,
    lev_nodem_std_wide = (leverage_win / sd(leverage_win, na.rm = TRUE)) * mpso,
    lev_nodem_std      = (leverage_win / sd(leverage_win, na.rm = TRUE)) * mps,
    d2d_nodem_std_gdp  = (dd_win      / sd(dd_win,       na.rm = TRUE)) * dlog_gdp,
    d2d_nodem_std_wide = (dd_win      / sd(dd_win,       na.rm = TRUE)) * mpso,
    d2d_nodem_std      = (dd_win      / sd(dd_win,       na.rm = TRUE)) * mps,
    aboveA_dummy_gdp   = aboveA_dummy * dlog_gdp,
    aboveA_dummy_wide  = aboveA_dummy * mpso,
    aboveA_dummy_raw   = aboveA_dummy * mps
  )

vars9 <- list(
  m1 = c("dlog_capital", "lev_nodem_std_gdp",  "lev_nodem_std_wide", "lev_nodem_std"),
  m2 = c("dlog_capital", "lev_nodem_std_gdp",  "lev_nodem_std_wide", "lev_nodem_std", controls_firm),
  m3 = c("dlog_capital", "aboveA_dummy_gdp",   "aboveA_dummy_wide",  "aboveA_dummy_raw", controls_firm),
  m4 = c("dlog_capital", "d2d_nodem_std_gdp",  "d2d_nodem_std_wide", "d2d_nodem_std", controls_firm),
  m5 = c("dlog_capital", "lev_nodem_std_gdp",  "aboveA_dummy_gdp",
         "lev_nodem_std_wide","lev_nodem_std", "aboveA_dummy_wide","aboveA_dummy_raw", controls_firm),
  m6 = c("dlog_capital", "lev_nodem_std_gdp",  "d2d_nodem_std_gdp",
         "lev_nodem_std_wide","lev_nodem_std",  "d2d_nodem_std_wide","d2d_nodem_std", controls_firm),
  m7 = c("dlog_capital", "mpso", "lev_nodem_std_gdp","d2d_nodem_std_gdp",
         "lev_nodem_std_wide","lev_nodem_std","d2d_nodem_std_wide","d2d_nodem_std",
         controls_firm,
         paste0("L",1:4,"_dlog_gdp"),
         paste0("L",1:4,"_dlog_cpi"),
         paste0("L",1:4,"_unemp")
  )
)

models9 <- list()
for (i in seq_along(vars9)) {
  vars_i <- vars9[[i]]
  df_i   <- df %>% select(name, Country_x, all_of(vars_i)) %>% drop_na()
  f_i    <- as.formula(paste0(
    "dlog_capital ~ ",
    paste(setdiff(vars_i, "dlog_capital"), collapse = " + "),
    " | name + Country_x"
  ))
  models9[[i]] <- feols(f_i, data = df_i, cluster = ~name + Country_x)
}

modelsummary(
  models9,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_nodem_std_wide"  = "Leverage (nodem,std) × shock (mpso)",
    "lev_nodem_std"       = "Leverage (nodem,std) × raw shock",
    "lev_nodem_std_gdp"   = "Leverage (nodem,std) × GDP",
    "d2d_nodem_std_wide"  = "DD (nodem,std) × shock (mpso)",
    "d2d_nodem_std"       = "DD (nodem,std) × raw shock",
    "d2d_nodem_std_gdp"   = "DD (nodem,std) × GDP",
    "aboveA_dummy_gdp"    = "Above-A × GDP",
    "aboveA_dummy_wide"   = "Above-A × shock (mpso)",
    "aboveA_dummy_raw"    = "Above-A × raw shock"
  ),
  gof_omit = ".*",
  title    = "Table 9: Heterogeneous responses (no demeaning)"
)

  


# =======================================
# Table 10: Main results, no cyclical sensitivities (sin embigl, con FE firma + país)
# =======================================

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# Controles firm-level (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# Controles agregados (lags macro)
controls_agg <- c(
  "L1_dlog_gdp", "L2_dlog_gdp", "L3_dlog_gdp", "L4_dlog_gdp",
  "L1_dlog_cpi", "L2_dlog_cpi", "L3_dlog_cpi", "L4_dlog_cpi",
  "L1_unemp",    "L2_unemp",    "L3_unemp",    "L4_unemp"
)

# 1) Modelo 1: leverage only (interacciones estandarizadas)
vars10_m1 <- c(
  "dlog_capital",
  "lev_std_wide",   # Leverage × shock (mpso)
  "lev_std_mps",    # Leverage × raw shock
  controls_firm
)
df10_m1 <- df %>%
  select(name, Country_x, all_of(vars10_m1)) %>%
  drop_na()

m10_1 <- feols(
  dlog_capital ~ 
    lev_std_wide +   # Leverage × shock (mpso)
    lev_std_mps +    # Leverage × raw shock
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df10_m1,
  cluster = ~name + Country_x
)

# 2) Modelo 2: distance-to-default only (interacciones estandarizadas)
vars10_m2 <- c(
  "dlog_capital",
  "dd_std_wide",  # DD × shock (mpso)
  "dd_std_mps",   # DD × raw shock
  controls_firm
)
df10_m2 <- df %>%
  select(name, Country_x, all_of(vars10_m2)) %>%
  drop_na()

m10_2 <- feols(
  dlog_capital ~ 
    dd_std_wide +   # DD × shock (mpso)
    dd_std_mps +    # DD × raw shock
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df10_m2,
  cluster = ~name + Country_x
)

# 3) Modelo 3: leverage + dd (ambos estandarizados)
vars10_m3 <- c(
  "dlog_capital",
  "lev_std_wide", "lev_std_mps",
  "dd_std_wide",  "dd_std_mps",
  controls_firm
)
df10_m3 <- df %>%
  select(name, Country_x, all_of(vars10_m3)) %>%
  drop_na()

m10_3 <- feols(
  dlog_capital ~ 
    lev_std_wide + lev_std_mps +
    dd_std_wide  + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df10_m3,
  cluster = ~name + Country_x
)

# 4) Modelo 4: incluye choque puro (mpso) y controles agregados
vars10_m4 <- c(
  "dlog_capital",
  "mpso",             # Monetary shock (wide)
  "lev_std_wide",     # Leverage × shock (mpso)
  "lev_std_mps",      # Leverage × raw shock
  "dd_std_wide",      # DD × shock (mpso)
  "dd_std_mps",       # DD × raw shock
  controls_firm,
  controls_agg
)
df10_m4 <- df %>%
  select(name, Country_x, all_of(vars10_m4)) %>%
  drop_na()

m10_4 <- feols(
  dlog_capital ~ 
    mpso +
    lev_std_wide + lev_std_mps +
    dd_std_wide  + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std +
    L1_dlog_gdp + L2_dlog_gdp + L3_dlog_gdp + L4_dlog_gdp +
    L1_dlog_cpi + L2_dlog_cpi + L3_dlog_cpi + L4_dlog_cpi +
    L1_unemp    + L2_unemp    + L3_unemp    + L4_unemp
  | name + Country_x,
  data    = df10_m4,
  cluster = ~name + Country_x
)

# Mostrar Table 10
modelsummary(
  list(`M1` = m10_1, `M2` = m10_2, `M3` = m10_3, `M4` = m10_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"   = "Leverage × shock (mpso)",
    "lev_std_mps"    = "Leverage × raw shock",
    "dd_std_wide"    = "DD × shock (mpso)",
    "dd_std_mps"     = "DD × raw shock",
    "mpso"           = "Monetary shock (wide)"
  ),
  gof_omit = ".*",
  title    = "Table 10: Main results (no cyclical sensitivities)"
)


# =======================================
# Table 11: Expansionary vs. Contractionary (sin embigl, con FE firma + país)
# =======================================

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# 1) Creamos choques asimétricos basados en mps (raw shock)
df <- df %>%
  mutate(
    # Leverage × shock (mpso) asimétrico
    lev_std_wide_p   = if_else(mps > 0, lev_std_wide, 0),
    lev_std_wide_n   = if_else(mps < 0, lev_std_wide, 0),
    # DD × shock (mpso) asimétrico
    dd_std_wide_p    = if_else(mps > 0, dd_std_wide, 0),
    dd_std_wide_n    = if_else(mps < 0, dd_std_wide, 0)
  )

# Variables de control a nivel firma (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 1) Modelo 1: leverage simétrico
vars11_m1 <- c(
  "dlog_capital",
  "lev_std_gdp",     # Leverage × GDP growth
  "lev_std_wide",    # Leverage × shock (mpso)
  "lev_std_mps",     # Leverage × raw shock
  controls_firm
)
df11_m1 <- df %>%
  select(name, Country_x, all_of(vars11_m1)) %>%
  drop_na()

m11_1 <- feols(
  dlog_capital ~ 
    lev_std_gdp + 
    lev_std_wide + 
    lev_std_mps + 
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df11_m1,
  cluster = ~name + Country_x
)

# 2) Modelo 2: leverage pos/neg
vars11_m2 <- c(
  "dlog_capital",
  "lev_std_gdp",       # Leverage × GDP growth
  "lev_std_wide_p",    # Leverage × shock (mpso) (+)
  "lev_std_wide_n",    # Leverage × shock (mpso) (–)
  "lev_std_mps",       # Leverage × raw shock
  controls_firm
)
df11_m2 <- df %>%
  select(name, Country_x, all_of(vars11_m2)) %>%
  drop_na()

m11_2 <- feols(
  dlog_capital ~ 
    lev_std_gdp + 
    lev_std_wide_p + 
    lev_std_wide_n + 
    lev_std_mps + 
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df11_m2,
  cluster = ~name + Country_x
)

# 3) Modelo 3: dd simétrico
vars11_m3 <- c(
  "dlog_capital",
  "dd_std_gdp",      # DD × GDP growth
  "dd_std_wide",     # DD × shock (mpso)
  "dd_std_mps",      # DD × raw shock
  controls_firm
)
df11_m3 <- df %>%
  select(name, Country_x, all_of(vars11_m3)) %>%
  drop_na()

m11_3 <- feols(
  dlog_capital ~ 
    dd_std_gdp + 
    dd_std_wide + 
    dd_std_mps + 
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df11_m3,
  cluster = ~name + Country_x
)

# 4) Modelo 4: dd pos/neg
vars11_m4 <- c(
  "dlog_capital",
  "dd_std_gdp",      # DD × GDP growth
  "dd_std_wide_p",   # DD × shock (mpso) (+)
  "dd_std_wide_n",   # DD × shock (mpso) (–)
  "dd_std_mps",      # DD × raw shock
  controls_firm
)
df11_m4 <- df %>%
  select(name, Country_x, all_of(vars11_m4)) %>%
  drop_na()

m11_4 <- feols(
  dlog_capital ~ 
    dd_std_gdp + 
    dd_std_wide_p + 
    dd_std_wide_n + 
    dd_std_mps + 
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df11_m4,
  cluster = ~name + Country_x
)

# Mostrar Table 11
modelsummary(
  list(`(1)` = m11_1, `(2)` = m11_2, `(3)` = m11_3, `(4)` = m11_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"      = "Leverage × shock (mpso)",
    "lev_std_wide_p"    = "Leverage × shock (mpso) (+)",
    "lev_std_wide_n"    = "Leverage × shock (mpso) (–)",
    "lev_std_gdp"       = "Leverage × GDP growth",
    "lev_std_mps"       = "Leverage × raw shock",
    "dd_std_wide"       = "DD × shock (mpso)",
    "dd_std_wide_p"     = "DD × shock (mpso) (+)",
    "dd_std_wide_n"     = "DD × shock (mpso) (–)",
    "dd_std_gdp"        = "DD × GDP growth",
    "dd_std_mps"        = "DD × raw shock"
  ),
  gof_omit = ".*",
  title    = "Table 11: Expansionary vs. Contractionary (sin embigl)"
)



# =======================================
# Table 15: Lagged Investment (sin embigl, con FE firma + país)
# =======================================

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# 0) Crear Ldl_capital (lag de dlog_capital) por firma
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    Ldl_capital = lag(dlog_capital, 1)
  ) %>%
  ungroup()

# 1) Definir controles
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")
controls_agg  <- c(
  "L1_dlog_gdp", "L2_dlog_gdp", "L3_dlog_gdp", "L4_dlog_gdp",
  "L1_dlog_cpi", "L2_dlog_cpi", "L3_dlog_cpi", "L4_dlog_cpi",
  "L1_unemp",    "L2_unemp",    "L3_unemp",    "L4_unemp"
)

# 2) Modelo 1: leverage only + lagged inv (interacciones estandarizadas)
vars15_m1 <- c(
  "dlog_capital",
  "lev_std_wide",  # Leverage × shock (mpso)
  "lev_std_mps",   # Leverage × raw shock
  "Ldl_capital",
  controls_firm
)
df15_m1 <- df %>%
  select(name, Country_x, all_of(vars15_m1)) %>%
  drop_na()

m15_1 <- feols(
  dlog_capital ~ 
    lev_std_wide + lev_std_mps +
    Ldl_capital +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df15_m1,
  cluster = ~name + Country_x
)

# 3) Modelo 2: dd only + lagged inv (interacciones estandarizadas)
vars15_m2 <- c(
  "dlog_capital",
  "dd_std_wide",  # DD × shock (mpso)
  "dd_std_mps",   # DD × raw shock
  "Ldl_capital",
  controls_firm
)
df15_m2 <- df %>%
  select(name, Country_x, all_of(vars15_m2)) %>%
  drop_na()

m15_2 <- feols(
  dlog_capital ~ 
    dd_std_wide + dd_std_mps +
    Ldl_capital +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df15_m2,
  cluster = ~name + Country_x
)

# 4) Modelo 3: leverage + dd + lagged inv (estandarizadas)
vars15_m3 <- c(
  "dlog_capital",
  "lev_std_wide", "lev_std_mps",
  "dd_std_wide",  "dd_std_mps",
  "Ldl_capital",
  controls_firm
)
df15_m3 <- df %>%
  select(name, Country_x, all_of(vars15_m3)) %>%
  drop_na()

m15_3 <- feols(
  dlog_capital ~ 
    lev_std_wide + lev_std_mps +
    dd_std_wide  + dd_std_mps +
    Ldl_capital +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df15_m3,
  cluster = ~name + Country_x
)

# 5) Modelo 4: leverage + dd + agregados + lagged inv
vars15_m4 <- c(
  "dlog_capital",
  "lev_std_wide", "lev_std_mps",
  "dd_std_wide",  "dd_std_mps",
  "Ldl_capital",
  controls_firm, controls_agg
)
df15_m4 <- df %>%
  select(name, Country_x, all_of(vars15_m4)) %>%
  drop_na()

m15_4 <- feols(
  dlog_capital ~ 
    lev_std_wide + lev_std_mps +
    dd_std_wide  + dd_std_mps +
    Ldl_capital +
    rsales_g_std + size_std + sh_current_a_std +
    L1_dlog_gdp + L2_dlog_gdp + L3_dlog_gdp + L4_dlog_gdp +
    L1_dlog_cpi + L2_dlog_cpi + L3_dlog_cpi + L4_dlog_cpi +
    L1_unemp    + L2_unemp    + L3_unemp    + L4_unemp
  | name + Country_x,
  data    = df15_m4,
  cluster = ~name + Country_x
)

# 6) Mostrar Table 15
modelsummary(
  list(M1 = m15_1, M2 = m15_2, M3 = m15_3, M4 = m15_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"   = "Leverage × shock (mpso)",
    "lev_std_mps"    = "Leverage × raw shock",
    "dd_std_wide"    = "DD × shock (mpso)",
    "dd_std_mps"     = "DD × raw shock",
    "Ldl_capital"    = "Lagged Δ log k"
  ),
  keep     = c("lev_std_wide", "Ldl_capital", "dd_std_wide"),
  order    = c("lev_std_wide", "Ldl_capital", "dd_std_wide"),
  gof_omit = ".*",
  title    = "Table 15: Lagged Investment (sin riesgo soberano)"
)



# ----------------------------------------------------------------
# Table 18: Controlling for differences in cyclical sensitivities
# ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# 1) Crear interacciones de sensibilidad a CPI y desempleo
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    lev_std_cpi = lev_std * dlog_cpi,  # Leverage × inflación
    dd_std_cpi  = dd_std  * dlog_cpi,  # DD × inflación
    lev_std_ur  = lev_std * unemp,     # Leverage × desempleo
    dd_std_ur   = dd_std  * unemp      # DD × desempleo
  ) %>%
  ungroup()

# Controles firm‐level (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 1) Modelo 1: sensibilidad al GDP (interacciones estandarizadas previas)
vars18_m1 <- c(
  "dlog_capital",
  "lev_std_gdp",    # Leverage × GDP growth
  "lev_std_wide",   # Leverage × shock (mpso)
  "lev_std_mps",    # Leverage × raw shock
  controls_firm
)
df18_m1 <- df %>%
  select(name, Country_x, all_of(vars18_m1)) %>%
  drop_na()

m18_1 <- feols(
  dlog_capital ~
    lev_std_gdp + lev_std_wide + lev_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df18_m1,
  cluster = ~name + Country_x
)

# 2) Modelo 2: sensibilidad a DD (interacciones estandarizadas previas)
vars18_m2 <- c(
  "dlog_capital",
  "dd_std_gdp",   # DD × GDP growth
  "dd_std_wide",  # DD × shock (mpso)
  "dd_std_mps",   # DD × raw shock
  controls_firm
)
df18_m2 <- df %>%
  select(name, Country_x, all_of(vars18_m2)) %>%
  drop_na()

m18_2 <- feols(
  dlog_capital ~
    dd_std_gdp + dd_std_wide + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df18_m2,
  cluster = ~name + Country_x
)

# 3) Modelo 3: sensibilidad a CPI (Leverage)
vars18_m3 <- c(
  "dlog_capital",
  "lev_std_gdp",    # Leverage × GDP growth
  "lev_std_cpi",    # Leverage × CPI growth
  "lev_std_wide",   # Leverage × shock (mpso)
  "lev_std_mps",    # Leverage × raw shock
  controls_firm
)
df18_m3 <- df %>%
  select(name, Country_x, all_of(vars18_m3)) %>%
  drop_na()

m18_3 <- feols(
  dlog_capital ~
    lev_std_gdp + lev_std_cpi + lev_std_wide + lev_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df18_m3,
  cluster = ~name + Country_x
)

# 4) Modelo 4: sensibilidad a DD + CPI
vars18_m4 <- c(
  "dlog_capital",
  "dd_std_gdp",   # DD × GDP growth
  "dd_std_cpi",   # DD × CPI growth
  "dd_std_wide",  # DD × shock (mpso)
  "dd_std_mps",   # DD × raw shock
  controls_firm
)
df18_m4 <- df %>%
  select(name, Country_x, all_of(vars18_m4)) %>%
  drop_na()

m18_4 <- feols(
  dlog_capital ~
    dd_std_gdp + dd_std_cpi + dd_std_wide + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df18_m4,
  cluster = ~name + Country_x
)

# 5) Modelo 5: sensibilidad a desempleo (Leverage)
vars18_m5 <- c(
  "dlog_capital",
  "lev_std_ur",    # Leverage × unemployment
  "lev_std_wide",  # Leverage × shock (mpso)
  "lev_std_mps",   # Leverage × raw shock
  controls_firm
)
df18_m5 <- df %>%
  select(name, Country_x, all_of(vars18_m5)) %>%
  drop_na()

m18_5 <- feols(
  dlog_capital ~
    lev_std_ur + lev_std_wide + lev_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df18_m5,
  cluster = ~name + Country_x
)

# 6) Modelo 6: sensibilidad a desempleo (DD)
vars18_m6 <- c(
  "dlog_capital",
  "dd_std_ur",     # DD × unemployment
  "dd_std_wide",   # DD × shock (mpso)
  "dd_std_mps",    # DD × raw shock
  controls_firm
)
df18_m6 <- df %>%
  select(name, Country_x, all_of(vars18_m6)) %>%
  drop_na()

m18_6 <- feols(
  dlog_capital ~
    dd_std_ur + dd_std_wide + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df18_m6,
  cluster = ~name + Country_x
)

# Mostrar Table 18
modelsummary(
  list(M1 = m18_1, M2 = m18_2, M3 = m18_3, M4 = m18_4, M5 = m18_5, M6 = m18_6),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"     = "Leverage × shock (mpso)",
    "dd_std_wide"      = "DD × shock (mpso)",
    "lev_std_gdp"      = "Leverage × GDP growth",
    "dd_std_gdp"       = "DD × GDP growth",
    "lev_std_cpi"      = "Leverage × CPI growth",
    "dd_std_cpi"       = "DD × CPI growth",
    "lev_std_ur"       = "Leverage × unemployment",
    "dd_std_ur"        = "DD × unemployment"
  ),
  keep     = c(
    "lev_std_wide", "dd_std_wide",
    "lev_std_gdp",  "dd_std_gdp",
    "lev_std_cpi",  "dd_std_cpi",
    "lev_std_ur",   "dd_std_ur"
  ),
  order    = c(
    "lev_std_wide", "dd_std_wide",
    "lev_std_gdp",  "dd_std_gdp",
    "lev_std_cpi",  "dd_std_cpi",
    "lev_std_ur",   "dd_std_ur"
  ),
  gof_omit = ".*",
  title    = "Table 18: Controlling for differences in cyclical sensitivities"
)




# =======================================
# Table 19: Target vs. Path Decomposition (sin embigl, con FE firma + país)
# =======================================

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# ————————
# 0) Creamos las interacciones Target/Path
# ————————
df <- df %>%
  mutate(
    # Interacciones estándar
    lev_std_gdp              = lev_std * dlog_gdp,
    lev_std_wide             = lev_std * mpso,     # en lugar de "wide"
    lev_std_mps              = lev_std * mps,      # en lugar de "mps_raw"
    # Target = componente positiva del choque
    lev_std_target           = lev_std * shock_pos,
    dd_std_target            = dd_std  * shock_pos,
    # Path = componente negativa del choque
    lev_std_path             = lev_std * shock_neg,
    dd_std_path              = dd_std  * shock_neg
  )

# ————————
# 1) Definir controles firm-level (sin embigl)
# ————————
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# ————————
# 2) Table 19: Target vs. Path Decomposition
# ————————

# Modelo 1: leverage × (GDP, wide, raw)
m19_1 <- feols(
  dlog_capital ~
    lev_std_gdp + lev_std_wide + lev_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(lev_std_gdp, lev_std_wide, lev_std_mps, rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# Modelo 2: leverage × (GDP, target, path, raw)
m19_2 <- feols(
  dlog_capital ~
    lev_std_gdp + lev_std_target + lev_std_path + lev_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(lev_std_gdp, lev_std_target, lev_std_path, lev_std_mps, rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# Modelo 3: DD × (GDP, wide, raw)
m19_3 <- feols(
  dlog_capital ~
    dd_std_gdp + dd_std_wide + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(dd_std_gdp, dd_std_wide, dd_std_mps, rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# Modelo 4: DD × (GDP, target, path, raw)
m19_4 <- feols(
  dlog_capital ~
    dd_std_gdp + dd_std_target + dd_std_path + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(dd_std_gdp, dd_std_target, dd_std_path, dd_std_mps, rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# ————————
# 3) Mostrar Table 19
# ————————
modelsummary(
  list(`(1)` = m19_1, `(2)` = m19_2, `(3)` = m19_3, `(4)` = m19_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_wide"       = "Leverage × shock (wide)",
    "lev_std_target"     = "Leverage × shock (target)",
    "lev_std_path"       = "Leverage × shock (path)",
    "lev_std_gdp"        = "Leverage × GDP growth",
    "lev_std_mps"        = "Leverage × raw shock",
    "dd_std_wide"        = "DD × shock (wide)",
    "dd_std_target"      = "DD × shock (target)",
    "dd_std_path"        = "DD × shock (path)",
    "dd_std_gdp"         = "DD × GDP growth",
    "dd_std_mps"         = "DD × raw shock"
  ),
  keep     = c(
    "lev_std_wide", "lev_std_target", "lev_std_path",
    "dd_std_wide",  "dd_std_target",  "dd_std_path"
  ),
  order    = c(
    "lev_std_wide", "lev_std_target", "lev_std_path",
    "dd_std_wide",  "dd_std_target",  "dd_std_path"
  ),
  gof_omit = ".*",
  title    = "Table 19: Target vs. Path Decomposition (sin embigl)"
)





# ------------------------------------------------
# Table 20: Alternative Time Aggregation (sin embigl, con FE firma + país)
# ------------------------------------------------

library(dplyr)
library(fixest)
library(modelsummary)

# 0) Crear wide_sum (suma de shocks mps en rezagos 0–4) y nuevas interacciones
df <- df %>%
  mutate(
    wide_sum    = L0_mps + L1_mps + L2_mps + L3_mps + L4_mps,
    lev_std_sum = lev_std * wide_sum,  # Leverage × shock (sum)
    dd_std_sum  = dd_std  * wide_sum   # DD × shock (sum)
  )

# 1) Definir controles (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 2) Modelos

# M1: Leverage only
m20_1 <- feols(
  dlog_capital ~ 
    lev_std_gdp +         # Leverage × GDP growth
    lev_std_sum +         # Leverage × shock (sum)
    lev_std_mps +         # Leverage × raw shock
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(lev_std_gdp, lev_std_sum, lev_std_mps,
                           rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# M2: DD only
m20_2 <- feols(
  dlog_capital ~ 
    dd_std_gdp +          # DD × GDP growth
    dd_std_sum +          # DD × shock (sum)
    dd_std_mps +          # DD × raw shock
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(dd_std_gdp, dd_std_sum, dd_std_mps,
                           rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# M3: Leverage + DD
m20_3 <- feols(
  dlog_capital ~ 
    lev_std_gdp + dd_std_gdp +
    lev_std_sum + lev_std_mps +
    dd_std_sum  + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std
  | name + Country_x,
  data    = df %>% drop_na(lev_std_gdp, dd_std_gdp,
                           lev_std_sum, lev_std_mps,
                           dd_std_sum,  dd_std_mps,
                           rsales_g_std, size_std, sh_current_a_std),
  cluster = ~name + Country_x
)

# M4: wide_sum + controles agregados
m20_4 <- feols(
  dlog_capital ~ 
    lev_std_gdp + 
    wide_sum +
    lev_std_sum + lev_std_mps +
    dd_std_sum  + dd_std_mps +
    rsales_g_std + size_std + sh_current_a_std +
    L1_dlog_gdp + L2_dlog_gdp + L3_dlog_gdp + L4_dlog_gdp +
    L1_dlog_cpi + L2_dlog_cpi + L3_dlog_cpi + L4_dlog_cpi +
    L1_unemp    + L2_unemp    + L3_unemp    + L4_unemp
  | name + Country_x,
  data    = df %>% drop_na(lev_std_gdp, wide_sum, lev_std_sum, lev_std_mps,
                           dd_std_sum,  dd_std_mps,
                           rsales_g_std, size_std, sh_current_a_std,
                           L1_dlog_gdp:L4_unemp),
  cluster = ~name + Country_x
)

# 3) Mostrar resultados
modelsummary(
  list(M1 = m20_1, M2 = m20_2, M3 = m20_3, M4 = m20_4),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  coef_map = c(
    "lev_std_sum"      = "Leverage × shock (sum)",
    "dd_std_sum"       = "DD × shock (sum)",
    "wide_sum"         = "Sum of shocks",
    "lev_std_gdp"      = "Leverage × GDP growth",
    "dd_std_gdp"       = "DD × GDP growth",
    "lev_std_mps"      = "Leverage × raw shock",
    "dd_std_mps"       = "DD × raw shock"
  ),
  keep     = c("lev_std_sum", "dd_std_sum", "wide_sum"),
  gof_omit = ".*",
  title    = "Table 20: Alternative Time Aggregation (sin embigl)"
)





# --------------------------------------------------------
# Table 21: Interaction with Other Firm-Level Covariates (sin embigl, con FE firma + país)
# --------------------------------------------------------

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# 0) Crear variables estandarizadas y de crecimiento futuro de ventas
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Estandarizar leverage y dd
    lev_std = (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
    dd_std  = (dd       - mean(dd,        na.rm = TRUE)) / sd(dd,        na.rm = TRUE),
    # Crecimiento futuro de ventas (4 trimestres adelante)
    rsales_Fg4     = log(lead(saleq, 4)) - log(lead(saleq, 3)),
    rsales_Fg4_std = (rsales_Fg4 - mean(rsales_Fg4, na.rm = TRUE)) / sd(rsales_Fg4, na.rm = TRUE)
  ) %>%
  ungroup()

# 1) Crear interacciones con 'mpso' (shock ortogonalizado) y 'mps' (raw shock)
df <- df %>%
  mutate(
    # Interacciones "standardized & demean"
    lev_std_gdp             = lev_std        * dlog_gdp,      # Leverage × GDP growth
    dd_std_gdp              = dd_std         * dlog_gdp,      # DD × GDP growth
    lev_std_wide            = lev_std        * mpso,          # Leverage × shock (wide)
    dd_std_wide             = dd_std         * mpso,          # DD × shock (wide)
    lev_std_mps             = lev_std        * mps,           # Leverage × raw shock
    dd_std_mps              = dd_std         * mps,           # DD × raw shock
    
    # Ventanas de ventas y tamaño con shock (wide)
    L0rsales_g_std_wide     = rsales_g_std   * mpso,          # Sales growth × shock (wide)
    L0rsales_Fg4_std_wide   = rsales_Fg4_std * mpso,          # Future sales growth × shock (wide)
    L0size_std_wide         = size_std       * mpso,          # Size × shock (wide)
    liq_std_wide            = sh_current_a_std * mpso,       # Liquidity × shock (wide)
    
    # Interacciones "raw" (winsorized & demean) con raw shock mps
    lev_std_mps_raw         = lev_std        * mps,           # Leverage (dem) × raw shock
    dd_std_mps_raw          = dd_std         * mps,           # DD (dem) × raw shock
    liq_std_mps_raw         = sh_current_a_std * mps         # Liquidity (dem) × raw shock
  )

# 2) Definir controles firm-level (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 3) Listas de variables para cada uno de los 8 modelos
vars21 <- list(
  m1 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide",
         "L0rsales_g_std_wide", "lev_std_mps",
         controls_firm),
  m2 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide",
         "L0rsales_g_std_wide", "dd_std_mps",
         controls_firm),
  m3 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide",
         "L0rsales_Fg4_std_wide", "lev_std_mps", "rsales_Fg4_std",
         controls_firm),
  m4 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide",
         "L0rsales_Fg4_std_wide", "dd_std_mps", "rsales_Fg4_std",
         controls_firm),
  m5 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide",
         "L0size_std_wide", "lev_std_mps",
         controls_firm),
  m6 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide",
         "L0size_std_wide", "dd_std_mps",
         controls_firm),
  m7 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide",
         "lev_std_mps", "liq_std_wide", "liq_std_mps_raw",
         controls_firm),
  m8 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide",
         "dd_std_mps", "liq_std_wide", "liq_std_mps_raw",
         controls_firm)
)

# 4) Estimar los 8 modelos con FE por firma y país, cluster por firma y país
models21 <- lapply(vars21, function(vars) {
  df_sub  <- df %>% select(name, Country_x, all_of(vars)) %>% drop_na()
  formula <- as.formula(paste("dlog_capital ~",
                              paste(setdiff(vars, "dlog_capital"), collapse = " + "),
                              "| name + Country_x"))
  feols(formula, data = df_sub, cluster = ~name + Country_x)
})

# 5) Mostrar Table 21
modelsummary(
  models21,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_std_wide", "dd_std_wide",
               "L0rsales_g_std_wide", "L0rsales_Fg4_std_wide",
               "L0size_std_wide", "liq_std_wide"),
  order    = c("lev_std_wide", "dd_std_wide",
               "L0rsales_g_std_wide", "L0rsales_Fg4_std_wide",
               "L0size_std_wide", "liq_std_wide"),
  coef_map = c(
    "lev_std_wide"          = "Leverage × shock (wide)",
    "dd_std_wide"           = "DD × shock (wide)",
    "L0rsales_g_std_wide"   = "Sales growth × shock (wide)",
    "L0rsales_Fg4_std_wide" = "Future sales growth × shock (wide)",
    "L0size_std_wide"       = "Size × shock (wide)",
    "liq_std_wide"          = "Liquidity × shock (wide)"
  ),
  gof_omit = ".*",
  title    = "Table 21: Interaction with Other Firm-Level Covariates"
)








# --------------------------------------------------------
# Table 22: Interaction with Other Measures of Financial Positions (sin embigl, con FE firma + país)
# --------------------------------------------------------

library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)

# 0) Estandarizar y preparar variables
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Estandarizaciones básicas
    lev_std        = (leverage       - mean(leverage,       na.rm = TRUE)) / sd(leverage,       na.rm = TRUE),
    dd_std         = (dd             - mean(dd,             na.rm = TRUE)) / sd(dd,             na.rm = TRUE),
    size_std       = (log_size       - mean(log_size,       na.rm = TRUE)) / sd(log_size,       na.rm = TRUE),
    rsales_g_std   = dln_saleq,
    sh_current_a_std = current_ratio,
    # Flujo de caja operativo
    cashop_std     = (cashop         - mean(cashop,         na.rm = TRUE)) / sd(cashop,         na.rm = TRUE)
  ) %>%
  ungroup()

# 1) Crear interacciones con 'mpso' (shock ortogonalizado) y 'mps' (raw shock)
df <- df %>%
  mutate(
    # Columnas 1–2: tamaño
    L0size_std_wide        = size_std       * mpso,
    # Columnas 3–4: flujo de caja
    cashop_wide            = cashop_std     * mpso,
    cashop_raw             = cashop_std     * mps,
    # Columnas 5–6: dividendos (dummy 'dividends')
    L0div_wide             = dividends      * mpso,
    div_raw                = dividends      * mps,
    # Columnas 7–8: liquidez
    liq_wide               = sh_current_a_std * mpso,
    liq_raw                = sh_current_a_std * mps,
    
    # Columnas 9–10: interactuar leverage
    lev_std_gdp            = lev_std        * dlog_gdp,
    lev_std_wide           = lev_std        * mpso,
    lev_std_mps            = lev_std        * mps,
    
    # Columnas 11–12: interactuar DD
    dd_std_gdp             = dd_std         * dlog_gdp,
    dd_std_wide            = dd_std         * mpso,
    dd_std_mps             = dd_std         * mps
  )

# 2) Controles a nivel firma (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 3) Listas de variables para cada modelo
vars22 <- list(
  m1 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide", "lev_std_mps",
         "L0size_std_wide",
         controls_firm),
  m2 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide", "dd_std_mps",
         "L0size_std_wide",
         controls_firm),
  m3 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide", "lev_std_mps",
         "cashop_wide", "cashop_raw",
         controls_firm),
  m4 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide", "dd_std_mps",
         "cashop_wide", "cashop_raw",
         controls_firm),
  m5 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide", "lev_std_mps",
         "L0div_wide", "div_raw",
         controls_firm),
  m6 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide", "dd_std_mps",
         "L0div_wide", "div_raw",
         controls_firm),
  m7 = c("dlog_capital",
         "lev_std_gdp", "lev_std_wide", "lev_std_mps",
         "liq_wide", "liq_raw",
         controls_firm),
  m8 = c("dlog_capital",
         "dd_std_gdp", "dd_std_wide", "dd_std_mps",
         "liq_wide", "liq_raw",
         controls_firm)
)

# 4) Estimar los 8 modelos con FE por firma y país, cluster por firma y país
models22 <- lapply(vars22, function(vs) {
  df_sub <- df %>% select(name, Country_x, all_of(vs)) %>% drop_na()
  fml    <- as.formula(paste("dlog_capital ~",
                             paste(setdiff(vs, "dlog_capital"), collapse = " + "),
                             "| name + Country_x"))
  feols(fml, data = df_sub, cluster = ~name + Country_x)
})

# 5) Mostrar Table 22
modelsummary(
  models22,
  stars = c('*' = .10, '**' = .05, '***' = .01),
  keep  = c("lev_std_wide", "dd_std_wide",
            "L0size_std_wide", "cashop_wide", "L0div_wide", "liq_wide"),
  order = c("lev_std_wide", "dd_std_wide",
            "L0size_std_wide", "cashop_wide", "L0div_wide", "liq_wide"),
  coef_map = c(
    "lev_std_wide"           = "Leverage × shock (wide)",
    "dd_std_wide"            = "DD × shock (wide)",
    "L0size_std_wide"        = "Size × shock (wide)",
    "cashop_wide"            = "Cash flows × shock (wide)",
    "L0div_wide"             = "Dividend indicator × shock (wide)",
    "liq_wide"               = "Liquidity × shock (wide)"
  ),
  gof_omit = ".*",
  title    = "Table 22: Interaction with Other Measures of Financial Positions"
)






# -------------------------------------------------------------------------
# Table 23: IV usando posiciones financieras rezagadas (sin embigl, con FE firma + país)
# -------------------------------------------------------------------------

library(dplyr)
library(fixest)
library(modelsummary)

# 0) Crear lags de las interacciones endógenas (usando nombres estandarizados)
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Interacciones estándar para leverage
    lev_std_gdp    = lev_std * dlog_gdp,
    lev_std_wide   = lev_std * mpso,
    lev_std_mps    = lev_std * mps,
    
    # Interacciones estándar para DD
    dd_std_gdp     = dd_std  * dlog_gdp,
    dd_std_wide    = dd_std  * mpso,
    dd_std_mps     = dd_std  * mps,
    
    # Lags para leverage
    levL1_std_gdp  = lag(lev_std_gdp, 1),
    levL1_std_wide = lag(lev_std_wide, 1),
    levL1_std_mps  = lag(lev_std_mps, 1),
    levL2_std_gdp  = lag(lev_std_gdp, 2),
    levL2_std_wide = lag(lev_std_wide, 2),
    levL2_std_mps  = lag(lev_std_mps, 2),
    levL4_std_gdp  = lag(lev_std_gdp, 4),
    levL4_std_wide = lag(lev_std_wide, 4),
    levL4_std_mps  = lag(lev_std_mps, 4),
    
    # Lags para DD
    ddL1_std_gdp   = lag(dd_std_gdp, 1),
    ddL1_std_wide  = lag(dd_std_wide, 1),
    ddL1_std_mps   = lag(dd_std_mps, 1),
    ddL2_std_gdp   = lag(dd_std_gdp, 2),
    ddL2_std_wide  = lag(dd_std_wide, 2),
    ddL2_std_mps   = lag(dd_std_mps, 2),
    ddL4_std_gdp   = lag(dd_std_gdp, 4),
    ddL4_std_wide  = lag(dd_std_wide, 4),
    ddL4_std_mps   = lag(dd_std_mps, 4)
  ) %>%
  ungroup()

# 1) Controles firm-level (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 2) Modelos IV con filtrado previo para evitar NAs excesivos

## Model 1: Leverage instrumentado con lag 1
df23_1 <- df %>%
  drop_na(
    dlog_capital,
    rsales_g_std, size_std, sh_current_a_std,
    lev_std_gdp, lev_std_wide, lev_std_mps,
    levL1_std_gdp, levL1_std_wide, levL1_std_mps
  )

m23_1 <- feols(
  dlog_capital ~ rsales_g_std + size_std + sh_current_a_std |
    name + Country_x |
    (lev_std_gdp + lev_std_wide + lev_std_mps) ~
    (levL1_std_gdp + levL1_std_wide + levL1_std_mps),
  data    = df23_1,
  cluster = ~name + Country_x
)

## Model 2: Leverage instrumentado con lag 2
df23_2 <- df %>%
  drop_na(
    dlog_capital,
    rsales_g_std, size_std, sh_current_a_std,
    lev_std_gdp, lev_std_wide, lev_std_mps,
    levL2_std_gdp, levL2_std_wide, levL2_std_mps
  )

m23_2 <- feols(
  dlog_capital ~ rsales_g_std + size_std + sh_current_a_std |
    name + Country_x |
    (lev_std_gdp + lev_std_wide + lev_std_mps) ~
    (levL2_std_gdp + levL2_std_wide + levL2_std_mps),
  data    = df23_2,
  cluster = ~name + Country_x
)

## Model 3: Leverage instrumentado con lag 4
df23_3 <- df %>%
  drop_na(
    dlog_capital,
    rsales_g_std, size_std, sh_current_a_std,
    lev_std_gdp, lev_std_wide, lev_std_mps,
    levL4_std_gdp, levL4_std_wide, levL4_std_mps
  )

m23_3 <- feols(
  dlog_capital ~ rsales_g_std + size_std + sh_current_a_std |
    name + Country_x |
    (lev_std_gdp + lev_std_wide + lev_std_mps) ~
    (levL4_std_gdp + levL4_std_wide + levL4_std_mps),
  data    = df23_3,
  cluster = ~name + Country_x
)

## Model 4: DD instrumentado con lag 1
df23_4 <- df %>%
  drop_na(
    dlog_capital,
    rsales_g_std, size_std, sh_current_a_std,
    dd_std_gdp, dd_std_wide, dd_std_mps,
    ddL1_std_gdp, ddL1_std_wide, ddL1_std_mps
  )

m23_4 <- feols(
  dlog_capital ~ rsales_g_std + size_std + sh_current_a_std |
    name + Country_x |
    (dd_std_gdp + dd_std_wide + dd_std_mps) ~
    (ddL1_std_gdp + ddL1_std_wide + ddL1_std_mps),
  data    = df23_4,
  cluster = ~name + Country_x
)

## Model 5: DD instrumentado con lag 2
df23_5 <- df %>%
  drop_na(
    dlog_capital,
    rsales_g_std, size_std, sh_current_a_std,
    dd_std_gdp, dd_std_wide, dd_std_mps,
    ddL2_std_gdp, ddL2_std_wide, ddL2_std_mps
  )

m23_5 <- feols(
  dlog_capital ~ rsales_g_std + size_std + sh_current_a_std |
    name + Country_x |
    (dd_std_gdp + dd_std_wide + dd_std_mps) ~
    (ddL2_std_gdp + ddL2_std_wide + ddL2_std_mps),
  data    = df23_5,
  cluster = ~name + Country_x
)

## Model 6: DD instrumentado con lag 4
df23_6 <- df %>%
  drop_na(
    dlog_capital,
    rsales_g_std, size_std, sh_current_a_std,
    dd_std_gdp, dd_std_wide, dd_std_mps,
    ddL4_std_gdp, ddL4_std_wide, ddL4_std_mps
  )

m23_6 <- feols(
  dlog_capital ~ rsales_g_std + size_std + sh_current_a_std |
    name + Country_x |
    (dd_std_gdp + dd_std_wide + dd_std_mps) ~
    (ddL4_std_gdp + ddL4_std_wide + ddL4_std_mps),
  data    = df23_6,
  cluster = ~name + Country_x
)

# 3) Mostrar Table 23
modelsummary(
  list(
    "Lev 1q" = m23_1,
    "Lev 2q" = m23_2,
    "Lev 4q" = m23_3,
    "DD 1q"  = m23_4,
    "DD 2q"  = m23_5,
    "DD 4q"  = m23_6
  ),
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c("lev_std_wide", "dd_std_wide"),
  gof_omit = ".*",
  title    = "Table 23: IV usando posiciones financieras rezagadas"
)





# --------------------------------------------------------
# Table 24: Decomposition of Leverage (sin embigl, con FE firma + país)
# --------------------------------------------------------

library(dplyr)
library(fixest)
library(modelsummary)

# 1) Crear y estandarizar medidas de posición financiera
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Leverage
    lev_std    = (leverage - mean(leverage,   na.rm = TRUE)) / sd(leverage,   na.rm = TRUE),
    # Net leverage = (deuda total – caja) / activos
    levnet     = (tdebt - cheq1) / atq,
    levnet_std = (levnet   - mean(levnet,   na.rm = TRUE)) / sd(levnet,   na.rm = TRUE),
    # ST debt share = deuda corriente / deuda total
    shstdt     = dlcq / tdebt,
    shstdt_std = (shstdt  - mean(shstdt,  na.rm = TRUE)) / sd(shstdt,  na.rm = TRUE),
    # LT debt share = deuda no corriente / deuda total
    shltdt     = dlttq / tdebt,
    shltdt_std = (shltdt  - mean(shltdt,  na.rm = TRUE)) / sd(shltdt,  na.rm = TRUE),
    # Other liabilities share = (activos – deuda total) / activos
    sh_ol      = (atq - tdebt) / atq,
    sh_ol_std  = (sh_ol   - mean(sh_ol,   na.rm = TRUE)) / sd(sh_ol,   na.rm = TRUE),
    # Total liabilities share = deuda total / activos = leverage
    sh_l       = leverage,
    sh_l_std   = lev_std
  ) %>%
  ungroup()

# 2) Construir interacciones con los shocks (usando mpso y mps)
df <- df %>%
  mutate(
    # Leverage
    lev_std_gdp               = lev_std    * dlog_gdp,
    lev_std_wide              = lev_std    * mpso,
    lev_std_mps               = lev_std    * mps,
    # Net leverage
    levnet_std_gdp            = levnet_std * dlog_gdp,
    levnet_std_wide           = levnet_std * mpso,
    levnet_std_mps            = levnet_std * mps,
    # ST debt share
    shstdt_std_gdp            = shstdt_std * dlog_gdp,
    shstdt_std_wide           = shstdt_std * mpso,
    shstdt_std_mps            = shstdt_std * mps,
    # LT debt share
    shltdt_std_gdp            = shltdt_std * dlog_gdp,
    shltdt_std_wide           = shltdt_std * mpso,
    shltdt_std_mps            = shltdt_std * mps,
    # Other liabilities share
    sh_ol_std_gdp             = sh_ol_std  * dlog_gdp,
    sh_ol_std_wide            = sh_ol_std  * mpso,
    sh_ol_std_mps             = sh_ol_std  * mps,
    # Total liabilities share
    sh_l_std_gdp              = sh_l_std   * dlog_gdp,
    sh_l_std_wide             = sh_l_std   * mpso,
    sh_l_std_mps              = sh_l_std   * mps
  )

# 3) Definir controles a nivel firma (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 4) Listas de variables para cada uno de los 7 modelos
#    EN M5: quitamos la colinealidad eliminando shltdt en la parte GDP,
#    ya que shstdt + shltdt + sh_ol = 1, lo que generaba colinealidad.
vars24 <- list(
  m1 = c("dlog_capital",
         "lev_std_gdp",   "lev_std_wide",   "lev_std_mps",
         controls_firm),
  m2 = c("dlog_capital",
         "levnet_std_gdp","levnet_std_wide","levnet_std_mps",
         controls_firm),
  m3 = c("dlog_capital",
         "shstdt_std_gdp","shstdt_std_wide","shstdt_std_mps",
         controls_firm),
  m4 = c("dlog_capital",
         "shltdt_std_gdp","shltdt_std_wide","shltdt_std_mps",
         controls_firm),
  m5 = c("dlog_capital",
         # Para evitar colinealidad, no incluimos shltdt_std_gdp junto a shstdt_std_gdp
         "shstdt_std_gdp",
         "shstdt_std_wide","shstdt_std_mps",
         "shltdt_std_wide","shltdt_std_mps",
         controls_firm),
  m6 = c("dlog_capital",
         "sh_ol_std_gdp","sh_ol_std_wide","sh_ol_std_mps",
         controls_firm),
  m7 = c("dlog_capital",
         "sh_l_std_gdp","sh_l_std_wide","sh_l_std_mps",
         controls_firm)
)

# 5) Estimar los 7 modelos con efectos fijos por firma y país, cluster por firma y país
models24 <- lapply(vars24, function(vs) {
  df_sub <- df %>% select(name, Country_x, all_of(vs)) %>% drop_na()
  fml    <- as.formula(paste(
    "dlog_capital ~",
    paste(setdiff(vs, "dlog_capital"), collapse = " + "),
    "| name + Country_x"
  ))
  feols(fml, data = df_sub, cluster = ~name + Country_x)
})

# 6) Mostrar Table 24 en consola
modelsummary(
  models24,
  stars    = c('*' = .10, '**' = .05, '***' = .01),
  keep     = c(
    "lev_std_wide",      "levnet_std_wide",
    "shstdt_std_wide",   "shltdt_std_wide",
    "sh_ol_std_wide",    "sh_l_std_wide"
  ),
  order    = c(
    "lev_std_wide",      "levnet_std_wide",
    "shstdt_std_wide",   "shltdt_std_wide",
    "sh_ol_std_wide",    "sh_l_std_wide"
  ),
  coef_map = c(
    "lev_std_wide"        = "Leverage × shock (wide)",
    "levnet_std_wide"     = "Net leverage × shock (wide)",
    "shstdt_std_wide"     = "ST debt share × shock (wide)",
    "shltdt_std_wide"     = "LT debt share × shock (wide)",
    "sh_ol_std_wide"      = "Other liabilities share × shock (wide)",
    "sh_l_std_wide"       = "Total liabilities share × shock (wide)"
  ),
  gof_omit = ".*",
  title    = "Table 24: Decomposition of Leverage"
)


