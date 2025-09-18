# --------------------------------------
# CALIBRATION TARGETS – CÁLCULO EN R (con log y nuevas categorías de edad)
# Fecha: 10/06/2025
# --------------------------------------

library(readxl)
library(dplyr)

# 1. Cargar solo las variables relevantes
df <- read_excel(
  "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.xlsx"
) %>%
  select(
    name,        # ID empresa
    dateq,       # periodo trimestral
    inv_rate,    # tasa de inversión real (proxy i/k trimestral)
    leverage,    # apalancamiento = tdebt/activos
    saleq,       # ventas trimestrales (proxy empleo)
    cashop,      # flujo operativo (proxy default)
    tdebt,       # deuda total
    age_sample   # edad de la firma en trimestres
  ) %>%
  rename(
    firm_id        = name,
    quarter_id     = dateq,
    i_rate         = inv_rate,
    gross_leverage = leverage,
    sales          = saleq,
    cash_flow      = cashop,
    debt           = tdebt,
    age_quarters   = age_sample
  ) %>%
  # 2. Construir variable "year" aproximada
  mutate(
    year = ceiling(age_quarters / 4)
  )

# 3. Agrupar por firma–año y calcular métricas anuales
df_ann <- df %>%
  group_by(firm_id, year) %>%
  summarise(
    # i/k anual: promedio trimestral * 4
    i_k_ann            = mean(i_rate, na.rm=TRUE) * 4,
    # Gross leverage anual
    gross_leverage_ann = mean(gross_leverage, na.rm=TRUE),
    # Default proxy (cualquier trimestre con flujo < 0 y deuda > 0)
    default_ann        = as.integer(any(cash_flow < 0 & debt > 0)),
    # Empleo proxy (suma de ventas en el año)
    employment_ann     = sum(sales, na.rm=TRUE),
    # Exit proxy (última aparición en el panel)
    exit_ann           = as.integer(max(age_quarters) == max(df$age_quarters[df$firm_id == firm_id])),
    .groups = "drop"
  )

# 4. Calcular edad en años dentro de cada firma
df_ann <- df_ann %>%
  group_by(firm_id) %>%
  mutate(
    first_year = min(year),
    age_year   = year - first_year + 1
  ) %>%
  ungroup()

# 5. Total agregado de “empleo” proxy
total_emp <- sum(df_ann$employment_ann, na.rm=TRUE)

# 6. Cálculo de los 10 Calibration targets corregidos
sd_log_i_k      <- sd(log(df_ann$i_k_ann), na.rm=TRUE)          # SD en log
mean_gross_lev  <- mean(df_ann$gross_leverage_ann, na.rm=TRUE)
default_rate    <- 100 * mean(df_ann$default_ann, na.rm=TRUE)
emp_share_age1  <- sum(df_ann$employment_ann[df_ann$age_year == 1],            na.rm=TRUE) / total_emp
emp_share_age2  <- sum(df_ann$employment_ann[df_ann$age_year == 2],            na.rm=TRUE) / total_emp
emp_share_age_gt2 <- sum(df_ann$employment_ann[df_ann$age_year >  2],          na.rm=TRUE) / total_emp
share_age1      <- mean(df_ann$age_year == 1, na.rm=TRUE)
share_age2      <- mean(df_ann$age_year == 2, na.rm=TRUE)
frac_pos_debt   <- mean(df_ann$gross_leverage_ann > 0, na.rm=TRUE)
exit_rate       <- 100 * mean(df_ann$exit_ann, na.rm=TRUE)

# 7. Tabla de resultados
calibration_targets <- tibble(
  Metric = c(
    "SD(log(i/k))",
    "Mean gross leverage ratio",
    "Annualized default rate (%)",
    "Employment share age1",
    "Employment share age2",
    "Employment share age>2",
    "Share firms age1",
    "Share firms age2",
    "Frac with positive debt",
    "Mean exit rate (%)"
  ),
  Data = c(
    sd_log_i_k,
    mean_gross_lev,
    default_rate,
    emp_share_age1,
    emp_share_age2,
    emp_share_age_gt2,
    share_age1,
    share_age2,
    frac_pos_debt,
    exit_rate
  )
)

print(calibration_targets)




##### PARTE 2.

# ------------------------------------------------
# CÁLCULO EN R – Untargetted statistics corregidas
# Fecha: 10/06/2025
# ------------------------------------------------

library(readxl)
library(dplyr)

# 1. Cargar variables relevantes
df <- read_excel(
  "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.xlsx"
) %>%
  select(
    Country,         # País de la empresa :contentReference[oaicite:0]{index=0}
    dateq,           # Trimestre (periodo) :contentReference[oaicite:1]{index=1}
    embigl,          # EMBI Global spread :contentReference[oaicite:2]{index=2}
    tdebt,           # Deuda total :contentReference[oaicite:3]{index=3}
    cheq1,           # Efectivo y equivalentes :contentReference[oaicite:4]{index=4}
    atq,             # Activos totales :contentReference[oaicite:5]{index=5}
    oiadpq,          # Ingreso operativo antes de D&A :contentReference[oaicite:6]{index=6}
    saleq,           # Ventas trimestrales :contentReference[oaicite:7]{index=7}
    current_ratio,   # Razón corriente :contentReference[oaicite:8]{index=8}
    cash_ratio       # Razón de caja :contentReference[oaicite:9]{index=9}
  ) %>%
  rename(
    country           = Country,
    quarter           = dateq,
    credit_spread_q   = embigl,
    debt              = tdebt,
    cash              = cheq1,
    assets            = atq,
    op_income         = oiadpq,
    sales             = saleq,
    curr_ratio        = current_ratio,
    cash_ratio        = cash_ratio
  )

# 2. Annualized credit spread:
#    — Primero extraemos unívocamente cada trimestre–país
cs <- df %>%
  distinct(country, quarter, credit_spread_q) %>%
  arrange(country, quarter)

#    — Luego sumamos todos los spreads y dividimos por el número total de observaciones
annualized_credit_spread <- sum(cs$credit_spread_q, na.rm = TRUE) /
  nrow(cs)

# 3. Mean net leverage ratio:
#    net_leverage = (debt – cash) / assets
mean_net_leverage <- df %>%
  mutate(net_leverage = (debt - cash) / assets) %>%
  summarise(mnl = mean(net_leverage, na.rm = TRUE)) %>%
  pull(mnl)

# 4. Aggregate TFP (proxy):
#    tfp_proxy = op_income / assets
mean_tfp_proxy <- df %>%
  mutate(tfp_proxy = op_income / assets) %>%
  summarise(mt = mean(tfp_proxy, na.rm = TRUE)) %>%
  pull(mt)

# 5. Estados crediticios (proxies de restricción):
stats_constr <- df %>%
  mutate(
    unconstrained_flag   = as.integer(curr_ratio >= 1),
    risky_constrained    = as.integer(curr_ratio < 1 & cash_ratio < 1),
    riskfree_constrained = as.integer(curr_ratio < 1 & cash_ratio >= 1)
  ) %>%
  summarise(
    frac_unconstrained    = mean(unconstrained_flag,   na.rm = TRUE),
    frac_risky_constr     = mean(risky_constrained,    na.rm = TRUE),
    frac_riskfree_constr  = mean(riskfree_constrained, na.rm = TRUE)
  )

# 6. Montar tabla final
untargetted_stats <- tibble(
  Metric = c(
    "Annualized credit spread",
    "Mean net leverage ratio",
    "Aggregate TFP (proxy)",
    "Fraction unconstrained",
    "Fraction Risky Constrained",
    "Fraction Risk-Free Constrained"
  ),
  Data = c(
    annualized_credit_spread,
    mean_net_leverage,
    mean_tfp_proxy,
    stats_constr$frac_unconstrained,
    stats_constr$frac_risky_constr,
    stats_constr$frac_riskfree_constr
  )
)

print(untargetted_stats)




#### Parte 3  : 

# --------------------------------------------------------
# CÁLCULO DE CALIBRACIÓN EN R
# Fecha: 10/06/2025
# Base: Data_Base_final.xlsx
# --------------------------------------------------------

library(readxl)
library(dplyr)

# 0. Leer datos y renombrar variables clave
df <- read_excel(
  "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.xlsx"
) %>%
  select(
    name,         # ticker
    dateq,        # trimestre
    tdebt,        # deuda total
    cheq1,        # efectivo
    atq,          # activos totales
    leverage,     # tdebt/atq
    capital,      # stock de capital
    netinv,       # inversión neta ΔPPE
    cpi,          # deflactor
    saleq,        # ventas trimestrales
    age_sample    # edad firma en trimestres
  ) %>%
  rename(
    firm_id      = name,
    quarter_id   = dateq,
    debt         = tdebt,
    cash         = cheq1,
    assets       = atq,
    gross_lev_q  = leverage,
    cap_q        = capital,
    netinv_q     = netinv,
    cpi          = cpi,
    sales_q      = saleq,
    age_q        = age_sample
  ) %>%
  mutate(
    # variables reales y anuales
    inv_real_q   = netinv_q / cpi,
    cap_real_q   = cap_q   / cpi,
    year         = ceiling(age_q / 4),
    in_compustat = quarter_id >= 28      # proxy: Q≥28
  )

# -----------------------------------------------------------------
# 3. Compustat sample
# -----------------------------------------------------------------
df_ann <- df %>%
  group_by(firm_id, year, in_compustat) %>%
  summarise(
    # sumas y promedios anuales
    sum_debt     = sum(debt,    na.rm=TRUE),
    sum_cash     = sum(cash,    na.rm=TRUE),
    sum_assets   = sum(assets,  na.rm=TRUE),
    sum_sales    = sum(sales_q, na.rm=TRUE),
    sum_cap      = sum(cap_real_q, na.rm=TRUE),
    # ratios
    gross_lev_ann= sum_debt / sum_assets,
    net_lev_ann  = (sum_debt - sum_cash) / sum_assets,
    emp_ann      = sum_sales,                    # proxy empleo = ventas
    age_year     = year - min(year) + 1,         # edad en años
    .groups="drop"
  )

comp <- df_ann %>%
  filter(in_compustat) 

compustat_stats <- tibble(
  Metric = c(
    "Frac with positive debt",
    "Mean gross leverage ratio",
    "Mean net leverage ratio",
    "Employment ratio (sales proxy)",
    "Age ratio"
  ),
  Data = c(
    mean(comp$gross_lev_ann > 0, na.rm=TRUE),
    mean(comp$gross_lev_ann,     na.rm=TRUE),
    mean(comp$net_lev_ann,       na.rm=TRUE),
    mean(comp$emp_ann) / mean(df_ann$emp_ann),        # empleo relativo vs muestra completa
    mean(comp$age_year, na.rm=TRUE) / mean(df_ann$age_year, na.rm=TRUE)
  )
)

print(compustat_stats)




# -----------------------------------------------------------------
# 4. Lifecycle dynamics (datos) – crecimiento de múltiples métricas por años de operación
# -----------------------------------------------------------------



# --------------------------------------
# 5. Summary statistics del panel real
# --------------------------------------
# --------------------------------------
# 5. Summary statistics del panel real (autocorrelación global)
# --------------------------------------

library(readxl)
library(dplyr)
library(lubridate)

# 1. Leer datos y renombrar variables reales
df <- read_excel(
  "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.xlsx"
) %>%
  select(
    firm_id      = name,         # Identificador de la empresa
    dateq,                        # Fecha trimestral
    real_inv,                     # Inversión real = netinv / cpi
    real_capital                  # Capital real = capital / cpi
  ) %>%
  mutate(
    dateq = as.Date(dateq),      # Asegurar clase Date
    year  = year(dateq)          # Extraer año calendario
  )

# 5.1 Calcular i/k anual (inversión real sobre capital real)
ik_ann <- df %>%
  group_by(firm_id, year) %>%
  summarise(
    sum_inv = sum(real_inv,        na.rm = TRUE),
    sum_cap = sum(real_capital,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(i_k = sum_inv / sum_cap) %>%
  arrange(firm_id, year)

# 5.2 Media y desviación estándar de i/k
mean_ik <- mean(ik_ann$i_k, na.rm = TRUE)
sd_ik   <- sd(ik_ann$i_k,   na.rm = TRUE)

# 5.3 Autocorrelación global de i/k (sin agrupar)
autocorr_ik <- cor(
  ik_ann$i_k,
  dplyr::lag(ik_ann$i_k),
  use = "complete.obs"
)

# 5.4 Construir tibble final
panel_stats <- tibble(
  Metric = c("E(i/k)", "SD(i/k)", "autocorr(i/k)"),
  Data   = c(mean_ik,   sd_ik,    autocorr_ik)
)

print(panel_stats)






