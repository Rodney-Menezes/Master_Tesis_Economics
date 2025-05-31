# Script en R para replicar las Tablas 1–2 y la Figura del paper de Ottonello & Winberry (2020)
# == Instalacion de Paquetes necesarios ==

install.packages("haven")
install.packages("dplyr")
install.packages("zoo")
install.packages("fixest")
install.packages("mFilter")
install.packages("xtable")
install.packages("openxlsx")
install.packages("tidyr")

# ===================================================
# Script en R para replicar “descriptive_statistics.do”
# Ottonello & Winberry (2020) adaptado a tu base
# ===================================================

# 0) Cargar librerías
library(haven)    # read_dta()
library(dplyr)    # manipulación
library(zoo)      # as.yearqtr()
library(fixest)   # feols(), resid()

# 1) Directorio base y carga de datos
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq = as.yearqtr(dateq)  # convierte tu variable trimestral
  )

# ---------------------------------------------------
# 2) Tabla 1: Estadísticas de shocks monetarios
# ---------------------------------------------------
# Renombra tus shocks (ya en %)
df1 <- df %>% rename(
  HF_shock = mps,   # Monetary Policy Surprise (raw)
  Q_shock  = mpso   # Otras sorpresas ortogonalizadas
)

# Calcula estadísticos
stats_HF <- df1 %>% summarise(
  Mean   = mean(HF_shock,   na.rm = TRUE),
  Median = median(HF_shock, na.rm = TRUE),
  SD     = sd(HF_shock,     na.rm = TRUE),
  Min    = min(HF_shock,    na.rm = TRUE),
  Max    = max(HF_shock,    na.rm = TRUE),
  N      = sum(!is.na(HF_shock))
)
stats_Q  <- df1 %>% summarise(
  Mean   = mean(Q_shock,   na.rm = TRUE),
  Median = median(Q_shock, na.rm = TRUE),
  SD     = sd(Q_shock,     na.rm = TRUE),
  Min    = min(Q_shock,    na.rm = TRUE),
  Max    = max(Q_shock,    na.rm = TRUE),
  N      = sum(!is.na(Q_shock))
)

# Imprime Tabla 1
tabla1 <- data.frame(
  Estadístico     = c("Mean","Median","S.D.","Min","Max","Observaciones"),
  High_Frequency  = unlist(stats_HF),
  Smoothed_Quarter = unlist(stats_Q)
)
cat("\n--- Tabla 1: Estadísticas de shocks monetarios ---\n")
print(tabla1)

# -------------------------------
# Tabla 2: Firm-Level Variables
# -------------------------------

# 1) Selección de variables disponibles
panel <- df %>% select(
  name,          # Identificador (ticker) :contentReference[oaicite:1]{index=1}
  dateq,         # Fecha trimestral :contentReference[oaicite:2]{index=2}
  dlog_capital,  # Δ log capital :contentReference[oaicite:3]{index=3}
  leverage,      # Apalancamiento = tdebt/atq :contentReference[oaicite:4]{index=4}
  dd             # Distance to default :contentReference[oaicite:5]{index=5}
)

# Panel A: distribuciones marginales de dlog_capital, leverage y dd
varsA <- c("dlog_capital","leverage","dd")
tabA <- panel %>%
  summarise(across(all_of(varsA),
                   list(
                     Mean   = ~mean(.x, na.rm=TRUE),
                     Median = ~median(.x, na.rm=TRUE),
                     SD     = ~sd(.x, na.rm=TRUE),
                     P95    = ~quantile(.x, .95, na.rm=TRUE),
                     N      = ~sum(!is.na(.x))
                   ),
                   .names = "{.col}_{.fn}"
  ))

# Reformatea para impresión
tabA_print <- do.call(rbind, lapply(varsA, function(v){
  data.frame(
    Variable      = v,
    Mean          = tabA[[paste0(v,"_Mean")]],
    Median        = tabA[[paste0(v,"_Median")]],
    `S.D.`        = tabA[[paste0(v,"_SD")]],
    `95th Pct.`   = tabA[[paste0(v,"_P95")]],
    Observaciones = tabA[[paste0(v,"_N")]]
  )
}))
cat("\n--- Tabla 2A: Distribuciones marginales ---\n")
print(tabA_print)

# Panel B: matriz de correlaciones (raw)
corrB <- panel %>%
  select(all_of(varsA)) %>%
  cor(use="pairwise.complete.obs")
cat("\n--- Tabla 2B: Correlaciones (raw) ---\n")
print(round(corrB, 3))

# -------------------------------
# Panel C: correlaciones residualizadas (solo FE por firma)
# -------------------------------

library(dplyr)
library(fixest)
library(tidyr)

# 1) Prepara los datos y elimina filas con NA
panel_fc <- df %>%
  select(
    name,
    leverage,
    dd,
    dln_saleq,
    log_size,
    current_ratio
  ) %>%
  drop_na(leverage, dd, dln_saleq, log_size, current_ratio)

# 2) Residualiza usando solo efectos fijos de firma
panel_resid <- panel_fc %>%
  mutate(
    leverage_resid = resid(
      feols(leverage ~ dln_saleq + log_size + current_ratio 
            | name, data = panel_fc)
    ),
    dd_resid = resid(
      feols(dd ~ dln_saleq + log_size + current_ratio 
            | name, data = panel_fc)
    )
  )

# 3) Calcula y muestra la matriz de correlaciones de los residuales
corrC <- panel_resid %>%
  select(leverage_resid, dd_resid) %>%
  cor(use = "pairwise.complete.obs")

cat("\n--- Tabla 2C: Correlaciones (residualizadas, FE firma) ---\n")
print(round(corrC, 3))

