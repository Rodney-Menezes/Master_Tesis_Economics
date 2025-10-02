# ==========================
# Script R corregido: Estimar “shock monetario” por país
# (Regresión panel con Efectos Fijos y reescalado fuerte del residuo)
# Corrige el problema de fechas Excel que se leen como números.
# Entrada: Data1.xlsx (panel trimestral con Periodo, gdp, ipc, Country_x, interest_rate)
# Salida: shocks_monetarios_por_pais.xlsx (residuos reescalados por país y periodo)
# ==========================

# Instalación de mFilter si hace falta
if (!requireNamespace("mFilter", quietly = TRUE)) {
  install.packages("mFilter")
}

# 1) Cargar librerías necesarias
library(dplyr)
library(openxlsx)   # Para leer/escribir Excel y dar formato
library(plm)        # Para regresión en panel con efectos fijos
library(lmtest)     # Para coeftest y pruebas de significancia
library(sandwich)   # Para vcovHC en errores estándar robustos
library(ggplot2)    # (opcional) Para graficar los shocks
library(mFilter)    # Para el filtro HP

# 2) Definir rutas de entrada y salida
input_path  <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Taylor Model/Data_Base_Taylor_Model.xlsx"
output_path <- "C:/Users/joser/Desktop/shocks_monetarios_por_pais.xlsx"

# Factor de escala para los residuos
scale_factor <- 1e14

# 3) Leer datos
data_raw <- read.xlsx(input_path, sheet = 1)

# 4) Convertir Periodo a Date
if (is.numeric(data_raw$Periodo)) {
  data_raw$Periodo <- as.Date(data_raw$Periodo, origin = "1899-12-30")
} else {
  data_raw$Periodo <- as.Date(data_raw$Periodo)
}
if (any(is.na(data_raw$Periodo))) {
  warning("Algunas fechas en 'Periodo' no pudieron convertirse correctamente.")
}

# 5) Preprocesar panel: factor país y orden
panel_data <- data_raw %>%
  mutate(country = as.factor(Country_x)) %>%
  arrange(country, Periodo)

# 6) Calcular rezago de tasa e IPC
panel_data <- panel_data %>%
  group_by(country) %>%
  arrange(Periodo) %>%
  mutate(
    lag_interest = lag(interest_rate),
    ipc          = ipc
  ) %>%
  ungroup()

# 7) Crear el output gap por país con filtro HP
panel_data <- panel_data %>%
  group_by(country) %>%
  arrange(Periodo) %>%
  mutate(
    log_gdp    = log(gdp),
    hp_trend   = hpfilter(log_gdp, freq = 1600)$trend,
    output_gap = 100 * (log_gdp - hp_trend)   # en porcentaje
  ) %>%
  ungroup()

# 8) Filtrar observaciones completas (sin NAs en variables clave)
panel_data <- panel_data %>%
  filter(
    !is.na(lag_interest),
    !is.na(ipc),
    !is.na(output_gap),
    !is.na(interest_rate)
  )

# 9) Convertir a pdata.frame para plm
panel_data <- pdata.frame(panel_data, index = c("country", "Periodo"))

# 10) Estimar el modelo “within” con efectos fijos de país y periodo
model_hp <- plm(
  interest_rate ~ lag_interest + ipc + output_gap,
  data   = panel_data,
  model  = "within",
  effect = "twoways"    # efectos fijos país + tiempo
)

# 11) Mostrar resultados con errores robustos HC1 cluster por país
cat("==== Regresión Taylor con output gap (efectos fijos país + tiempo) ====\n")
coefs_hp <- coeftest(
  model_hp,
  vcov = vcovHC(model_hp, type = "HC1", cluster = "group")
)
print(coefs_hp)

# 12) Extraer y escalar residuos (proxy de shock monetario)
panel_data$shock_monetario <- residuals(model_hp)
panel_data$shock_scaled    <- scale_factor * panel_data$shock_monetario

# 13) Preparar la tabla final
shocks_panel <- as.data.frame(panel_data) %>%
  select(country, Periodo, shock_scaled)

# 14) Exportar a Excel con estilo de fecha en la columna Periodo
wb <- createWorkbook()
addWorksheet(wb, "Shocks")
writeData(wb, "Shocks", shocks_panel, withFilter = FALSE)

# Definir estilo de fecha yyyy-mm-dd
dateStyle <- createStyle(numFmt = "yyyy-mm-dd")

# Aplicar estilo a la columna Periodo (columna 2, filas 2 hasta n+1)
addStyle(
  wb, sheet = "Shocks", style = dateStyle,
  cols = 2,
  rows = 2:(nrow(shocks_panel) + 1),
  gridExpand = TRUE
)

# Guardar el archivo
saveWorkbook(wb, output_path, overwrite = TRUE)
cat("\n> Shocks monetarios guardados en:\n", output_path, "\n")

