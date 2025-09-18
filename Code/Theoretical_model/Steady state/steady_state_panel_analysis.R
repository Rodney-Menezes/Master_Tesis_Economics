# --------------------------------------
# AUTORES: Tú
# PROPÓSITO: Calcular momentos empíricos con tus datos panel
# FECHA: 10/06/2025
# --------------------------------------

# 1. Cargar librerías
library(readxl)    # Para leer Excel
library(dplyr)     # Manipulación de datos
library(haven)     # Leer/escribir archivos .dta

# 2. Leer y renombrar variables
df <- read_excel(
  "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.xlsx"
) %>%
  select(
    name,         # Identificador del activo/empresa
    dateq,        # Fecha trimestral
    real_inv,     # Inversión real = netinv / cpi
    cashop,       # Flujo operativo de caja
    firm_value,   # Valor de mercado de la firma
    real_capital, # Capital real = capital / cpi
    tdebt,        # Deuda total
    unemp         # Empleo
  ) %>%
  rename(
    firm_id      = name,
    quarter      = dateq,
    investment   = real_inv,
    cash_flow    = cashop,
    market_value = firm_value,
    capital      = real_capital,
    debt         = tdebt,
    employment   = unemp
  )

# 3. Ordenar y crear las variables de nivel dentro de cada firma
df <- df %>%
  arrange(firm_id, quarter) %>%
  group_by(firm_id) %>%
  mutate(
    i_k       = investment / capital,
    tobin_q   = market_value / capital,
    cf_k      = dplyr::lag(cash_flow, 1) / capital,
    lev       = ifelse(debt < 0,
                       debt / (capital - debt),
                       debt / capital),
    pos_debt  = as.integer(debt > 0),
    gross_lev = debt * pos_debt / capital,
    emp_dhs   = (employment - dplyr::lag(employment, 1)) /
      ((employment + dplyr::lag(employment, 1)) / 2)
  ) %>%
  ungroup()

# 4. En un segundo paso, crear los rezagos explícitos
df <- df %>%
  group_by(firm_id) %>%
  mutate(
    lev_lag       = dplyr::lag(lev,       1),
    gross_lev_lag = dplyr::lag(gross_lev, 1)
  ) %>%
  ungroup()

# 5. Guardar en .dta en la ruta especificada
haven::write_dta(
  df,
  "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/simulated_data.dta"
)

# 6. Filtrar infinitos y NAs antes de las correlaciones
df_corr <- df %>%
  filter(
    is.finite(lev),       is.finite(lev_lag),
    is.finite(gross_lev), is.finite(gross_lev_lag)
  )

cat("Observaciones totales:      ", nrow(df),     "\n")
cat("Pares útiles para corr():  ", nrow(df_corr), "\n\n")

# 7. Correlaciones condicionales
if (nrow(df_corr) >= 2) {
  # Corr. “filtrada” lev vs. lev_lag
  corr_lev       <- cor(df_corr$lev,       df_corr$lev_lag)
  # Corr. “filtrada” gross_lev vs. gross_lev_lag
  corr_gross     <- cor(df_corr$gross_lev, df_corr$gross_lev_lag)
  
  # Además, replicar directamente las dos correlaciones originales:
  corr_full    <- with(
    df %>% filter(is.finite(lev)),
    cor(lev, dplyr::lag(lev), use = "complete.obs")
  )
  corr_ik_lev  <- with(
    df %>% filter(is.finite(lev) & is.finite(i_k)),
    cor(i_k, lev, use = "complete.obs")
  )
  
  # Imprimir resultados
  cat("corr(lev, lev_lag) (filtrada):             ", round(corr_lev,   4), "\n")
  cat("corr(gross_lev, gross_lev_lag) (filtrada): ", round(corr_gross, 4), "\n\n")
  cat("corr_full (lev vs. lev rezagado):          ", round(corr_full,  4), "\n")
  cat("corr_ik_lev (i_k vs. lev):                 ", round(corr_ik_lev,4), "\n")
} else {
  cat("No hay suficientes pares completos para calcular correlaciones.\n")
}

