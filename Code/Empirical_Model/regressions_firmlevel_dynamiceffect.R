# --------------------------------------------------------
# Dynamic Effects (horizonte 0–12 trimestres), usando mps
# --------------------------------------------------------

# Instalar paquetes necesarios (si aún no lo están)
install.packages("kableExtra")

library(haven)      # read_dta()
library(zoo)        # as.yearqtr()
library(purrr)      # map(), map_dfc()
library(fixest)     # feols()
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# ------------------------------------------------------------------
# Figure 1: Dynamics of Differential Response to Monetary Shocks
# ------------------------------------------------------------------

# 1) Leer datos y convertir dateq a trimestral
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(dateq = as.yearqtr(dateq))

# 2) Generar controles firm-level si no existen
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Crecimiento de ventas y estandarizado
    rsales_g     = log(saleq) - log(lag(saleq)),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    # Tamaño
    size_std     = (log(atq) - mean(log(atq), na.rm = TRUE)) / sd(log(atq), na.rm = TRUE),
    # Liquidez
    sh_current_a_std = (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE)
  ) %>%
  ungroup()

# 3) Winsorizar, de-medial y estandarizar leverage y dd, luego crear interacciones con mps
winsorize <- function(x, p = 0.005) {
  lo <- quantile(x, p, na.rm = TRUE)
  hi <- quantile(x, 1 - p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

df <- df %>%
  mutate(
    leverage_win = winsorize(leverage),
    dd_win       = winsorize(dd)
  ) %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    lev_dev = leverage_win - mean(leverage_win, na.rm = TRUE),
    dd_dev  = dd_win       - mean(dd_win,       na.rm = TRUE),
    lev_std = lev_dev / sd(lev_dev, na.rm = TRUE),
    dd_std  = dd_dev  / sd(dd_dev,  na.rm = TRUE),
    # Interacciones con mps (raw shock)
    lev_wins_dem_std = lev_std * mps,
    d2d_wins_dem_std = dd_std  * mps
  ) %>%
  ungroup()

# 4) Construir cumFh_dlog_capital para h = 0...12
df_dyn <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  do({
    tmp <- .
    # Para cada horizonte h de 0 a 12, acumular sumas de dlog_capital
    for (h in 0:12) {
      tmp[[paste0("cumF", h, "_dlog_capital")]] <-
        rowSums(
          map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)),
          na.rm = TRUE
        )
    }
    tmp
  }) %>%
  ungroup()

# 5) Definir controles firm-level para la regresión
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 6) Estimar dinámicas para h = 0…12

# Leverage
res_lev <- map(0:12, function(h) {
  fml <- as.formula(paste0(
    "cumF", h, "_dlog_capital ~ ",
    "lev_wins_dem_std + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~name + Country_x)
})

# Distance-to-default
res_dd <- map(0:12, function(h) {
  fml <- as.formula(paste0(
    "cumF", h, "_dlog_capital ~ ",
    "d2d_wins_dem_std + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~name + Country_x)
})

# 7) Extraer coeficientes y errores
lev_coefs <- tibble(
  horizon = 0:12,
  beta     = map_dbl(res_lev, ~ coef(.x)["lev_wins_dem_std"]),
  se       = map_dbl(res_lev, ~ sqrt(vcov(.x)["lev_wins_dem_std", "lev_wins_dem_std"]))
)
dd_coefs <- tibble(
  horizon = 0:12,
  beta     = map_dbl(res_dd, ~ coef(.x)["d2d_wins_dem_std"]),
  se       = map_dbl(res_dd, ~ sqrt(vcov(.x)["d2d_wins_dem_std", "d2d_wins_dem_std"]))
)

# 8) Guardar resultados (opcional) y mostrar
write.csv(lev_coefs,
          file.path(base_dir, "dynamics_lev_0_12.csv"),
          row.names = FALSE)
write.csv(dd_coefs,
          file.path(base_dir, "dynamics_dd_0_12.csv"),
          row.names = FALSE)

print(lev_coefs)
print(dd_coefs)

# 9) Preparar tabla combinada para presentación
table_dyn <- lev_coefs %>%
  rename(beta_lev = beta, se_lev = se) %>%
  left_join(
    dd_coefs %>% rename(beta_dd = beta, se_dd = se),
    by = "horizon"
  ) %>%
  mutate(
    Leverage = sprintf("%.3f (%.3f)", beta_lev, se_lev),
    `Dist-to-Default` = sprintf("%.3f (%.3f)", beta_dd, se_dd)
  ) %>%
  select(Horizon = horizon, Leverage, `Dist-to-Default`)

# 10) Imprimir con kableExtra
table_dyn %>%
  kable(
    caption = "Efecto dinámico: coeficiente (error estándar) para h = 0–12 trimestres",
    align   = c("c", "c", "c")
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))




# --------------------------------------------------------
# Figure 1: Dinámica de Respuesta - Gráficos en R (horizonte 0–12)
# --------------------------------------------------------

library(ggplot2)
library(readr)
library(gridExtra)
library(grid)      # para grid.newpage() y grid.draw()

# 1) Leer resultados para horizonte 0–12 (suprimir avisos de spec)
base_dir  <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
lev_coefs <- read_csv(file.path(base_dir, "dynamics_lev_0_12.csv"),     show_col_types = FALSE)
dd_coefs  <- read_csv(file.path(base_dir, "dynamics_dd_0_12.csv"),      show_col_types = FALSE)

# 2) Gráfico para Leverage × Shock
p_lev <- ggplot(lev_coefs, aes(x = horizon, y = beta)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
              alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    x        = "Horizonte (trimestres)",
    y        = "Efecto acumulado de inversión",
    title    = "Dinámica de Respuesta: Leverage × Shock",
    subtitle = "Coeficiente ± 1.96 × E.E."
  ) +
  theme_minimal()

# 3) Gráfico para Distance-to-Default × Shock
p_dd <- ggplot(dd_coefs, aes(x = horizon, y = beta)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
              alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    x        = "Horizonte (trimestres)",
    y        = "Efecto acumulado de inversión",
    title    = "Dinámica de Respuesta: Distance-to-Default × Shock",
    subtitle = "Coeficiente ± 1.96 × E.E."
  ) +
  theme_minimal()

# 4) Mostrar ambos gráficos uno debajo del otro

# Opción 1: con gridExtra sin print()
grid.arrange(p_lev, p_dd, ncol = 1)

# Crear nueva “página” gráfica
grid.newpage()
# Dibujar ambos plots juntos
g <- grid.arrange(p_lev, p_dd, ncol = 1)
grid.draw(g)


# Opción 2: si prefieres patchwork, descomenta esto:
# install.packages("patchwork")
library(patchwork)
p_lev / p_dd




# --------------------------------------------------------
# Figure 6a: Respuestas dinámicas de los pagos de interés (sin embigl, horizonte 0–12)
# --------------------------------------------------------

# Instalar patchwork si no está instalado
install.packages("patchwork")

library(haven)       # read_dta()
library(zoo)         # as.yearqtr()
library(dplyr)
library(purrr)
library(fixest)      # feols()
library(knitr)
library(kableExtra)
library(ggplot2)
library(patchwork)   # para apilar plots

# 1) Leer datos y convertir dateq a trimestral
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(dateq = as.yearqtr(dateq)) %>%
  arrange(name, dateq)

# 2) Crear/renombrar variable de pagos de interés
df <- df %>%
  mutate(
    int_exp = oiadpq  # Ajusta si tu variable de pagos de interés tiene otro nombre
  )

# 3) Generar controles firm-level y choques
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    rsales_g_std          = { x <- log(saleq) - log(lag(saleq)); (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    size_std              = { x <- log(atq); (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    sh_current_a_std      = { x <- current_ratio; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    lev_std               = { x <- leverage; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    dd_std                = { x <- dd;       (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Interacción con shock raw (mps)
    lev_wins_dem_std_wide = lev_std * mps,
    d2d_wins_dem_std_wide = dd_std  * mps
  ) %>%
  ungroup()

# 4) Construir variables dinámicas cumFh_int_exp para h = 0…12
df_dyn <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    cumF0_int_exp  = int_exp,
    cumF1_int_exp  = int_exp + lead(int_exp, 1),
    cumF2_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2),
    cumF3_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3),
    cumF4_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4),
    cumF5_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5),
    cumF6_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6),
    cumF7_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6) + lead(int_exp, 7),
    cumF8_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6) + lead(int_exp, 7) + lead(int_exp, 8),
    cumF9_int_exp  = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6) + lead(int_exp, 7) + lead(int_exp, 8) + lead(int_exp, 9),
    cumF10_int_exp = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6) + lead(int_exp, 7) + lead(int_exp, 8) + lead(int_exp, 9) + lead(int_exp, 10),
    cumF11_int_exp = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6) + lead(int_exp, 7) + lead(int_exp, 8) + lead(int_exp, 9) + lead(int_exp, 10) + lead(int_exp, 11),
    cumF12_int_exp = int_exp + lead(int_exp, 1) + lead(int_exp, 2) + lead(int_exp, 3) + lead(int_exp, 4) + lead(int_exp, 5) + lead(int_exp, 6) + lead(int_exp, 7) + lead(int_exp, 8) + lead(int_exp, 9) + lead(int_exp, 10) + lead(int_exp, 11) + lead(int_exp, 12)
  ) %>%
  ungroup()

# 5) Definir controles para la regresión (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 6) Estimar dinámicas para h = 0…12 (solo coeficiente × choque wide)
res_IR_lev <- map(0:12, function(h) {
  form <- as.formula(paste0(
    "cumF", h, "_int_exp ~ lev_wins_dem_std_wide + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(form, data = df_dyn, cluster = ~name + Country_x)
})

res_IR_dd <- map(0:12, function(h) {
  form <- as.formula(paste0(
    "cumF", h, "_int_exp ~ d2d_wins_dem_std_wide + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(form, data = df_dyn, cluster = ~name + Country_x)
})

# 7) Extraer coeficientes y errores
IR_lev <- tibble(
  horizon = 0:12,
  beta_lev = map_dbl(res_IR_lev, ~ coef(.x)["lev_wins_dem_std_wide"]),
  se_lev   = map_dbl(res_IR_lev, ~ sqrt(vcov(.x)["lev_wins_dem_std_wide", "lev_wins_dem_std_wide"]))
)

IR_dd <- tibble(
  horizon = 0:12,
  beta_dd  = map_dbl(res_IR_dd, ~ coef(.x)["d2d_wins_dem_std_wide"]),
  se_dd    = map_dbl(res_IR_dd, ~ sqrt(vcov(.x)["d2d_wins_dem_std_wide", "d2d_wins_dem_std_wide"]))
)

# 8) Combinar y presentar en tabla con kableExtra
table_IR <- IR_lev %>%
  left_join(IR_dd, by = "horizon") %>%
  mutate(
    Leverage          = sprintf("%.3f (%.3f)", beta_lev, se_lev),
    `Dist-to-Default` = sprintf("%.3f (%.3f)", beta_dd, se_dd)
  ) %>%
  select(Horizon = horizon, Leverage, `Dist-to-Default`)

table_IR %>%
  kable(
    caption = "Figure 6a: Respuestas dinámicas de pagos de interés (horizonte 0–12)",
    align   = c("c", "c", "c")
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))

# 9) Crear gráficos con ggplot2 y apilarlos usando patchwork

p_lev <- ggplot(IR_lev, aes(x = horizon, y = beta_lev)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta_lev - 1.96 * se_lev,
                  ymax = beta_lev + 1.96 * se_lev),
              alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Dinámica IR: Leverage × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Δ pagos de interés acumulados"
  ) +
  theme_minimal()

p_dd <- ggplot(IR_dd, aes(x = horizon, y = beta_dd)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta_dd - 1.96 * se_dd,
                  ymax = beta_dd + 1.96 * se_dd),
              alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Dinámica IR: Distance-to-Default × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Δ pagos de interés acumulados"
  ) +
  theme_minimal()

# 10) Apilar y mostrar los dos gráficos con patchwork
(p_lev / p_dd) +
  plot_annotation(
    title = "Figure 6a: Respuestas dinámicas de pagos de interés (horizonte 0–12)"
  )







# --------------------------------------------------------
# Figure X: Dinámica de Respuesta de Inversión (Δ log k) ante Choques Monetarios
# --------------------------------------------------------

# --------------------------------------------------------------------------------
# 0) Instalar paquetes necesarios (si no están)
# --------------------------------------------------------------------------------
install.packages("patchwork")   # Para apilar los gráficos

# --------------------------------------------------------------------------------
# 1) Cargar librerías
# --------------------------------------------------------------------------------
library(haven)       # read_dta()
library(zoo)         # as.yearqtr()
library(dplyr)       # manipulaciones
library(purrr)       # map()
library(fixest)      # feols()
library(ggplot2)     # para graficar
library(patchwork)   # para apilar plots

# --------------------------------------------------------------------------------
# 2) Leer datos y preparar variables básicas
# --------------------------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq = as.yearqtr(dateq)
  ) %>%
  arrange(name, dateq)

# Calcular Δ log capital (Δ log k) si no existe:
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    dlog_capital = log(capital) - log(lag(capital))  # Ajusta “capital” al nombre real de tu variable de capital
  ) %>%
  ungroup()

# --------------------------------------------------------------------------------
# 3) Generar controles firm-level y choques estandarizados
# --------------------------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Crecimiento de ventas y estandarizado
    rsales_g_std      = { x <- log(saleq) - log(lag(saleq)); (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Tamaño
    size_std          = { x <- log(atq); (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Liquidez
    sh_current_a_std  = { x <- current_ratio; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Leverage estándar
    lev_std           = { x <- leverage; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # DD estándar
    dd_std            = { x <- dd;       (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Choques raw (mps)
    lev_wins_dem_std_wide = lev_std * mps,
    d2d_wins_dem_std_wide = dd_std  * mps
  ) %>%
  ungroup()

# --------------------------------------------------------------------------------
# 4) Construir variables dinámicas cumFh para Δ log capital, h = 0…12
# --------------------------------------------------------------------------------
df_dyn <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  do({
    tmp <- .
    for (h in 0:12) {
      # suma acumulada de dlog_capital desde t hasta t+h
      tmp[[paste0("cumF", h, "_dlog_capital")]] <-
        rowSums(
          map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)),
          na.rm = TRUE
        )
    }
    tmp
  }) %>%
  ungroup()

# --------------------------------------------------------------------------------
# 5) Definir controles para las regresiones (sin embigl)
# --------------------------------------------------------------------------------
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# --------------------------------------------------------------------------------
# 6) Estimar regresiones dinámicas para h = 0…12
#    (a) Leverage × Shock
#    (b) Distance‐to‐Default × Shock
# --------------------------------------------------------------------------------
res_lev <- map(0:12, function(h) {
  form <- as.formula(paste0(
    "cumF", h, "_dlog_capital ~ lev_wins_dem_std_wide + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(form, data = df_dyn, cluster = ~name + Country_x)
})

res_dd <- map(0:12, function(h) {
  form <- as.formula(paste0(
    "cumF", h, "_dlog_capital ~ d2d_wins_dem_std_wide + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(form, data = df_dyn, cluster = ~name + Country_x)
})

# --------------------------------------------------------------------------------
# 7) Extraer coeficientes y errores
# --------------------------------------------------------------------------------
lev_coefs <- tibble(
  horizon = 0:12,
  beta     = map_dbl(res_lev, ~ coef(.x)["lev_wins_dem_std_wide"]),
  se       = map_dbl(res_lev, ~ sqrt(vcov(.x)["lev_wins_dem_std_wide", "lev_wins_dem_std_wide"]))
)

dd_coefs <- tibble(
  horizon = 0:12,
  beta     = map_dbl(res_dd, ~ coef(.x)["d2d_wins_dem_std_wide"]),
  se       = map_dbl(res_dd, ~ sqrt(vcov(.x)["d2d_wins_dem_std_wide", "d2d_wins_dem_std_wide"]))
)

# (Opcional) Guardar resultados
write.csv(lev_coefs, file.path(base_dir, "dynamics_invcap_lev_0_12.csv"), row.names = FALSE)
write.csv(dd_coefs,  file.path(base_dir, "dynamics_invcap_dd_0_12.csv"),  row.names = FALSE)

# --------------------------------------------------------------------------------
# 8) Crear gráficos con ggplot2 y apilarlos con patchwork
# --------------------------------------------------------------------------------
p_lev_invcap <- ggplot(lev_coefs, aes(x = horizon, y = beta)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
              alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Dinámica IR: Leverage × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Δ inversión acumulada (Δ log k)"
  ) +
  theme_minimal()

p_dd_invcap <- ggplot(dd_coefs, aes(x = horizon, y = beta)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
              alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Dinámica IR: DD × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Δ inversión acumulada (Δ log k)"
  ) +
  theme_minimal()

# 9) Apilar y mostrar ambos gráficos con patchwork
(p_lev_invcap / p_dd_invcap) +
  plot_annotation(
    title = "Figure X: Respuesta Dinámica de la Inversión (Δ log k) ante Choques Monetarios (horizonte 0–12)"
  )





# --------------------------------------------------------
# Figure 9: Dynamics, Not Controlling for Differences in Cyclical Sensitivities (horizonte 0–12)
# --------------------------------------------------------

# Instalar patchwork si no está instalado
install.packages("patchwork")

library(haven)       # read_dta()
library(zoo)         # as.yearqtr()
library(dplyr)       # manipulaciones
library(purrr)       # map()
library(fixest)      # feols()
library(ggplot2)     # para graficar
library(patchwork)   # apilar plots

# 1) Leer panel y convertir dateq a trimestral
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(dateq = as.yearqtr(dateq)) %>%
  arrange(name, dateq)

# 2) Asegurar que existen las interacciones y controles firm-level
#    (Si ya fueron creadas antes, puedes omitir esta sección)
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # Controles firm-level
    rsales_g_std          = { x <- log(saleq) - log(lag(saleq)); (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    size_std              = { x <- log(atq); (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    sh_current_a_std      = { x <- current_ratio; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Estandarizar leverage y dd
    lev_std               = { x <- leverage; (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    dd_std                = { x <- dd;       (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) },
    # Interacciones con choque raw (mps)
    lev_wins_dem_std_wide = lev_std * mps,
    lev_wins_dem_std      = lev_std * mps,   # si quieres incluir raw
    d2d_wins_dem_std_wide = dd_std  * mps,
    d2d_wins_dem_std      = dd_std  * mps
  ) %>%
  ungroup()

# 3) Construir cumFh_dlog_capital para h = 0…12
df_dyn <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  do({
    tmp <- .
    for (h in 0:12) {
      tmp[[paste0("cumF", h, "_dlog_capital")]] <-
        rowSums(
          map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)),
          na.rm = TRUE
        )
    }
    tmp
  }) %>%
  ungroup()

# 4) Definir controles de firma (sin embigl)
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")

# 5) Estimar dinámicas h = 0…12

# (a) Leverage × Shock (no incluye componante ciclo)
res_lev_nocy <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ ",
      "lev_wins_dem_std_wide + lev_wins_dem_std + ",
      paste(controls_firm, collapse = " + "),
      " | name + Country_x"
    )),
    data    = df_dyn,
    cluster = ~name + Country_x
  )
})

# (b) Distance‐to‐Default × Shock (no incluye componante ciclo)
res_dd_nocy <- map(0:12, function(h) {
  feols(
    as.formula(paste0(
      "cumF", h, "_dlog_capital ~ ",
      "d2d_wins_dem_std_wide + d2d_wins_dem_std + ",
      paste(controls_firm, collapse = " + "),
      " | name + Country_x"
    )),
    data    = df_dyn,
    cluster = ~name + Country_x
  )
})

# 6) Extraer coeficientes y errores
dyn_lev_nocy <- tibble(
  horizon = 0:12,
  beta     = map_dbl(res_lev_nocy, ~ coef(.x)["lev_wins_dem_std_wide"]),
  se       = map_dbl(res_lev_nocy, ~ sqrt(vcov(.x)["lev_wins_dem_std_wide", "lev_wins_dem_std_wide"]))
)

dyn_dd_nocy <- tibble(
  horizon = 0:12,
  beta     = map_dbl(res_dd_nocy, ~ coef(.x)["d2d_wins_dem_std_wide"]),
  se       = map_dbl(res_dd_nocy, ~ sqrt(vcov(.x)["d2d_wins_dem_std_wide", "d2d_wins_dem_std_wide"]))
)

# 7) Guardar resultados (opcional)
write.csv(dyn_lev_nocy, file.path(base_dir, "dynamics_lev_nocyclicalsensitivities.csv"), row.names = FALSE)
write.csv(dyn_dd_nocy,  file.path(base_dir, "dynamics_dd_nocyclicalsensitivities.csv"),  row.names = FALSE)

# 8) Graficar (horizonte 0–12)

p_nocy_lev <- ggplot(dyn_lev_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figure 9(a): Dynamics (No Cyclical Sens.) – Leverage × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

p_nocy_dd <- ggplot(dyn_dd_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), alpha = 0.2, fill = "grey70") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figure 9(b): Dynamics (No Cyclical Sens.) – DD × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

# 9) Mostrar con patchwork
p_nocy_lev / p_nocy_dd +
  plot_annotation(
    title = "Figure 9: Dynamics (No Cyclical Sensitivities) (horizonte 0–12)"
  )





library(haven); library(zoo); library(dplyr); library(purrr)
library(fixest)

# Ruta del .dta
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),               # convierte a yearqtr
    quarter = as.factor(dateq)                  # factor trimestre
  ) %>%
  arrange(gvkey, dateq)

# 1) Generar dlog_capital
df <- df %>%
  group_by(gvkey) %>%
  mutate(
    dlog_capital = log(capital) - log(lag(capital))
  ) %>%
  ungroup()

# 2) Generar balances firm‐level y estandarizarlos (artefactos similares a tu script)
df <- df %>%
  group_by(gvkey) %>%
  arrange(dateq) %>%
  mutate(
    # Crecimiento de ventas:
    rsales_g     = log(saleq) - log(lag(saleq)),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    # Tamaño:
    size_std     = (log(atq) - mean(log(atq), na.rm = TRUE)) / sd(log(atq), na.rm = TRUE),
    # Liquidez:
    sh_current_a_std = (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE),
    # Leverage estandarizado y demeaned:
    lev_std = (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
    # Distance‐to‐default estandarizado y demeaned:
    dd_std  = (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE)
  ) %>%
  ungroup()

# 3) Suponiendo que ya tienes mps (raw shock) y mpso (ortogonalizado),
#    los renombramos para facilitar (wide=mps_ortho)
df <- df %>%
  rename(
    wide = mpso    # "wide" = shock ortogonalizado (mpso)
    # y asumimos “mps” existe ya en df
  )

# 4) Crear las interacciones: 
df <- df %>%
  mutate(
    # Con GDP (horizonte estático, aunque no se usarán aquí)
    lev_wins_dem_std_gdp  = lev_std * dlog_gdp,
    d2d_wins_dem_std_gdp  = dd_std  * dlog_gdp,
    # Con shocks (wide=mpso, mps=raw)
    lev_wins_dem_std_wide = lev_std * wide,
    lev_wins_dem_std      = lev_std * mps,
    d2d_wins_dem_std_wide = dd_std  * wide,
    d2d_wins_dem_std      = dd_std  * mps
  )

# 5) Definir controles firm‐level y agregados (igual que en tu script original)
controls_firm     <- c("rsales_g_std", "size_std", "sh_current_a_std")
controls_aggregate <- c(
  paste0("L", 1:4, "_dlog_gdp"),
  paste0("L", 1:4, "_dlog_cpi"),
  paste0("L", 1:4, "_unemp")
)

# 6) Generar variables dinámicas acumuladas cumF0...cumF20 de dlog_capital
df_dyn <- df %>%
  group_by(gvkey) %>%
  arrange(dateq) %>%
  do({
    tmp <- .
    for (h in 0:20) {
      tmp[[paste0("cumF", h, "_dlog_capital")]] <-
        rowSums(
          map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)),
          na.rm = TRUE
        )
    }
    tmp
  }) %>%
  ungroup()



# ------------------------------------------------------------
# Figura 10: Efecto promedio de inversión ante choque monetario
# ------------------------------------------------------------


# ------------------------------------------------------------
# Figura 10 (ajustada): Respuesta dinámica de inversión (0–12) usando mps
# ------------------------------------------------------------

# 0) Cargar librerías necesarias
# Si hace falta, instala antes con:
# install.packages(c("haven","zoo","dplyr","purrr","fixest","ggplot2"))
library(haven)     # Para leer archivos .dta
library(zoo)       # Para as.yearqtr()
library(dplyr)     # Para manipulaciones de datos
library(purrr)     # Para map() y map_dfc()
library(fixest)    # Para feols()
library(ggplot2)   # Para graficar

# 1) Leer panel desde .dta y crear variable de trimestre
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),   # convierte a objeto yearqtr
    quarter = as.factor(dateq)     # factor para trimestre (si hace falta más adelante)
  ) %>%
  arrange(name, dateq)             # orden por firma ("name") y fecha

# 2) Preprocesar variables para evitar log(0) / log(NA)
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) Variables auxiliares que ponen en NA cualquier valor ≤ 0
    capital2 = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    saleq2   = if_else(!is.na(saleq)   & saleq   > 0, saleq,   NA_real_),
    atq2     = if_else(!is.na(atq)     & atq     > 0, atq,     NA_real_),
    
    # 2.2) Calcular logaritmos seguros (log(NA)→NA sin warning)
    log_capital   = log(capital2),
    lag_log_capital = lag(log_capital),
    log_saleq     = log(saleq2),
    lag_log_saleq   = lag(log_saleq),
    log_atq       = log(atq2),
    
    # 2.3) Calcular diferencias de logaritmos (dlog_capital y rsales_g)
    dlog_capital = log_capital - lag_log_capital,
    rsales_g     = log_saleq   - lag_log_saleq
  ) %>%
  # 2.4) Estandarizar rsales_g
  mutate(
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE)
  ) %>%
  # 2.5) Estandarizar tamaño (size) usando log_atq
  mutate(
    size_std = if_else(
      !is.na(log_atq),
      (log_atq - mean(log_atq, na.rm = TRUE)) / sd(log_atq, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  # 2.6) Estandarizar liquidez (current_ratio)
  mutate(
    sh_current_a_std = if_else(
      !is.na(current_ratio),
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  # 2.7) Estandarizar leverage (lev) y distance‐to‐default (dd)
  mutate(
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  # 2.8) Eliminar variables auxiliares ya no necesarias
  select(-capital2, -saleq2, -atq2,
         -log_capital, -lag_log_capital,
         -log_saleq, -lag_log_saleq,
         -log_atq)

# 3) Crear las interacciones necesarias con 'mps' (raw shock)
df <- df %>%
  mutate(
    # 3.1) Interacciones con crecimiento de GDP (parte estática)
    lev_wins_dem_std_gdp  = lev_std * dlog_gdp,
    d2d_wins_dem_std_gdp  = dd_std  * dlog_gdp,
    
    # 3.2) Interacciones con choque raw (mps):
    #       'lev_wins_dem_std_mps' = lev_std × mps
    #       'd2d_wins_dem_std_mps' = dd_std  × mps
    lev_wins_dem_std_mps  = lev_std * mps,
    d2d_wins_dem_std_mps  = dd_std  * mps
  )

# 4) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
df_dyn <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  do({
    tmp <- .
    for (h in 0:12) {
      tmp[[paste0("cumF", h, "_dlog_capital")]] <-
        rowSums(
          map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)),
          na.rm = TRUE
        )
    }
    tmp
  }) %>%
  ungroup()

# 5) Definir controles firm‐level y de país (solo efectos fijos)
controls_firm   <- c("rsales_g_std", "size_std", "sh_current_a_std")
# No hay controles agregados en esta versión: solo FE a nivel firma y país

# 6) Verificar que las columnas existan
stopifnot(
  all(paste0("cumF", 0:12, "_dlog_capital") %in% colnames(df_dyn)),
  all(c(
    "lev_wins_dem_std_gdp", "d2d_wins_dem_std_gdp",
    "mps", "lev_wins_dem_std_mps", "d2d_wins_dem_std_mps"
  ) %in% colnames(df_dyn))
)

# 7) Estimar los modelos para h = 0…12
modelo_list <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  
  # Fórmula: 
  # cumFh_dlog_capital ~ lev_wins_dem_std_gdp + d2d_wins_dem_std_gdp + 
  #                     mps + lev_wins_dem_std_mps + d2d_wins_dem_std_mps +
  #                     controles firm + FE(name + country)
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "lev_wins_dem_std_gdp + d2d_wins_dem_std_gdp + ",
    "mps + lev_wins_dem_std_mps + d2d_wins_dem_std_mps + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  
  feols(
    fml,
    data    = df_dyn,
    cluster = ~dateq + name  # cluster ambas dimensiones
  )
})

# Asignar nombres m0…m12
names(modelo_list) <- paste0("m", 0:12)

# 8) Extraer coeficientes y errores estándar de 'mps' y 'd2d_wins_dem_std_mps'
coefs_mps <- tibble(
  horizon   = 0:12,
  beta_mps  = map_dbl(modelo_list, ~ coef(.x)["mps"]),
  se_mps    = map_dbl(modelo_list, ~ sqrt(vcov(.x)["mps", "mps"]))
)

coefs_dd_mps <- tibble(
  horizon       = 0:12,
  beta_dd_mps   = map_dbl(modelo_list, ~ coef(.x)["d2d_wins_dem_std_mps"]),
  se_dd_mps     = map_dbl(modelo_list, ~ sqrt(vcov(.x)["d2d_wins_dem_std_mps", "d2d_wins_dem_std_mps"]))
)

# Unir en un solo data.frame final
tabla_res <- coefs_mps %>%
  left_join(coefs_dd_mps, by = "horizon")

# 9) Mostrar la tabla en consola
print(tabla_res)

# 10) Graficar las respuestas dinámicas en R

# Gráfico para coef_mps
p1 <- ggplot(tabla_res, aes(x = horizon, y = beta_mps)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta_mps - 1.96 * se_mps, 
                  ymax = beta_mps + 1.96 * se_mps),
              alpha = 0.2, fill = "steelblue") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Dinámica de Respuesta: mps (choque raw)",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

# Gráfico para coef_dd_mps
p2 <- ggplot(tabla_res, aes(x = horizon, y = beta_dd_mps)) +
  geom_line(size = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_ribbon(aes(ymin = beta_dd_mps - 1.96 * se_dd_mps,
                  ymax = beta_dd_mps + 1.96 * se_dd_mps),
              alpha = 0.2, fill = "firebrick") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Dinámica de Respuesta: DD × mps",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

# Mostrar los dos gráficos uno debajo del otro
print(p1)
print(p2)



# ------------------------------------------------------------
# Figura 11: Dynamics Controlando por Inversión Rezagada (horizonte 0–12)
# ------------------------------------------------------------
# ------------------------------------------------------------
# Figura 11: Dynamics Controlando por Inversión Rezagada (horizonte 0–12)
# ------------------------------------------------------------

# 0) Cargar librerías necesarias
# Si falta alguna, instala con:
# install.packages(c("haven","zoo","dplyr","purrr","fixest","ggplot2","tidyr"))
library(haven)     # Para read_dta()
library(zoo)       # Para as.yearqtr()
library(dplyr)     # Para manipulación de datos
library(tidyr)     # Para replace_na()
library(purrr)     # Para map() y map_dfc()
library(fixest)    # Para feols()
library(ggplot2)   # Para graficar

# ------------------------------------------------------------
# 1) Leer el panel y crear variable trimestral
# ------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),    # convierte a yearqtr
    quarter = as.factor(dateq)      # factor para trimestre (opcional)
  ) %>%
  arrange(name, dateq)              # orden por firma ("name") y fecha

# ------------------------------------------------------------
# 2) Preprocesar variables para evitar log(0)/log(NA) y obtener dlog_capital, etc.
# ------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) Crear versiones “limpias” para logaritmos (si ≤0 → NA)
    capital2 = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    saleq2   = if_else(!is.na(saleq)   & saleq   > 0, saleq,   NA_real_),
    atq2     = if_else(!is.na(atq)     & atq     > 0, atq,     NA_real_),
    
    # 2.2) Calcular logaritmos de las versiones limpias (log(NA)→NA sin warning)
    log_capital    = log(capital2),
    lag_log_capital = lag(log_capital),
    log_saleq      = log(saleq2),
    lag_log_saleq   = lag(log_saleq),
    log_atq        = log(atq2),
    
    # 2.3) Diferencias de logaritmos: dlog_capital y rsales_g
    dlog_capital = log_capital - lag_log_capital,
    rsales_g     = log_saleq   - lag_log_saleq,
    
    # 2.4) Lag de dlog_capital para Ldl_capital
    Ldl_capital = lag(dlog_capital, 1),
    
    # 2.5) Estandarizar rsales_g
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.6) Estandarizar size (log_atq)
    size_std = if_else(
      !is.na(log_atq),
      (log_atq - mean(log_atq, na.rm = TRUE)) / sd(log_atq, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.7) Estandarizar liquidez (current_ratio)
    sh_current_a_std = if_else(
      !is.na(current_ratio),
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.8) Estandarizar leverage (lev) y distance‐to‐default (dd)
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  # 2.9) Eliminar columnas intermedias que ya no se usan
  select(-capital2, -saleq2, -atq2,
         -log_capital, -lag_log_capital,
         -log_saleq, -lag_log_saleq,
         -log_atq)

# ------------------------------------------------------------
# 3) Crear interacciones con mps (raw shock)
# ------------------------------------------------------------
df <- df %>%
  mutate(
    # 3.1) Interacciones con crecimiento de GDP (estático)
    lev_wins_dem_std_gdp  = lev_std * dlog_gdp,
    d2d_wins_dem_std_gdp  = dd_std  * dlog_gdp,
    
    # 3.2) Interacciones con mps (raw shock)
    lev_wins_dem_std_mps  = lev_std * mps,
    d2d_wins_dem_std_mps  = dd_std  * mps
  )

# ------------------------------------------------------------
# 4) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
# ------------------------------------------------------------
#    Para evitar problemas con do() y nombres anónimos, usamos un bucle for sobre h = 0…12.
# ------------------------------------------------------------

# Partimos de df y, dentro de cada h, añadimos la columna cumFh_dlog_capital.
df_dyn <- df %>% arrange(name, dateq)

for (h in 0:12) {
  # Nombre de la nueva columna
  new_col <- paste0("cumF", h, "_dlog_capital")
  
  # Dentro de cada grupo “name”, sumar los siguientes h términos de dlog_capital
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      # rowSums(map_dfc(0:h, ~ lead(dlog_capital, .x))) acumula dlog_capital desde t hasta t+h
      !!new_col := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# ------------------------------------------------------------
# 5) Después de crear las dinámicas, reemplazar NA en regressores por 0
#    para no eliminar filas enteras en las regresiones.
# ------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        rsales_g_std, size_std, sh_current_a_std,
        Ldl_capital,
        lev_wins_dem_std_gdp, lev_wins_dem_std_mps,
        d2d_wins_dem_std_gdp, d2d_wins_dem_std_mps
      ),
      ~ replace_na(.x, 0)
    )
  )

# ------------------------------------------------------------
# 6) Definir controles firm‐level (solo se usan en el RHS)
#    y efectos fijos a nivel firma + país (Country_x)
# ------------------------------------------------------------
controls_firm <- c("rsales_g_std", "size_std", "sh_current_a_std")
# Efectos fijos: name (firma) + Country_x (país)

# ------------------------------------------------------------
# 7) Verificar que existan las columnas que necesitamos
# ------------------------------------------------------------
stopifnot(
  all(paste0("cumF", 0:12, "_dlog_capital") %in% colnames(df_dyn)),
  all(c(
    "lev_wins_dem_std_gdp", "lev_wins_dem_std_mps",
    "d2d_wins_dem_std_gdp", "d2d_wins_dem_std_mps",
    "Ldl_capital"
  ) %in% colnames(df_dyn))
)

# ------------------------------------------------------------
# 8) (a) Dinámica para Leverage: bucle h = 0…12
# ------------------------------------------------------------
models_lev <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  
  # Fórmula para Leverage (incluye Ldl_capital):
  # cumFh_dlog_capital ~ lev_wins_dem_std_gdp + lev_wins_dem_std_mps + Ldl_capital + controles firm
  # FE a nivel: name + Country_x
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "lev_wins_dem_std_gdp + lev_wins_dem_std_mps + Ldl_capital + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  
  feols(
    fml,
    data    = df_dyn,
    cluster = ~dateq + name
  )
})

names(models_lev) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 9) Extraer coeficientes y errores estándar de ‘lev_wins_dem_std_mps’
# ------------------------------------------------------------
coefs_lev_mps <- tibble(
  horizon   = 0:12,
  beta_lev  = map_dbl(models_lev, ~ coef(.x)["lev_wins_dem_std_mps"]),
  se_lev    = map_dbl(models_lev, ~ sqrt(vcov(.x)["lev_wins_dem_std_mps", "lev_wins_dem_std_mps"]))
)

# ------------------------------------------------------------
# 10) Mostrar tabla de dinámica Leverage
# ------------------------------------------------------------
cat("\nTabla 11(a): Dinámica con Leverage (coeficiente de lev_wins_dem_std_mps)\n")
print(coefs_lev_mps)

# ------------------------------------------------------------
# 11) Graficar respuesta dinámica de Leverage
# ------------------------------------------------------------
p_lev_dyn <- ggplot(coefs_lev_mps, aes(x = horizon, y = beta_lev)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta_lev - 1.96 * se_lev,
                  ymax = beta_lev + 1.96 * se_lev),
              alpha = 0.2, fill = "steelblue") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 11(a): Dinámica IR – Leverage × mps (con Ldl_capital)",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

print(p_lev_dyn)

# ------------------------------------------------------------
# 12) (b) Dinámica para Distance‐to‐Default: bucle h = 0…12
# ------------------------------------------------------------
models_dd <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  
  # Fórmula para DD (incluye Ldl_capital):
  # cumFh_dlog_capital ~ d2d_wins_dem_std_gdp + d2d_wins_dem_std_mps + Ldl_capital + controles firm
  # FE a nivel: name + Country_x
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "d2d_wins_dem_std_gdp + d2d_wins_dem_std_mps + Ldl_capital + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  
  feols(
    fml,
    data    = df_dyn,
    cluster = ~dateq + name
  )
})

names(models_dd) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 13) Extraer coeficientes y errores estándar de ‘d2d_wins_dem_std_mps’
# ------------------------------------------------------------
coefs_dd_mps <- tibble(
  horizon      = 0:12,
  beta_dd      = map_dbl(models_dd, ~ coef(.x)["d2d_wins_dem_std_mps"]),
  se_dd        = map_dbl(models_dd, ~ sqrt(vcov(.x)["d2d_wins_dem_std_mps", "d2d_wins_dem_std_mps"]))
)

# ------------------------------------------------------------
# 14) Mostrar tabla de dinámica DD
# ------------------------------------------------------------
cat("\nTabla 11(b): Dinámica con Distance-to-Default (coeficiente de d2d_wins_dem_std_mps)\n")
print(coefs_dd_mps)

# ------------------------------------------------------------
# 15) Graficar respuesta dinámica de Distance-to-Default
# ------------------------------------------------------------
p_dd_dyn <- ggplot(coefs_dd_mps, aes(x = horizon, y = beta_dd)) +
  geom_line(size = 1, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_ribbon(aes(ymin = beta_dd - 1.96 * se_dd,
                  ymax = beta_dd + 1.96 * se_dd),
              alpha = 0.2, fill = "firebrick") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 11(b): Dinámica IR – DD × mps (con Ldl_capital)",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  theme_minimal()

print(p_dd_dyn)



# ------------------------------------------------------------
# Fin del script para Figura 11
# ------------------------------------------------------------

# *******************************************************************************
# Figure 12 (usando toda la muestra, variable “cash_ratio” en lugar de Size)
# *******************************************************************************

# 0) Cargar librerías necesarias
# Si falta alguna, instálala con install.packages(...)
library(haven)     # read_dta()
library(zoo)       # as.yearqtr()
library(dplyr)     # manipulación de datos
library(tidyr)     # replace_na()
library(purrr)     # map() y map_dfc()
library(fixest)    # feols()
library(ggplot2)   # graficar

# ------------------------------------------------------------
# 1) Leer panel completo (sin recortes de periodo) y crear variable trimestral
# ------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df_full <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),    # convierte a yearqtr
    quarter = as.factor(dateq)      # factor trimestral (opcional)
  ) %>%
  arrange(name, dateq)

# ------------------------------------------------------------
# 2) Preprocesar variables básicas: dlog_capital, Ldl_capital, controles firma‐nivel
# ------------------------------------------------------------
df_full <- df_full %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) “Limpiar” capital para logaritmos (si ≤0 → NA)
    capital2     = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    log_capital  = log(capital2),
    lag_log_cap  = lag(log_capital),
    dlog_capital = log_capital - lag_log_cap,
    Ldl_capital  = lag(dlog_capital, 1),
    
    # 2.2) Venta trimestral para controlar crecimiento de ventas
    saleq2       = if_else(!is.na(saleq) & saleq > 0, saleq, NA_real_),
    log_sales    = log(saleq2),
    rsales_g     = log_sales - lag(log_sales),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.3) Estandarizar apalancamiento y distance‐to‐default
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  # 2.4) Eliminar variables intermedias ya no necesarias
  select(-capital2, -lag_log_cap, -saleq2, -rsales_g)

# ------------------------------------------------------------
# 3) Generar proxies de liquidez (“cash_ratio”) y sus interacciones
# ------------------------------------------------------------
df_full <- df_full %>%
  group_by(dateq) %>%
  mutate(
    # 3.1) Estandarizar cash_ratio para cada trimestre
    cash_ratio_std = if_else(
      !is.na(cash_ratio),
      (cash_ratio - mean(cash_ratio, na.rm = TRUE)) / sd(cash_ratio, na.rm = TRUE),
      NA_real_
    ),
    
    # 3.2) Dummy HighCash = 1 si cash_ratio ≥ mediana(trimestre), 0 si < mediana
    median_cr   = median(cash_ratio, na.rm = TRUE),
    HighCash    = as.numeric(cash_ratio >= median_cr)
  ) %>%
  ungroup() %>%
  mutate(
    # 3.3) Interacciones con choque raw (mps) y con crecimiento de GDP (dlog_gdp)
    lev_gdp    = lev_std * dlog_gdp,
    dd_gdp     = dd_std  * dlog_gdp,
    lev_mps    = lev_std * mps,
    dd_mps     = dd_std  * mps,
    
    # 3.4) Interacciones liquidez × mps:
    cash_mps      = cash_ratio_std * mps,   # liquidez continua
    highcash_mps  = HighCash * mps          # liquidez dummy
  ) %>%
  # 3.5) Eliminar variable intermedia
  select(-median_cr)

# ------------------------------------------------------------
# 4) Crear variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
# ------------------------------------------------------------
df_dyn <- df_full %>% arrange(name, dateq)

for (h in 0:12) {
  new_col <- paste0("cumF", h, "_dlog_capital")
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      !!new_col := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# ------------------------------------------------------------
# 5) Reemplazar NA en regressors (para no perder filas salvo si falta cumFh)
# ------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        rsales_g_std, Ldl_capital,
        lev_gdp, lev_mps,
        dd_gdp, dd_mps,
        cash_ratio_std, cash_mps, highcash_mps
      ),
      ~ replace_na(.x, 0)
    )
  )

# ------------------------------------------------------------
# 6) Definir controles firma‐nivel y efectos fijos
#    - Controles: rsales_g_std  
#    - FEs: name (firma) + Country_x (país)
# ------------------------------------------------------------
controls_firm <- c("rsales_g_std")
# Efectos fijos se colocarán en la fórmula con “| name + Country_x”

# ------------------------------------------------------------
# 7) Verificar que existan las columnas que necesitamos
# ------------------------------------------------------------
stopifnot(
  all(paste0("cumF", 0:12, "_dlog_capital") %in% colnames(df_dyn)),
  all(c(
    "lev_gdp", "lev_mps",
    "dd_gdp", "dd_mps",
    "Ldl_capital",
    "cash_mps", "highcash_mps"
  ) %in% colnames(df_dyn))
)

# ------------------------------------------------------------
# 8) Estimar dinámicas para Toda la muestra (h = 0…12)
#    (a) Liquidez continua: “cash_mps”
# ------------------------------------------------------------
models_cashc <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "cash_mps + lev_gdp + lev_mps + dd_gdp + dd_mps + Ldl_capital + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_cashc) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 9) Extraer coeficientes y errores de “cash_mps” (liquidez continua)
# ------------------------------------------------------------
coefs_cashc <- tibble(
  horizon    = 0:12,
  beta_cashc = map_dbl(models_cashc, ~ coef(.x)["cash_mps"]),
  se_cashc   = map_dbl(models_cashc, ~ sqrt(vcov(.x)["cash_mps", "cash_mps"]))
)

cat("\nFigura 12(a): Toda la muestra – coeficiente cash_mps (liquidez continua × mps)\n")
print(coefs_cashc)

# ------------------------------------------------------------
# 10) Graficar respuesta dinámica de “cash_mps”
# ------------------------------------------------------------
p_cashc <- ggplot(coefs_cashc, aes(x = horizon, y = beta_cashc)) +
  geom_line(size = 1, color = "forestgreen") +
  geom_point(size = 2, color = "forestgreen") +
  geom_ribbon(aes(ymin = beta_cashc - 1.96 * se_cashc,
                  ymax = beta_cashc + 1.96 * se_cashc),
              alpha = 0.2, fill = "forestgreen") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Fig 12(a): Dinámica – cash_mps (liquidez continua × mps)",
    x     = "Horizonte (trimestres)",
    y     = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_cashc)

# ------------------------------------------------------------
# 11) Estimar dinámicas (b) Liquidez dummy: “highcash_mps”
# ------------------------------------------------------------
models_highc <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "highcash_mps + lev_gdp + lev_mps + dd_gdp + dd_mps + Ldl_capital + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_highc) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 12) Extraer coeficientes y errores de “highcash_mps” (liquidez dummy)
# ------------------------------------------------------------
coefs_highc <- tibble(
  horizon     = 0:12,
  beta_highc  = map_dbl(models_highc, ~ coef(.x)["highcash_mps"]),
  se_highc    = map_dbl(models_highc, ~ sqrt(vcov(.x)["highcash_mps", "highcash_mps"]))
)

cat("\nFigura 12(b): Toda la muestra – coeficiente highcash_mps (liquidez dummy × mps)\n")
print(coefs_highc)

# ------------------------------------------------------------
# 13) Graficar respuesta dinámica de “highcash_mps”
# ------------------------------------------------------------
p_highc <- ggplot(coefs_highc, aes(x = horizon, y = beta_highc)) +
  geom_line(size = 1, color = "darkorange") +
  geom_point(size = 2, color = "darkorange") +
  geom_ribbon(aes(ymin = beta_highc - 1.96 * se_highc,
                  ymax = beta_highc + 1.96 * se_highc),
              alpha = 0.2, fill = "darkorange") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Fig 12(b): Dinámica – highcash_mps (liquidez dummy × mps)",
    x     = "Horizonte (trimestres)",
    y     = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_highc)

# ------------------------------------------------------------
# Fin del script para Figure 12 con toda la muestra (sin recortes de periodo)
# ------------------------------------------------------------





# *******************************************************************************
# Figura 13: Joint Dynamics en Posición Financiera (Leverage o DD) y Liquidez
# (Se usa TODA la muestra, sin recortes de periodos antiguos)
# *******************************************************************************

# 0) Cargar librerías necesarias
# Si falta alguna, instálala con: install.packages(c("haven","zoo","dplyr","tidyr","purrr","fixest","ggplot2"))
library(haven)     # read_dta()
library(zoo)       # as.yearqtr()
library(dplyr)     # Manipulación de datos
library(tidyr)     # replace_na()
library(purrr)     # map(), map_dfc()
library(fixest)    # feols()
library(ggplot2)   # Graficar

# ------------------------------------------------------------
# 1) Leer el panel “construct_panel_data_final.dta” (incluye cash_ratio)
# ------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),    # convierte a trimestral
    quarter = as.factor(dateq)      # factor trimestral (opcional)
  ) %>%
  arrange(name, dateq)

# ------------------------------------------------------------
# 2) Preprocesar variables básicas:
#    - dlog_capital, Ldl_capital
#    - Controles firm‐level (rsales_g_std)
#    - Apalancamiento (lev_std), DD (dd_std)
#    - Liquidez (cash_ratio_std, HighCash)
# ------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) Evitar log(0) en “capital”
    capital2     = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    log_capital  = log(capital2),
    lag_log_cap  = lag(log_capital),
    dlog_capital = log_capital - lag_log_cap,
    Ldl_capital  = lag(dlog_capital, 1),
    
    # 2.2) Crecimiento de ventas para control rsales_g_std
    saleq2       = if_else(!is.na(saleq) & saleq > 0, saleq, NA_real_),
    log_sales    = log(saleq2),
    rsales_g     = log_sales - lag(log_sales),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.3) Estandarizar apalancamiento y DD
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.4) Estandarizar liquidez (cash_ratio) por trimestre
    cash_ratio2    = if_else(!is.na(cash_ratio), cash_ratio, NA_real_),
    cash_ratio_std = (cash_ratio2 - mean(cash_ratio2, na.rm = TRUE)) / sd(cash_ratio2, na.rm = TRUE),
    
    # 2.5) Dummy HighCash = 1 si cash_ratio ≥ mediana(cash_ratio) de ese trimestre
    median_cr = median(cash_ratio2, na.rm = TRUE),
    HighCash  = as.numeric(cash_ratio2 >= median_cr)
  ) %>%
  ungroup() %>%
  # 2.6) Eliminar variables intermedias no necesarias
  select(-capital2, -lag_log_cap, -saleq2, -rsales_g, -log_sales, -cash_ratio2, -median_cr)

# ------------------------------------------------------------
# 3) Interacciones con shocks:
#    - Financieras: lev_gdp, lev_mps, dd_gdp, dd_mps
#    - Liquidez: cash_mps (continua), highcash_mps (dummy)
# ------------------------------------------------------------
df <- df %>%
  mutate(
    # Choques macro y financieros
    lev_gdp = lev_std * dlog_gdp,
    dd_gdp  = dd_std  * dlog_gdp,
    lev_mps = lev_std * mps,
    dd_mps  = dd_std  * mps,
    
    # Interacciones liquidez × raw shock (mps)
    cash_mps     = cash_ratio_std * mps,    # liquidez continua
    highcash_mps = HighCash * mps           # liquidez dummy
  )

# ------------------------------------------------------------
# 4) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
# ------------------------------------------------------------
df_dyn <- df %>% arrange(name, dateq)
for (h in 0:12) {
  new_col <- paste0("cumF", h, "_dlog_capital")
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      !!new_col := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# ------------------------------------------------------------
# 5) Reemplazar NA en todos los regresores de RHS (para no eliminar filas salvo faltante de cumFh)
# ------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        rsales_g_std, Ldl_capital,
        lev_gdp, lev_mps,
        dd_gdp, dd_mps,
        cash_ratio_std, cash_mps, highcash_mps
      ),
      ~ replace_na(.x, 0)
    )
  )

# ------------------------------------------------------------
# 6) Definir controles y FEs:
#    - Controles de firma: rsales_g_std
#    - FEs: name (firma) + Country_x (país)
# ------------------------------------------------------------
controls_firm <- c("rsales_g_std")
# En la fórmula R: se escribirá “| name + Country_x”

# ------------------------------------------------------------
# 7) (a) Joint Dynamics: Leverage & Liquidez (h = 0…12)
#    Especificación:
#      cumFh_dlog_capital ~ 
#         cash_ratio_std + cash_mps + highcash_mps +
#         lev_gdp + lev_mps + Ldl_capital + rsales_g_std
#      | name + Country_x, cluster(~dateq + name)
# 
#    NOTA: incluimos tanto cash_ratio_std (nivel) como sus interacciones cash_mps / highcash_mps
# ------------------------------------------------------------
models_lev_liq <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "cash_ratio_std + cash_mps + highcash_mps + ",
    "lev_gdp + lev_mps + Ldl_capital + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_lev_liq) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 8) Extraer coeficientes y errores para “cash_mps” y “lev_mps”
# ------------------------------------------------------------
coefs_lev_liq <- tibble(
  horizon      = 0:12,
  beta_cash    = map_dbl(models_lev_liq, ~ coef(.x)["cash_mps"]),
  se_cash      = map_dbl(models_lev_liq, ~ sqrt(vcov(.x)["cash_mps", "cash_mps"])),
  beta_highc   = map_dbl(models_lev_liq, ~ coef(.x)["highcash_mps"]),
  se_highc     = map_dbl(models_lev_liq, ~ sqrt(vcov(.x)["highcash_mps", "highcash_mps"])),
  beta_lev     = map_dbl(models_lev_liq, ~ coef(.x)["lev_mps"]),
  se_lev       = map_dbl(models_lev_liq, ~ sqrt(vcov(.x)["lev_mps", "lev_mps"]))
)

cat("\nFigura 13(a): Dynamics – Leverage & Liquidez\n")
print(coefs_lev_liq)

# ------------------------------------------------------------
# 9) Graficar respuestas dinámicas combinadas:
#    - cash_mps (línea azul), highcash_mps (línea naranja punteada)
#    - lev_mps (línea roja)
# ------------------------------------------------------------
p_lev_liq <- ggplot(coefs_lev_liq, aes(x = horizon)) +
  # Interacción cash_mps (liquidez continua x mps)
  geom_line(aes(y = beta_cash), size = 1, color = "blue") +
  geom_ribbon(aes(ymin = beta_cash - 1.96 * se_cash,
                  ymax = beta_cash + 1.96 * se_cash),
              fill = "blue", alpha = 0.2) +
  # Interacción highcash_mps (liquidez dummy x mps)
  geom_line(aes(y = beta_highc), size = 1, linetype = "dashed", color = "orange") +
  geom_ribbon(aes(ymin = beta_highc - 1.96 * se_highc,
                  ymax = beta_highc + 1.96 * se_highc),
              fill = "orange", alpha = 0.2) +
  # Interacción lev_mps (apalancamiento x mps)
  geom_line(aes(y = beta_lev), size = 1, color = "red") +
  geom_ribbon(aes(ymin = beta_lev - 1.96 * se_lev,
                  ymax = beta_lev + 1.96 * se_lev),
              fill = "red", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 13(a): Joint Dynamics – Liquidez (blue/orange) y Leverage (red)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_lev_liq)

# ------------------------------------------------------------
# 10) (b) Joint Dynamics: DD & Liquidez (h = 0…12)
#     Especificación:
#       cumFh_dlog_capital ~ 
#          cash_ratio_std + cash_mps + highcash_mps +
#          dd_gdp + dd_mps + Ldl_capital + rsales_g_std
#       | name + Country_x, cluster(~dateq + name)
# ------------------------------------------------------------
models_dd_liq <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "cash_ratio_std + cash_mps + highcash_mps + ",
    "dd_gdp + dd_mps + Ldl_capital + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_dd_liq) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 11) Extraer coeficientes y errores para “cash_mps” y “dd_mps”
# ------------------------------------------------------------
coefs_dd_liq <- tibble(
  horizon     = 0:12,
  beta_cash   = map_dbl(models_dd_liq, ~ coef(.x)["cash_mps"]),
  se_cash     = map_dbl(models_dd_liq, ~ sqrt(vcov(.x)["cash_mps", "cash_mps"])),
  beta_highc  = map_dbl(models_dd_liq, ~ coef(.x)["highcash_mps"]),
  se_highc    = map_dbl(models_dd_liq, ~ sqrt(vcov(.x)["highcash_mps", "highcash_mps"])),
  beta_dd     = map_dbl(models_dd_liq, ~ coef(.x)["dd_mps"]),
  se_dd       = map_dbl(models_dd_liq, ~ sqrt(vcov(.x)["dd_mps", "dd_mps"]))
)

cat("\nFigura 13(b): Dynamics – DD & Liquidez\n")
print(coefs_dd_liq)

# ------------------------------------------------------------
# 12) Graficar respuestas dinámicas combinadas:
#    - cash_mps (línea azul), highcash_mps (línea naranja punteada)
#    - dd_mps (línea verde)
# ------------------------------------------------------------
p_dd_liq <- ggplot(coefs_dd_liq, aes(x = horizon)) +
  # Interacción cash_mps (liquidez continua x mps)
  geom_line(aes(y = beta_cash), size = 1, color = "blue") +
  geom_ribbon(aes(ymin = beta_cash - 1.96 * se_cash,
                  ymax = beta_cash + 1.96 * se_cash),
              fill = "blue", alpha = 0.2) +
  # Interacción highcash_mps (liquidez dummy x mps)
  geom_line(aes(y = beta_highc), size = 1, linetype = "dashed", color = "orange") +
  geom_ribbon(aes(ymin = beta_highc - 1.96 * se_highc,
                  ymax = beta_highc + 1.96 * se_highc),
              fill = "orange", alpha = 0.2) +
  # Interacción dd_mps (distance-to-default x mps)
  geom_line(aes(y = beta_dd), size = 1, color = "darkgreen") +
  geom_ribbon(aes(ymin = beta_dd - 1.96 * se_dd,
                  ymax = beta_dd + 1.96 * se_dd),
              fill = "darkgreen", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 13(b): Joint Dynamics – Liquidez (blue/orange) y DD (verde)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_dd_liq)

# ------------------------------------------------------------
# Fin del script para Figura 13 adaptada (Leverage/DD y Liquidez)
# ------------------------------------------------------------



# *******************************************************************************
# Figura 14: Joint Dynamics en Posición Financiera (Leverage o DD) y Edad
# (Se usa TODA la muestra, sin recortes de periodos antiguos)
# *******************************************************************************

# 0) Cargar librerías necesarias
# Si falta alguna, instálala con:
# install.packages(c("haven","zoo","dplyr","tidyr","purrr","fixest","ggplot2"))
library(haven)     # read_dta()
library(zoo)       # as.yearqtr()
library(dplyr)     # Manipulación de datos
library(tidyr)     # replace_na()
library(purrr)     # map(), map_dfc()
library(fixest)    # feols()
library(ggplot2)   # Graficar

# ------------------------------------------------------------
# 1) Leer el panel “construct_panel_data_final.dta” (incluye variables de edad)
# ------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),   # convierte a trimestral
    quarter = as.factor(dateq)     # (opcional) factor para trimestre
  ) %>%
  arrange(name, dateq)






# *******************************************************************************
# Figura 14: Joint Dynamics en Posición Financiera (Leverage o DD) y Edad
# (Usando TODA la muestra, sin recortes de periodos antiguos ni Gran Depresión)
# *******************************************************************************

# 0) Cargar librerías necesarias
library(haven)     # read_dta()
library(zoo)       # as.yearqtr()
library(dplyr)     # Manipulaciones de datos
library(tidyr)     # replace_na()
library(purrr)     # map(), map_dfc()
library(fixest)    # feols()
library(ggplot2)   # Graficar

# ------------------------------------------------------------
# 1) Leer el panel “construct_panel_data_final.dta”
# ------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),   # convierte la variable `dateq` a trimestral
    quarter = as.factor(dateq)     # (opcional) factor para trimestre
  ) %>%
  arrange(name, dateq)

# ------------------------------------------------------------
# 2) Preprocesar variables firm‐level:
#    - dlog_capital, Ldl_capital
#    - rsales_g_std (crecimiento ventas estandarizado)
#    - lev_std, dd_std (apalancamiento y DD estandarizados)
#    - Construir dummies de edad: AgeMiddle / AgeOld
# ------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) Evitar log(0) en “capital”
    capital2     = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    log_capital  = log(capital2),
    lag_log_cap  = lag(log_capital),
    dlog_capital = log_capital - lag_log_cap,
    Ldl_capital  = lag(dlog_capital, 1),
    
    # 2.2) Crecimiento de ventas para control
    saleq2       = if_else(!is.na(saleq) & saleq > 0, saleq, NA_real_),
    log_sales    = log(saleq2),
    rsales_g     = log_sales - lag(log_sales),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.3) Estandarizar apalancamiento y DD
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  # 2.4) Eliminar variables temporales intermedias
  select(-capital2, -lag_log_cap, -saleq2, -log_sales, -rsales_g)

# ------------------------------------------------------------
# 3) Construir percentiles de “age_sample” por trimestre y dummies de edad
#    (AgeMiddle / AgeOld), usando siempre todo el panel
# ------------------------------------------------------------
df <- df %>%
  group_by(dateq) %>%
  mutate(
    age_p25   = quantile(age_sample, 0.25, na.rm = TRUE),
    age_p75   = quantile(age_sample, 0.75, na.rm = TRUE),
    
    # 3.1) Dummy AgeMiddle: entre percentil25 y percentil75
    AgeMiddle = as.numeric(age_sample >= age_p25 & age_sample < age_p75),
    # 3.2) Dummy AgeOld: en percentil75 o más
    AgeOld    = as.numeric(age_sample >= age_p75)
  ) %>%
  ungroup() %>%
  # 3.3) No necesitamos conservar age_p25/age_p75 tras construir dummies
  select(-age_p25, -age_p75)

# ------------------------------------------------------------
# 4) Crear interacciones financieras y de edad:
#    - lev_gdp, lev_wide, lev_mps
#    - dd_gdp,  dd_wide,  dd_mps
#    - AgeMiddle_gdp, AgeOld_gdp
#    - AgeMiddle_wide, AgeOld_wide
# ------------------------------------------------------------
df <- df %>%
  mutate(
    # 4.1) Interacciones de apalancamiento y DD
    lev_gdp  = lev_std * dlog_gdp,
    lev_wide = lev_std * mpso,   # shock ortogonalizado
    lev_mps  = lev_std * mps,    # shock raw
    dd_gdp   = dd_std  * dlog_gdp,
    dd_wide  = dd_std  * mpso,
    dd_mps   = dd_std  * mps,
    
    # 4.2) Interacciones de edad
    AgeMiddle_gdp  = AgeMiddle * dlog_gdp,
    AgeOld_gdp     = AgeOld    * dlog_gdp,
    AgeMiddle_wide = AgeMiddle * mpso,
    AgeOld_wide    = AgeOld    * mpso
  )

# ------------------------------------------------------------
# 5) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
# ------------------------------------------------------------
df_dyn <- df %>% arrange(name, dateq)
for (h in 0:12) {
  colname <- paste0("cumF", h, "_dlog_capital")
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      # Suma acumulada de dlog_capital desde t hasta t+h dentro de cada firma
      !!colname := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# ------------------------------------------------------------
# 6) Reemplazar NA en regresores de RHS por 0
#    (para no eliminar filas salvo que falte cumFh_dlog_capital)
# ------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        # Controles y rezagos
        "rsales_g_std", "Ldl_capital",
        # Interacciones financieras
        "lev_gdp", "lev_wide", "lev_mps",
        "dd_gdp",  "dd_wide",  "dd_mps",
        # Dummies y sus interacciones de edad
        "AgeMiddle", "AgeOld",
        "AgeMiddle_gdp", "AgeOld_gdp",
        "AgeMiddle_wide", "AgeOld_wide"
      ),
      ~ replace_na(.x, 0)
    )
  )

# ------------------------------------------------------------
# 7) Definir controles firma‐level y efectos fijos
#    - Controles de firma: rsales_g_std
#    - FEs fijos: name (firma) + Country_x (país)
# ------------------------------------------------------------
controls_firm <- c("rsales_g_std")
# En cada fórmula: “| name + Country_x”

# ------------------------------------------------------------
# 8) (a) Joint Dynamics: Leverage & Edad (h = 0…12)
#    Especificación:
#      cumFh_dlog_capital ~
#        AgeMiddle_gdp + AgeOld_gdp +
#        AgeMiddle_wide + AgeOld_wide +
#        AgeMiddle + AgeOld +
#        lev_gdp + lev_wide + lev_mps +
#        Ldl_capital + rsales_g_std
#      | name + Country_x, cluster(~dateq + name)
# ------------------------------------------------------------
models_lev_age <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    # Variables de edad
    "AgeMiddle_gdp + AgeOld_gdp + ",
    "AgeMiddle_wide + AgeOld_wide + ",
    "AgeMiddle + AgeOld + ",
    # Interacciones financieras
    "lev_gdp + lev_wide + lev_mps + ",
    # Rezago de inversión y control
    "Ldl_capital + ", paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_lev_age) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 9) Extraer coeficientes y errores para:
#    - AgeMiddle_wide, AgeOld_wide
#    - lev_wide
# ------------------------------------------------------------
coefs_lev_age <- tibble(
  horizon          = 0:12,
  beta_mid_wide    = map_dbl(models_lev_age, ~ coef(.x)["AgeMiddle_wide"]),
  se_mid_wide      = map_dbl(models_lev_age, ~ sqrt(vcov(.x)["AgeMiddle_wide", "AgeMiddle_wide"])),
  beta_old_wide    = map_dbl(models_lev_age, ~ coef(.x)["AgeOld_wide"]),
  se_old_wide      = map_dbl(models_lev_age, ~ sqrt(vcov(.x)["AgeOld_wide", "AgeOld_wide"])),
  beta_lev_wide    = map_dbl(models_lev_age, ~ coef(.x)["lev_wide"]),
  se_lev_wide      = map_dbl(models_lev_age, ~ sqrt(vcov(.x)["lev_wide", "lev_wide"]))
)

cat("\nFigura 14(a): Joint Dynamics – Leverage & Edad\n")
print(coefs_lev_age)

# ------------------------------------------------------------
# 10) Graficar “AgeMiddle_wide” (azul), “AgeOld_wide” (naranja) y “Leverage × wide” (rojo)
# ------------------------------------------------------------
p_lev_age <- ggplot(coefs_lev_age, aes(x = horizon)) +
  # AgeMiddle_wide
  geom_line(aes(y = beta_mid_wide), size = 1, color = "blue") +
  geom_ribbon(aes(ymin = beta_mid_wide - 1.96 * se_mid_wide,
                  ymax = beta_mid_wide + 1.96 * se_mid_wide),
              fill = "blue", alpha = 0.2) +
  # AgeOld_wide
  geom_line(aes(y = beta_old_wide), size = 1, color = "orange") +
  geom_ribbon(aes(ymin = beta_old_wide - 1.96 * se_old_wide,
                  ymax = beta_old_wide + 1.96 * se_old_wide),
              fill = "orange", alpha = 0.2) +
  # Leverage × wide
  geom_line(aes(y = beta_lev_wide), size = 1, color = "red") +
  geom_ribbon(aes(ymin = beta_lev_wide - 1.96 * se_lev_wide,
                  ymax = beta_lev_wide + 1.96 * se_lev_wide),
              fill = "red", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 14(a): Joint Dynamics – Edad (AgeMiddle_wide azul, AgeOld_wide naranja) y Leverage_wide (rojo)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_lev_age)

# ------------------------------------------------------------
# 11) (b) Joint Dynamics: DD & Edad (h = 0…12)
#     Especificación:
#       cumFh_dlog_capital ~
#         AgeMiddle_gdp + AgeOld_gdp +
#         AgeMiddle_wide + AgeOld_wide +
#         AgeMiddle + AgeOld +
#         dd_gdp + dd_wide + dd_mps +
#         Ldl_capital + rsales_g_std
#       | name + Country_x, cluster(~dateq + name)
# ------------------------------------------------------------
models_dd_age <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    # Edad
    "AgeMiddle_gdp + AgeOld_gdp + ",
    "AgeMiddle_wide + AgeOld_wide + ",
    "AgeMiddle + AgeOld + ",
    # DD
    "dd_gdp + dd_wide + dd_mps + ",
    # Rezago y control
    "Ldl_capital + ", paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_dd_age) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 12) Extraer coeficientes y errores para:
#    - AgeMiddle_wide, AgeOld_wide
#    - dd_wide
# ------------------------------------------------------------
coefs_dd_age <- tibble(
  horizon          = 0:12,
  beta_mid_wide    = map_dbl(models_dd_age, ~ coef(.x)["AgeMiddle_wide"]),
  se_mid_wide      = map_dbl(models_dd_age, ~ sqrt(vcov(.x)["AgeMiddle_wide", "AgeMiddle_wide"])),
  beta_old_wide    = map_dbl(models_dd_age, ~ coef(.x)["AgeOld_wide"]),
  se_old_wide      = map_dbl(models_dd_age, ~ sqrt(vcov(.x)["AgeOld_wide", "AgeOld_wide"])),
  beta_dd_wide     = map_dbl(models_dd_age, ~ coef(.x)["dd_wide"]),
  se_dd_wide       = map_dbl(models_dd_age, ~ sqrt(vcov(.x)["dd_wide", "dd_wide"]))
)

cat("\nFigura 14(b): Joint Dynamics – DD & Edad\n")
print(coefs_dd_age)

# ------------------------------------------------------------
# 13) Graficar “AgeMiddle_wide” (azul), “AgeOld_wide” (naranja) y “DD × wide” (verde)
# ------------------------------------------------------------
p_dd_age <- ggplot(coefs_dd_age, aes(x = horizon)) +
  # AgeMiddle_wide
  geom_line(aes(y = beta_mid_wide), size = 1, color = "blue") +
  geom_ribbon(aes(ymin = beta_mid_wide - 1.96 * se_mid_wide,
                  ymax = beta_mid_wide + 1.96 * se_mid_wide),
              fill = "blue", alpha = 0.2) +
  # AgeOld_wide
  geom_line(aes(y = beta_old_wide), size = 1, color = "orange") +
  geom_ribbon(aes(ymin = beta_old_wide - 1.96 * se_old_wide,
                  ymax = beta_old_wide + 1.96 * se_old_wide),
              fill = "orange", alpha = 0.2) +
  # DD × wide (verde)
  geom_line(aes(y = beta_dd_wide), size = 1, color = "darkgreen") +
  geom_ribbon(aes(ymin = beta_dd_wide - 1.96 * se_dd_wide,
                  ymax = beta_dd_wide + 1.96 * se_dd_wide),
              fill = "darkgreen", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 14(b): Joint Dynamics – Edad (AgeMiddle_wide azul, AgeOld_wide naranja) y DD_wide (verde)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_dd_age)

# ------------------------------------------------------------
# Fin del script para Figura 14 (Leverage/DD y Edad)
# ------------------------------------------------------------





# *******************************************************************************
# Figura 16: Joint Dynamics en Posición Financiera (Leverage o DD) y Liquidez
# (Usando TODO el panel sin recortes de periodos antiguos)
# *******************************************************************************

# 0) Cargar librerías necesarias
#    Si falta alguna, instalar con: install.packages(c("haven","zoo","dplyr","tidyr","purrr","fixest","ggplot2"))
library(haven)     # read_dta()
library(zoo)       # as.yearqtr()
library(dplyr)     # Manipulación de datos
library(tidyr)     # replace_na()
library(purrr)     # map(), map_dfc()
library(fixest)    # feols()
library(ggplot2)   # Graficar

# ------------------------------------------------------------
# 1) Leer el panel “construct_panel_data_final.dta”
#    (asegúrate de que incluye “current_ratio” para liquidez)
# ------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq   = as.yearqtr(dateq),   # convierte a trimestral
    quarter = as.factor(dateq)     # (opcional) factor para trimestre
  ) %>%
  arrange(name, dateq)

# ------------------------------------------------------------
# 2) Preprocesar variables firm‐level:
#    - dlog_capital, Ldl_capital
#    - rsales_g_std  (crecimiento de ventas estandarizado)
#    - lev_std, dd_std  (apalancamiento y DD estandarizados)
#    - liq_std = current_ratio estandarizada
# ------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) Evitar log(0) en “capital”
    capital2     = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    log_capital  = log(capital2),
    lag_log_cap  = lag(log_capital),
    dlog_capital = log_capital - lag_log_cap,
    Ldl_capital  = lag(dlog_capital, 1),
    
    # 2.2) Crecimiento de ventas para rsales_g_std
    saleq2       = if_else(!is.na(saleq) & saleq > 0, saleq, NA_real_),
    log_sales    = log(saleq2),
    rsales_g     = log_sales - lag(log_sales),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.3) Estandarizar apalancamiento y DD
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.4) Estandarizar liquidez a partir de current_ratio
    liq_std = if_else(
      !is.na(current_ratio),
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  # 2.5) Eliminar variables temporales intermedias
  select(-capital2, -lag_log_cap, -saleq2, -log_sales, -rsales_g)

# ------------------------------------------------------------
# 3) Crear interacciones financieras y de liquidez:
#    - Leverage:   lev_gdp, lev_wide, lev_mps
#    - DD:         dd_gdp, dd_wide, dd_mps
#    - Liquidez:   liq_gdp = liq_std × dlog_gdp
#                  liq_wide = liq_std × mpso
#                  liq_mps  = liq_std × mps
# ------------------------------------------------------------
df <- df %>%
  mutate(
    # 3.1) Interacciones de apalancamiento
    lev_gdp   = lev_std * dlog_gdp,
    lev_wide  = lev_std * mpso,    # choque ortogonalizado
    lev_mps   = lev_std * mps,     # choque raw
    
    # 3.2) Interacciones de DD
    dd_gdp    = dd_std  * dlog_gdp,
    dd_wide   = dd_std  * mpso,
    dd_mps    = dd_std  * mps,
    
    # 3.3) Interacciones de liquidez
    liq_gdp   = liq_std * dlog_gdp,
    liq_wide  = liq_std * mpso,
    liq_mps   = liq_std * mps
  )

# ------------------------------------------------------------
# 4) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
# ------------------------------------------------------------
df_dyn <- df %>% arrange(name, dateq)
for (h in 0:12) {
  new_col <- paste0("cumF", h, "_dlog_capital")
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      # Suma acumulada de dlog_capital desde t hasta t+h
      !!new_col := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# ------------------------------------------------------------
# 5) Reemplazar NA en regresores de RHS por 0
#    (para no eliminar filas salvo que falte cumFh_dlog_capital)
# ------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        # Controles y rezagos
        "rsales_g_std", "Ldl_capital",
        # Interacciones leverage
        "lev_gdp",   "lev_wide",   "lev_mps",
        # Interacciones DD
        "dd_gdp",    "dd_wide",    "dd_mps",
        # Interacciones liquidez
        "liq_std",   "liq_gdp",    "liq_wide",   "liq_mps"
      ),
      ~ replace_na(.x, 0)
    )
  )

# ------------------------------------------------------------
# 6) Definir controles firma‐level y efectos fijos
#    - Controles: rsales_g_std
#    - Efectos fijos: name (firma) + Country_x (país)
# ------------------------------------------------------------
controls_firm <- c("rsales_g_std")
# En cada fórmula: “| name + Country_x”

# ------------------------------------------------------------
# 7) (a) Joint Dynamics: Leverage & Liquidez (h = 0…12)
#    Especificación:
#      cumFh_dlog_capital ~
#        lev_gdp  + lev_wide  + lev_mps  +
#        liq_gdp  + liq_wide  + liq_mps  +
#        Ldl_capital + rsales_g_std
#      | name + Country_x, cluster(~dateq + name)
# ------------------------------------------------------------
models_lev_liq <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    # Interacciones leverage
    "lev_gdp + lev_wide + lev_mps + ",
    # Interacciones liquidez
    "liq_gdp + liq_wide + liq_mps + ",
    # Rezago y control
    "Ldl_capital + ", paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_lev_liq) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 8) Extraer coeficientes y errores de:
#    • lev_wide   (Leverage × mpso)
#    • liq_wide   (Liquidez × mpso)
# ------------------------------------------------------------
coefs_lev_liq <- tibble(
  horizon         = 0:12,
  beta_lev_wide   = map_dbl(models_lev_liq, ~ coef(.x)["lev_wide"]),
  se_lev_wide     = map_dbl(models_lev_liq, ~ sqrt(vcov(.x)["lev_wide", "lev_wide"])),
  beta_liq_wide   = map_dbl(models_lev_liq, ~ coef(.x)["liq_wide"]),
  se_liq_wide     = map_dbl(models_lev_liq, ~ sqrt(vcov(.x)["liq_wide", "liq_wide"]))
)

cat("\nFigura 16(a): Joint Dynamics – Leverage & Liquidez\n")
print(coefs_lev_liq)

# ------------------------------------------------------------
# 9) Graficar coeficientes sobre “lev_wide” (rojo) y “liq_wide” (azul)
# ------------------------------------------------------------
p_lev_liq <- ggplot(coefs_lev_liq, aes(x = horizon)) +
  # Leverage × wide (rojo)
  geom_line(aes(y = beta_lev_wide), size = 1, color = "red") +
  geom_ribbon(aes(ymin = beta_lev_wide - 1.96 * se_lev_wide,
                  ymax = beta_lev_wide + 1.96 * se_lev_wide),
              fill = "red", alpha = 0.2) +
  # Liquidez × wide (azul)
  geom_line(aes(y = beta_liq_wide), size = 1, color = "blue") +
  geom_ribbon(aes(ymin = beta_liq_wide - 1.96 * se_liq_wide,
                  ymax = beta_liq_wide + 1.96 * se_liq_wide),
              fill = "blue", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 16(a): Joint Dynamics – Leverage (rojo) y Liquidez (azul)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_lev_liq)

# ------------------------------------------------------------
# 10) (b) Joint Dynamics: DD & Liquidez (h = 0…12)
#     Especificación:
#       cumFh_dlog_capital ~
#         dd_gdp  + dd_wide  + dd_mps  +
#         liq_gdp + liq_wide + liq_mps +
#         Ldl_capital + rsales_g_std
#       | name + Country_x, cluster(~dateq + name)
# ------------------------------------------------------------
models_dd_liq <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    # Interacciones DD
    "dd_gdp + dd_wide + dd_mps + ",
    # Interacciones liquidez
    "liq_gdp + liq_wide + liq_mps + ",
    # Rezago y control
    "Ldl_capital + ", paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_dd_liq) <- paste0("m", 0:12)

# ------------------------------------------------------------
# 11) Extraer coeficientes y errores de:
#    • dd_wide   (DD × mpso)
#    • liq_wide  (Liquidez × mpso)
# ------------------------------------------------------------
coefs_dd_liq <- tibble(
  horizon         = 0:12,
  beta_dd_wide    = map_dbl(models_dd_liq, ~ coef(.x)["dd_wide"]),
  se_dd_wide      = map_dbl(models_dd_liq, ~ sqrt(vcov(.x)["dd_wide", "dd_wide"])),
  beta_liq_wide   = map_dbl(models_dd_liq, ~ coef(.x)["liq_wide"]),
  se_liq_wide     = map_dbl(models_dd_liq, ~ sqrt(vcov(.x)["liq_wide", "liq_wide"]))
)

cat("\nFigura 16(b): Joint Dynamics – DD & Liquidez\n")
print(coefs_dd_liq)

# ------------------------------------------------------------
# 12) Graficar coeficientes sobre “dd_wide” (verde) y “liq_wide” (azul)
# ------------------------------------------------------------
p_dd_liq <- ggplot(coefs_dd_liq, aes(x = horizon)) +
  # DD × wide (verde)
  geom_line(aes(y = beta_dd_wide), size = 1, color = "darkgreen") +
  geom_ribbon(aes(ymin = beta_dd_wide - 1.96 * se_dd_wide,
                  ymax = beta_dd_wide + 1.96 * se_dd_wide),
              fill = "darkgreen", alpha = 0.2) +
  # Liquidez × wide (azul)
  geom_line(aes(y = beta_liq_wide), size = 1, color = "blue") +
  geom_ribbon(aes(ymin = beta_liq_wide - 1.96 * se_liq_wide,
                  ymax = beta_liq_wide + 1.96 * se_liq_wide),
              fill = "blue", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 16(b): Joint Dynamics – DD (verde) y Liquidez (azul)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_dd_liq)

# ------------------------------------------------------------
# Fin del script para Figura 16 (Leverage/DD y Liquidez)
# ------------------------------------------------------------






# *******************************************************************************
# Figura 21: Dinámicas de las Respuestas Diferenciales a Choques Monetarios por Volatilidad
#*******************************************************************************

# *******************************************************************************
# Figura 21: Dinámica de la Respuesta Diferencial a Choques Monetarios por Volatilidad (5‐Year ÚNICAMENTE)
# (Usando TODO el panel sin recortes de periodos antiguos)
#*******************************************************************************

### DEBE SER A 5 AÑOS DE VOLATILIDAD Y 10, SI NO SE PUEDE A 10, A 7.


# 0) Cargar librerías necesarias
#    (Instalar si hiciera falta: haven, zoo, dplyr, tidyr, purrr, fixest, ggplot2)
library(haven)     # read_dta()
library(zoo)       # as.yearqtr(), rollapply()
library(dplyr)     # Manipulación de datos
library(tidyr)     # replace_na()
library(purrr)     # map(), map_dfc()
library(fixest)    # feols()
library(ggplot2)   # Graficar

# --------------------------------------------------------------------------------
# 1) Leer el panel “Data_Base_final.dta”
#    (asegúrate de que contiene: current_ratio, saleq, leverage, dd, dlog_gdp, mpso, mps)
# --------------------------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq = as.yearqtr(dateq)
  ) %>%
  arrange(name, dateq)

# --------------------------------------------------------------------------------
# 2) Preprocesar variables firm‐level:
#    - dlog_capital y Ldl_capital
#    - rsales_g_std  (crecimiento de ventas estandarizado)
#    - lev_std, dd_std  (apalancamiento y DD estandarizados)
#    - liq_std  (liquidez estandarizada a partir de current_ratio)
#    - Crear “yoy_sales” = Δ log(saleq) año‐año
#    - Calcular volatilidad de ventas 5‐Year:
#          • vol_5yr  = sd de yoy_sales en ventana móvil de 20 trimestres (5 años)
#      Luego estandarizar vol_5yr dentro de cada firma.
# --------------------------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) Evitar log(0) en “capital”
    capital2     = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    log_capital  = log(capital2),
    lag_log_cap  = lag(log_capital),
    dlog_capital = log_capital - lag_log_cap,
    Ldl_capital  = lag(dlog_capital, 1),
    
    # 2.2) Crecimiento de ventas para rsales_g_std
    saleq2       = if_else(!is.na(saleq) & saleq > 0, saleq, NA_real_),
    log_sales    = log(saleq2),
    rsales_g     = log_sales - lag(log_sales),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.3) Estandarizar apalancamiento y DD
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.4) Estandarizar liquidez a partir de current_ratio
    liq_std = if_else(
      !is.na(current_ratio),
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.5) Preparar “yoy_sales” = Δ log(saleq) año‐año
    yoy_sales = if_else(
      !is.na(log_sales) & !is.na(lag(log_sales, 4)),
      log_sales - lag(log_sales, 4),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.6) Volatilidad 5‐Year: sd de yoy_sales en ventana móvil de 20 trimestres
    vol_5yr = zoo::rollapply(
      data  = yoy_sales,
      width = 20,
      FUN   = function(x) sd(x, na.rm = TRUE),
      align = "right",
      fill  = NA
    ),
    
    # 2.7) Estandarizar vol_5yr dentro de cada firma
    vol_5yr_std = (vol_5yr - mean(vol_5yr, na.rm = TRUE)) / sd(vol_5yr, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # 2.8) Eliminar variables intermedias que ya no necesitamos
  select(-capital2, -lag_log_cap, -saleq2, -log_sales, -rsales_g, -yoy_sales)

# --------------------------------------------------------------------------------
# 3) Crear interacciones de volatilidad 5‐Year con choques monetarios:
#      • vol_5yr_std_gdp  = vol_5yr_std × dlog_gdp
#      • vol_5yr_std_wide = vol_5yr_std × mpso
#      • vol_5yr_std_mps  = vol_5yr_std × mps
# --------------------------------------------------------------------------------
df <- df %>%
  mutate(
    vol_5yr_std_gdp  = vol_5yr_std * dlog_gdp,
    vol_5yr_std_wide = vol_5yr_std * mpso,
    vol_5yr_std_mps  = vol_5yr_std * mps
  )

# --------------------------------------------------------------------------------
# 4) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
#    (local projections dentro de cada firma)
# --------------------------------------------------------------------------------
df_dyn <- df %>% arrange(name, dateq)
for (h in 0:12) {
  new_col <- paste0("cumF", h, "_dlog_capital")
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      !!new_col := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# --------------------------------------------------------------------------------
# 5) Reemplazar NA en regresores de RHS por 0
#    (para conservar filas a menos que falte cumFh_dlog_capital)
# --------------------------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        "rsales_g_std", "Ldl_capital",
        "vol_5yr_std_gdp", "vol_5yr_std_wide", "vol_5yr_std_mps"
      ),
      ~ replace_na(.x, 0)
    )
  )

# --------------------------------------------------------------------------------
# 6) Definir controles firma‐level (sin size):
#      – “rsales_g_std” y “liq_std”
#    Efectos fijos: name + Country_x
# --------------------------------------------------------------------------------
controls_firm_nosize <- c("rsales_g_std", "liq_std")

# --------------------------------------------------------------------------------
# 7) Dinámica 5‐Year: bucle h = 0…12
#      Especificación:
#        cumFh_dlog_capital ~
#          vol_5yr_std_gdp + vol_5yr_std_wide + vol_5yr_std_mps +
#          rsales_g_std + liq_std
#        | name + Country_x, cluster(~dateq + name)
# --------------------------------------------------------------------------------
models_vol_5yr <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "vol_5yr_std_gdp + vol_5yr_std_wide + vol_5yr_std_mps + ",
    paste(controls_firm_nosize, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_vol_5yr) <- paste0("m", 0:12)

# --------------------------------------------------------------------------------
# 8) Extraer coeficiente y error estándar de “vol_5yr_std_wide”
# --------------------------------------------------------------------------------
coefs_vol_5yr <- tibble(
  horizon           = 0:12,
  beta_vol_5yr_wide = map_dbl(models_vol_5yr, function(mod) {
    coefs <- coef(mod)
    if ("vol_5yr_std_wide" %in% names(coefs)) {
      coefs["vol_5yr_std_wide"]
    } else {
      0
    }
  }),
  se_vol_5yr_wide = map_dbl(models_vol_5yr, function(mod) {
    V <- try(vcov(mod), silent = TRUE)
    if (!inherits(V, "try-error") &&
        "vol_5yr_std_wide" %in% rownames(V) &&
        "vol_5yr_std_wide" %in% colnames(V)) {
      sqrt(V["vol_5yr_std_wide", "vol_5yr_std_wide"])
    } else {
      0
    }
  })
)

cat("\nFigura 21(a): Dinámica – Volatilidad 5-Year × Shock (mpso)\n")
print(coefs_vol_5yr)

# --------------------------------------------------------------------------------
# 9) Graficar coeficientes “vol_5yr_std_wide” con bandas ±1.96 × SE
# --------------------------------------------------------------------------------
p_vol_5yr <- ggplot(coefs_vol_5yr, aes(x = horizon, y = beta_vol_5yr_wide)) +
  geom_line(size = 1, color = "purple") +
  geom_point(size = 2, color = "purple") +
  geom_ribbon(aes(ymin = beta_vol_5yr_wide - 1.96 * se_vol_5yr_wide,
                  ymax = beta_vol_5yr_wide + 1.96 * se_vol_5yr_wide),
              fill = "purple", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title    = "Figura 21(a): Dinámica – Volatilidad 5-Year × Shock (mpso)",
    subtitle = "Horizonte 0–12 trimestres",
    x        = "Horizonte (trimestres)",
    y        = "Coeficiente acumulado"
  ) +
  theme_minimal()
print(p_vol_5yr)

# --------------------------------------------------------------------------------
# 10) Fin del script para Figura 21 (solo Volatilidad 5‐Year)
# --------------------------------------------------------------------------------







# *******************************************************************************
# Figura 22: Dinámica Conjunta de Posición Financiera (DD / Leverage) y Volatilidad 5‐Year
# (Usando TODO el panel sin recortes de periodos antiguos; horizonte h = 0…12)
#*******************************************************************************

# 0) Cargar librerías necesarias
library(haven)     # read_dta()
library(zoo)       # as.yearqtr(), rollapply()
library(dplyr)     # Manipulación de datos
library(tidyr)     # replace_na()
library(purrr)     # map(), map_dfc()
library(fixest)    # feols()
library(ggplot2)   # Graficar (si se quisiera)

# --------------------------------------------------------------------------------
# 1) Leer el panel “Data_Base_final.dta” y ordenar por firma + fecha
# --------------------------------------------------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(dateq = as.yearqtr(dateq)) %>%
  arrange(name, dateq)

# --------------------------------------------------------------------------------
# 2) Preprocesar variables firm‐level (idéntico a Figura 21):
#    - dlog_capital, Ldl_capital
#    - rsales_g_std  (crecimiento de ventas estandarizado)
#    - lev_std, dd_std  (apalancamiento y DD estandarizados)
#    - liq_std  (liquidez estandarizada a partir de current_ratio)
#    - “yoy_sales” = Δ log(saleq) año‐año
#    - vol_5yr = sd de yoy_sales en ventana móvil de 20 trimestres (5 años), estandarizada dentro de cada firma
#    - Interacciones:
#        • vol_5yr_std_gdp   = vol_5yr_std × dlog_gdp
#        • vol_5yr_std_wide  = vol_5yr_std × mpso
#        • vol_5yr_std_mps   = vol_5yr_std × mps
#        • d2d_wins_dem_std_gdp  = dd_std × dlog_gdp
#        • d2d_wins_dem_std_wide = dd_std × mpso
#        • d2d_wins_dem_std      = dd_std × mps
#        • lev_wins_dem_std_gdp   = lev_std × dlog_gdp
#        • lev_wins_dem_std_wide  = lev_std × mpso
#        • lev_wins_dem_std       = lev_std × mps
# --------------------------------------------------------------------------------
df <- df %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.1) dlog_capital y su rezago
    capital2     = if_else(!is.na(capital) & capital > 0, capital, NA_real_),
    log_capital  = log(capital2),
    lag_log_cap  = lag(log_capital),
    dlog_capital = log_capital - lag_log_cap,
    Ldl_capital  = lag(dlog_capital, 1),
    
    # 2.2) Crecimiento de ventas para rsales_g_std
    saleq2       = if_else(!is.na(saleq) & saleq > 0, saleq, NA_real_),
    log_sales    = log(saleq2),
    rsales_g     = log_sales - lag(log_sales),
    rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE),
    
    # 2.3) Estandarizar apalancamiento y DD
    lev_std = if_else(
      !is.na(leverage),
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE),
      NA_real_
    ),
    dd_std  = if_else(
      !is.na(dd),
      (dd - mean(dd, na.rm = TRUE)) / sd(dd, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.4) Estandarizar liquidez a partir de current_ratio
    liq_std = if_else(
      !is.na(current_ratio),
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE),
      NA_real_
    ),
    
    # 2.5) Preparar “yoy_sales” = Δ log(saleq) año‐año
    yoy_sales = if_else(
      !is.na(log_sales) & !is.na(lag(log_sales, 4)),
      log_sales - lag(log_sales, 4),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  group_by(name) %>%
  arrange(dateq) %>%
  mutate(
    # 2.6) Volatilidad 5‐Year: sd de yoy_sales en ventana móvil de 20 trimestres
    vol_5yr = zoo::rollapply(
      data  = yoy_sales,
      width = 20,
      FUN   = function(x) sd(x, na.rm = TRUE),
      align = "right",
      fill  = NA
    ),
    
    # 2.7) Estandarizar vol_5yr dentro de cada firma
    vol_5yr_std = (vol_5yr - mean(vol_5yr, na.rm = TRUE)) / sd(vol_5yr, na.rm = TRUE),
    
    # 2.8) Interacciones para volatilidad 5‐Year con choques monetarios
    vol_5yr_std_gdp  = vol_5yr_std * dlog_gdp,
    vol_5yr_std_wide = vol_5yr_std * mpso,
    vol_5yr_std_mps  = vol_5yr_std * mps,
    
    # 2.9) Interacciones para DD con choques monetarios
    d2d_wins_dem_std_gdp  = dd_std * dlog_gdp,
    d2d_wins_dem_std_wide = dd_std * mpso,
    d2d_wins_dem_std      = dd_std * mps,
    
    # 2.10) Interacciones para Leverage con choques monetarios
    lev_wins_dem_std_gdp  = lev_std * dlog_gdp,
    lev_wins_dem_std_wide = lev_std * mpso,
    lev_wins_dem_std      = lev_std * mps
  ) %>%
  ungroup() %>%
  # 2.11) Eliminar variables intermedias ya no necesarias
  select(-capital2, -lag_log_cap, -saleq2, -log_sales, -rsales_g, -yoy_sales)

# --------------------------------------------------------------------------------
# 3) Generar variables dinámicas acumuladas cumF0…cumF12 de dlog_capital
#    (local projections dentro de cada firma)
# --------------------------------------------------------------------------------
df_dyn <- df %>% arrange(name, dateq)
for (h in 0:12) {
  new_col <- paste0("cumF", h, "_dlog_capital")
  df_dyn <- df_dyn %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      !!new_col := rowSums(
        map_dfc(0:h, ~ lead(dlog_capital, .x)),
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# --------------------------------------------------------------------------------
# 4) Reemplazar NA en regresores de RHS por 0
#    (para conservar filas salvo que falte cumFh_dlog_capital)
# --------------------------------------------------------------------------------
df_dyn <- df_dyn %>%
  mutate(
    across(
      c(
        "rsales_g_std", "liq_std", "Ldl_capital",
        # Interacciones volatilidad 5‐Year
        "vol_5yr_std_gdp", "vol_5yr_std_wide", "vol_5yr_std_mps",
        # Interacciones DD
        "d2d_wins_dem_std_gdp", "d2d_wins_dem_std_wide", "d2d_wins_dem_std",
        # Interacciones Leverage
        "lev_wins_dem_std_gdp", "lev_wins_dem_std_wide", "lev_wins_dem_std"
      ),
      ~ replace_na(.x, 0)
    )
  )

# --------------------------------------------------------------------------------
# 5) Definir controles firma‐level:
#      – “rsales_g_std” y “liq_std”
#    Efectos fijos: name + Country_x
# --------------------------------------------------------------------------------
controls_firm <- c("rsales_g_std", "liq_std")

# --------------------------------------------------------------------------------
# 6) (a) Dinámica conjunta: Distance‐to‐Default y Volatilidad 5‐Year
#      Bucle h = 0…12:
#        cumFh_dlog_capital ~
#          vol_5yr_std_gdp + vol_5yr_std_wide + vol_5yr_std_mps +
#          d2d_wins_dem_std_gdp + d2d_wins_dem_std_wide + d2d_wins_dem_std +
#          rsales_g_std + liq_std
#        | name + Country_x, cluster(~dateq + name)
# --------------------------------------------------------------------------------
models_dd_volatility <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "vol_5yr_std_gdp + vol_5yr_std_wide + vol_5yr_std_mps + ",
    "d2d_wins_dem_std_gdp + d2d_wins_dem_std_wide + d2d_wins_dem_std + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_dd_volatility) <- paste0("m", 0:12)

# --------------------------------------------------------------------------------
# 7) Extraer coeficientes y errores de “vol_5yr_std_wide” y “d2d_wins_dem_std_wide”
# --------------------------------------------------------------------------------
coefs_dd_volatility <- tibble(
  horizon                  = 0:12,
  beta_vol_5yr_wide_dd     = map_dbl(models_dd_volatility, function(mod) {
    cfs <- coef(mod)
    if ("vol_5yr_std_wide" %in% names(cfs)) cfs["vol_5yr_std_wide"] else NA_real_
  }),
  se_vol_5yr_wide_dd       = map_dbl(models_dd_volatility, function(mod) {
    V <- try(vcov(mod), silent = TRUE)
    if (!inherits(V, "try-error") &&
        "vol_5yr_std_wide" %in% rownames(V) &&
        "vol_5yr_std_wide" %in% colnames(V)) {
      sqrt(V["vol_5yr_std_wide", "vol_5yr_std_wide"])
    } else {
      NA_real_
    }
  }),
  beta_d2d_wide            = map_dbl(models_dd_volatility, function(mod) {
    cfs <- coef(mod)
    if ("d2d_wins_dem_std_wide" %in% names(cfs)) cfs["d2d_wins_dem_std_wide"] else NA_real_
  }),
  se_d2d_wide              = map_dbl(models_dd_volatility, function(mod) {
    V <- try(vcov(mod), silent = TRUE)
    if (!inherits(V, "try-error") &&
        "d2d_wins_dem_std_wide" %in% rownames(V) &&
        "d2d_wins_dem_std_wide" %in% colnames(V)) {
      sqrt(V["d2d_wins_dem_std_wide", "d2d_wins_dem_std_wide"])
    } else {
      NA_real_
    }
  })
)

cat("\nFigura 22(a): Dinámica conjunta – DD × Volatilidad 5-Year\n")
print(coefs_dd_volatility)

#---------------------------------------------------------------------------------
# 8) (b) Dinámica conjunta: Leverage y Volatilidad 5‐Year
#      Bucle h = 0…12:
#        cumFh_dlog_capital ~
#          vol_5yr_std_gdp + vol_5yr_std_wide + vol_5yr_std_mps +
#          lev_wins_dem_std_gdp + lev_wins_dem_std_wide + lev_wins_dem_std +
#          rsales_g_std + liq_std
#        | name + Country_x, cluster(~dateq + name)
# --------------------------------------------------------------------------------
models_lev_volatility <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var, " ~ ",
    "vol_5yr_std_gdp + vol_5yr_std_wide + vol_5yr_std_mps + ",
    "lev_wins_dem_std_gdp + lev_wins_dem_std_wide + lev_wins_dem_std + ",
    paste(controls_firm, collapse = " + "),
    " | name + Country_x"
  ))
  feols(fml, data = df_dyn, cluster = ~dateq + name)
})
names(models_lev_volatility) <- paste0("m", 0:12)

# --------------------------------------------------------------------------------
# 9) Extraer coeficientes y errores de “vol_5yr_std_wide” y “lev_wins_dem_std_wide”
# --------------------------------------------------------------------------------
coefs_lev_volatility <- tibble(
  horizon                  = 0:12,
  beta_vol_5yr_wide_lev    = map_dbl(models_lev_volatility, function(mod) {
    cfs <- coef(mod)
    if ("vol_5yr_std_wide" %in% names(cfs)) cfs["vol_5yr_std_wide"] else NA_real_
  }),
  se_vol_5yr_wide_lev      = map_dbl(models_lev_volatility, function(mod) {
    V <- try(vcov(mod), silent = TRUE)
    if (!inherits(V, "try-error") &&
        "vol_5yr_std_wide" %in% rownames(V) &&
        "vol_5yr_std_wide" %in% colnames(V)) {
      sqrt(V["vol_5yr_std_wide", "vol_5yr_std_wide"])
    } else {
      NA_real_
    }
  }),
  beta_lev_wide            = map_dbl(models_lev_volatility, function(mod) {
    cfs <- coef(mod)
    if ("lev_wins_dem_std_wide" %in% names(cfs)) cfs["lev_wins_dem_std_wide"] else NA_real_
  }),
  se_lev_wide              = map_dbl(models_lev_volatility, function(mod) {
    V <- try(vcov(mod), silent = TRUE)
    if (!inherits(V, "try-error") &&
        "lev_wins_dem_std_wide" %in% rownames(V) &&
        "lev_wins_dem_std_wide" %in% colnames(V)) {
      sqrt(V["lev_wins_dem_std_wide", "lev_wins_dem_std_wide"])
    } else {
      NA_real_
    }
  })
)

cat("\nFigura 22(b): Dinámica conjunta – Leverage × Volatilidad 5-Year\n")
print(coefs_lev_volatility)

# --------------------------------------------------------------------------------
# 10) (Opcional) Graficar cada serie con bandas ±1.96 × SE
#     Ejemplo para Figura 22(a):
# --------------------------------------------------------------------------------
p_dd_vol <- ggplot(coefs_dd_volatility, aes(x = horizon)) +
  geom_line(aes(y = beta_vol_5yr_wide_dd), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = beta_vol_5yr_wide_dd - 1.96 * se_vol_5yr_wide_dd,
                  ymax = beta_vol_5yr_wide_dd + 1.96 * se_vol_5yr_wide_dd),
              fill = "blue", alpha = 0.2) +
  geom_line(aes(y = beta_d2d_wide), color = "red", size = 1) +
  geom_ribbon(aes(ymin = beta_d2d_wide - 1.96 * se_d2d_wide,
                  ymax = beta_d2d_wide + 1.96 * se_d2d_wide),
              fill = "red", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 22(a): Dinámica conjunta – DD (rojo) vs Volatilidad 5-Year (azul)",
    subtitle = "Horizonte 0–12 trimestres",
    x = "Horizonte (trimestres)",
    y = "Coeficiente acumulado"
  ) +
  theme_minimal()
p_dd_vol  # Mostrar gráfico para 22(a)

# Ejemplo para Figura 22(b):
p_lev_vol <- ggplot(coefs_lev_volatility, aes(x = horizon)) +
  geom_line(aes(y = beta_vol_5yr_wide_lev), color = "purple", size = 1) +
  geom_ribbon(aes(ymin = beta_vol_5yr_wide_lev - 1.96 * se_vol_5yr_wide_lev,
                  ymax = beta_vol_5yr_wide_lev + 1.96 * se_vol_5yr_wide_lev),
              fill = "purple", alpha = 0.2) +
  geom_line(aes(y = beta_lev_wide), color = "darkgreen", size = 1) +
  geom_ribbon(aes(ymin = beta_lev_wide - 1.96 * se_lev_wide,
                  ymax = beta_lev_wide + 1.96 * se_lev_wide),
              fill = "darkgreen", alpha = 0.2) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 22(b): Dinámica conjunta – Leverage (verde) vs Volatilidad 5-Year (púrpura)",
    subtitle = "Horizonte 0–12 trimestres",
    x = "Horizonte (trimestres)",
    y = "Coeficiente acumulado"
  ) +
  theme_minimal()
p_lev_vol  # Mostrar gráfico para 22(b)

# --------------------------------------------------------------------------------
# 11) Fin del script para Figura 22
# --------------------------------------------------------------------------------










