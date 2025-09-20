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
install.packages("ggplot2")

# ===================================================
# Script en R para replicar “descriptive_statistics.do”
# Ottonello & Winberry (2020) adaptado a tu base
# ===================================================

# 0) Cargar librerías
library(haven)    # read_dta()
library(dplyr)    # manipulación
library(zoo)      # as.yearqtr()
library(fixest)   # feols(), resid()
library(tidyr)    # pivot_longer()
library(ggplot2)  # gráficos

# Funciones auxiliares -----------------------------------------------------
winsorize_vec <- function(x, probs = c(0.005, 0.995)) {
  if (!is.numeric(x)) {
    return(x)
  }
  if (all(is.na(x))) {
    return(x)
  }
  qs <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
  lower <- qs[1]
  upper <- qs[2]
  pmax(pmin(x, upper), lower)
}

standardize_vec <- function(x) {
  if (!is.numeric(x)) {
    return(rep(NA_real_, length(x)))
  }
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  mu <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mu) / sd_val
}

mean_if_any <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

calc_distribution <- function(data, var,
                              probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)) {
  x <- data[[var]]
  result <- data.frame(Variable = var, stringsAsFactors = FALSE)
  if (is.null(x) || !is.numeric(x)) {
    result$Mean <- NA_real_
    result$SD <- NA_real_
    for (p in probs) {
      col_name <- paste0("P", as.integer(round(p * 100)))
      result[[col_name]] <- NA_real_
    }
    result$N <- 0L
    return(result)
  }
  non_missing <- x[!is.na(x)]
  result$Mean <- if (length(non_missing) == 0) NA_real_ else mean(non_missing)
  result$SD <- if (length(non_missing) <= 1) NA_real_ else sd(non_missing)
  for (p in probs) {
    col_name <- paste0("P", as.integer(round(p * 100)))
    result[[col_name]] <- if (length(non_missing) == 0) {
      NA_real_
    } else {
      quantile(non_missing, probs = p, names = FALSE, type = 7)
    }
  }
  result$N <- length(non_missing)
  result
}

# 1) Directorio base y carga de datos --------------------------------------
base_dir <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final"
df <- read_dta(file.path(base_dir, "Data_Base_final.dta")) %>%
  mutate(
    dateq = as.yearqtr(dateq)  # convierte tu variable trimestral
  )

# 1.1) Winsorización (0.5%) y variables derivadas --------------------------
vars_to_winsorize <- c(
  "mps", "mpso", "dlog_capital", "leverage", "dd",
  "dln_saleq", "log_size", "current_ratio", "atq", "saleq"
)
vars_to_winsorize <- intersect(vars_to_winsorize, names(df))

if (length(vars_to_winsorize) > 0) {
  df <- df %>%
    mutate(across(all_of(vars_to_winsorize), winsorize_vec))
}

if ("Country_x" %in% names(df) && !("country" %in% names(df))) {
  df <- df %>% mutate(country = Country_x)
} else if (!("country" %in% names(df)) && "Country_y" %in% names(df)) {
  df <- df %>% mutate(country = Country_y)
}

sector_vars <- grep("^sec_", names(df), value = TRUE)
if (length(sector_vars) > 0) {
  sector_mat <- df %>% select(all_of(sector_vars)) %>% as.matrix()
  sector_mat_na <- sector_mat
  sector_mat_na[is.na(sector_mat_na)] <- -Inf
  sector_idx <- max.col(sector_mat_na, ties.method = "first")
  sector_mat_zero <- sector_mat
  sector_mat_zero[is.na(sector_mat_zero)] <- 0
  has_sector <- rowSums(sector_mat_zero) > 0
  sector_clean <- gsub("^sec_", "", sector_vars)
  sector_vec <- rep(NA_character_, nrow(df))
  sector_vec[has_sector] <- sector_clean[sector_idx[has_sector]]
  df$sector <- sector_vec
}

if (all(c("atq", "saleq") %in% names(df))) {
  df <- df %>%
    mutate(
      size_raw = if_else(atq > 0 & saleq > 0,
                         (log(atq) + log(saleq)) / 2,
                         NA_real_)
    ) %>%
    mutate(size_raw = winsorize_vec(size_raw))
} else {
  df$size_raw <- NA_real_
}

if ("leverage" %in% names(df)) {
  df <- df %>% mutate(leverage_std = standardize_vec(leverage))
}
if ("dd" %in% names(df)) {
  df <- df %>% mutate(dd_std = standardize_vec(dd))
}
# -------------------------------
# Tabla 2: Firm-Level Variables
# -------------------------------

panel <- df %>% select(
  name,          # Identificador (ticker) :contentReference[oaicite:1]{index=1}
  dateq,         # Fecha trimestral :contentReference[oaicite:2]{index=2}
  dlog_capital,  # Δ log capital :contentReference[oaicite:3]{index=3}
  leverage,      # Apalancamiento = tdebt/atq :contentReference[oaicite:4]{index=4}
  dd             # Distance to default :contentReference[oaicite:5]{index=5}
)

# Panel A: distribuciones marginales de dlog_capital, leverage y dd (winsorizadas)
varsA <- c("dlog_capital", "leverage", "dd")
tabA <- panel %>%
  summarise(across(all_of(varsA),
                   list(
                     Mean   = ~mean_if_any(.x),
                     Median = ~if (all(is.na(.x))) NA_real_ else median(.x, na.rm = TRUE),
                     SD     = ~if (sum(!is.na(.x)) <= 1) NA_real_ else sd(.x, na.rm = TRUE),
                     P95    = ~if (all(is.na(.x))) NA_real_ else quantile(.x, 0.95, na.rm = TRUE, names = FALSE),
                     N      = ~sum(!is.na(.x))
                   ),
                   .names = "{.col}_{.fn}"))

tabA_print <- do.call(rbind, lapply(varsA, function(v) {
  data.frame(
    Variable      = v,
    Mean          = tabA[[paste0(v, "_Mean")]],
    Median        = tabA[[paste0(v, "_Median")]],
    `S.D.`        = tabA[[paste0(v, "_SD")]],
    `95th Pct.`   = tabA[[paste0(v, "_P95")]],
    Observaciones = tabA[[paste0(v, "_N")]]
  )
}))
cat("\n--- Tabla 2A: Distribuciones marginales (winsorizadas) ---\n")
print(tabA_print)

# Panel B: matriz de correlaciones (winsorizada)
corrB <- panel %>%
  select(all_of(varsA)) %>%
  cor(use = "pairwise.complete.obs")
cat("\n--- Tabla 2B: Correlaciones (winsorizadas) ---\n")
print(round(corrB, 3))

# -------------------------------
# Panel C: correlaciones residualizadas (solo FE por firma)
# -------------------------------

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

panel_resid <- panel_fc %>%
  mutate(
    leverage_resid = resid(
      feols(leverage ~ dln_saleq + log_size + current_ratio | name, data = panel_fc)
    ),
    dd_resid = resid(
      feols(dd ~ dln_saleq + log_size + current_ratio | name, data = panel_fc)
    )
  )

corrC <- panel_resid %>%
  select(leverage_resid, dd_resid) %>%
  cor(use = "pairwise.complete.obs")

cat("\n--- Tabla 2C: Correlaciones (residualizadas, FE firma, winsorizadas) ---\n")
print(round(corrC, 3))

# ---------------------------------------------------
# 3) Conteo de firmas y observaciones por país y sector
# ---------------------------------------------------

cat("\n--- Conteo por país (winsorizado) ---\n")
conteo_pais <- df %>%
  filter(!is.na(country)) %>%
  group_by(country) %>%
  summarise(
    Firms = n_distinct(name),
    Observaciones = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Observaciones))
print(conteo_pais)

cat("\n--- Conteo por sector (winsorizado) ---\n")
conteo_sector <- df %>%
  filter(!is.na(sector)) %>%
  group_by(sector) %>%
  summarise(
    Firms = n_distinct(name),
    Observaciones = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Observaciones))
print(conteo_sector)

cat("\n--- Principales combinaciones país-sector (top 10) ---\n")
conteo_pais_sector <- df %>%
  filter(!is.na(country), !is.na(sector)) %>%
  group_by(country, sector) %>%
  summarise(
    Firms = n_distinct(name),
    Observaciones = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Observaciones)) %>%
  slice_head(n = 10)
print(conteo_pais_sector)

# ---------------------------------------------------
# 4) Top firmas por tamaño (size_raw winsorizado)
# ---------------------------------------------------

if ("size_raw" %in% names(df)) {
  top_firmas_tamano <- df %>%
    filter(!is.na(size_raw)) %>%
    group_by(name, country, sector) %>%
    summarise(
      Promedio_size = mean(size_raw, na.rm = TRUE),
      Observaciones = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Promedio_size)) %>%
    slice_head(n = 15)
  cat("\n--- Top firmas por tamaño (size_raw) ---\n")
  print(top_firmas_tamano)
}

# ---------------------------------------------------
# 5) Distribuciones y percentiles de leverage y dd
# ---------------------------------------------------

distribution_vars <- intersect(c("leverage", "dd"), names(df))
if (length(distribution_vars) > 0) {
  distribucion_tabla <- do.call(
    rbind,
    lapply(distribution_vars, function(v) calc_distribution(df, v))
  )
  cat("\n--- Tabla de percentiles (winsorizados) ---\n")
  print(distribucion_tabla)
}

std_cols <- intersect(c("leverage_std", "dd_std"), names(df))
if (length(std_cols) > 0) {
  dist_long <- df %>%
    select(all_of(std_cols)) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "valor") %>%
    drop_na(valor) %>%
    mutate(variable = recode(variable,
                             leverage_std = "Leverage (std)",
                             dd_std = "Distance to Default (std)"))
} else {
  dist_long <- data.frame(variable = character(), valor = numeric(),
                          stringsAsFactors = FALSE)
}

if (nrow(dist_long) > 0) {
  hist_plot <- ggplot(dist_long, aes(x = valor, fill = variable)) +
    geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
    facet_wrap(~variable, scales = "free_y") +
    labs(
      title = "Distribuciones estandarizadas de leverage y dd",
      x = "Valor estandarizado",
      y = "Frecuencia"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  print(hist_plot)

  violin_plot <- ggplot(dist_long, aes(x = variable, y = valor, fill = variable)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.6) +
    labs(
      title = "Distribuciones estandarizadas de leverage y dd",
      x = "",
      y = "Valor estandarizado"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  print(violin_plot)
} else {
  cat("\nNo hay datos suficientes para graficar las distribuciones estandarizadas.\n")
}

# ---------------------------------------------------
# 6) Matriz de correlaciones entre leverage, dd y tamaño
# ---------------------------------------------------

corr_vars <- intersect(c("leverage", "dd", "size_raw"), names(df))
if (length(corr_vars) >= 2) {
  corr_matrix <- df %>%
    select(all_of(corr_vars)) %>%
    cor(use = "pairwise.complete.obs")
  cat("\n--- Matriz de correlaciones (winsorizadas) ---\n")
  print(round(corr_matrix, 3))
}

# ---------------------------------------------------
# 7) Evolución temporal del shock
# ---------------------------------------------------

shock_ts <- df1 %>%
  group_by(dateq) %>%
  summarise(
    HF_shock = mean_if_any(HF_shock),
    Q_shock = mean_if_any(Q_shock),
    .groups = "drop"
  ) %>%
  arrange(dateq)

cat("\n--- Evolución temporal promedio de los shocks (winsorizados) ---\n")
print(shock_ts)

if (nrow(shock_ts) > 0) {
  shock_ts_long <- shock_ts %>%
    pivot_longer(cols = c(HF_shock, Q_shock),
                 names_to = "Shock",
                 values_to = "Valor")

  shock_plot <- ggplot(shock_ts_long, aes(x = dateq, y = Valor, color = Shock)) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    labs(
      title = "Evolución temporal de los shocks monetarios",
      x = "Trimestre",
      y = "Shock (winsorizado, %)",
      color = "Serie"
    ) +
    theme_minimal()
  print(shock_plot)
} else {
  cat("\nNo hay observaciones para graficar la evolución temporal de los shocks.\n")
}