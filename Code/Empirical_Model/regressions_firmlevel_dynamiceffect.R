# --------------------------------------------------------
# 0) Cargar librerías necesarias
# --------------------------------------------------------
library(haven)
library(zoo)
library(dplyr)
library(purrr)
library(fixest)
library(readr)
library(knitr)
library(kableExtra)
library(ggplot2)
library(patchwork)
library(tidyr)

# --------------------------------------------------------
# 1) Leer base de datos y convertir fecha a trimestral
# --------------------------------------------------------
ruta_dta <- "C:/Users/joser/Downloads/Tesis Master/Data Tesis/Data Final/Data_Base_final.dta"
base_dir <- dirname(ruta_dta)

df <- read_dta(ruta_dta) %>%
  mutate(dateq = as.yearqtr(dateq)) %>%
  arrange(name, dateq)


# ================================================
# 4) Preprocesamiento: dlog_capital
# ================================================
df <- df %>%
  arrange(name, dateq) %>%
  group_by(name) %>%
  mutate(
    ratio_cap     = real_capital / lag(real_capital),
    dlog_capital  = suppressWarnings(ifelse(ratio_cap > 0, log(ratio_cap), NA_real_)),
    shock = -shock
  ) %>%
  select(-ratio_cap) %>%
  ungroup()

#--------------------------------------------------------
# Figure 6a: Respuestas dinámicas de los pagos de interés (horizonte 0–12) con controles macro
# --------------------------------------------------------

# --------------------------------------------------------
# 2) Crear variable int_exp (pagos de interés) si no existe
# --------------------------------------------------------
if (!"int_exp" %in% names(df)) {
  df <- df %>% mutate(int_exp = oiadpq)
} else {
  message("Variable 'int_exp' ya existe; omitiendo creación.")
}

# --------------------------------------------------------
# 3) Controles firm-level y choque puro (verificar existencia)
# --------------------------------------------------------
if (!"rsales_g_std" %in% names(df)) {
  df <- df %>%
    group_by(name) %>% arrange(dateq) %>%
    mutate(rsales_g_std = {
      x <- log(saleq) - log(lag(saleq))
      (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
    }) %>% ungroup()
} else {
  message("Variable 'rsales_g_std' ya existe; omitiendo.")
}

if (!"size_index" %in% names(df)) {
  df <- df %>%
    group_by(name) %>% arrange(dateq) %>%
    mutate(size_index = {
      tmp <- (
        (log(atq) - mean(log(atq), na.rm=TRUE)) / sd(log(atq), na.rm=TRUE) +
          (saleq - mean(saleq, na.rm=TRUE)) / sd(saleq, na.rm=TRUE)
      ) / 2
      (tmp - mean(tmp, na.rm=TRUE)) / sd(tmp, na.rm=TRUE)
    }) %>% ungroup()
} else {
  message("Variable 'size_index' ya existe; omitiendo.")
}

if (!"sh_current_a_std" %in% names(df)) {
  df <- df %>%
    group_by(name) %>% arrange(dateq) %>%
    mutate(sh_current_a_std = {
      x <- current_ratio
      (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
    }) %>% ungroup()
} else {
  message("Variable 'sh_current_a_std' ya existe; omitiendo.")
}

# Estandarizar controles macro
macro_vars <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")
existing_macro <- intersect(macro_vars, names(df))
missing_std <- setdiff(paste0(existing_macro, "_std"), names(df))
if (length(existing_macro) > 0 && length(missing_std) > 0) {
  df <- df %>% mutate(
    across(all_of(existing_macro),
           ~ (.-mean(., na.rm=TRUE)) / sd(., na.rm=TRUE),
           .names = "{.col}_std")
  )
} else {
  message("Variables macro estandarizadas ya existen o no hay macros disponibles; omitiendo.")
}

if (!all(c("lev_std","dd_std") %in% names(df))) {
  df <- df %>%
    group_by(name) %>% arrange(dateq) %>%
    mutate(
      lev_std = {
        x <- leverage
        qs <- quantile(x, c(0.01, 0.99), na.rm=TRUE)
        x_w <- pmin(pmax(x, qs[1]), qs[2])
        x_dm <- x_w - mean(x_w, na.rm=TRUE)
        x_dm / sd(x_dm, na.rm=TRUE)
      },
      dd_std = {
        x <- dd
        qs <- quantile(x, c(0.01, 0.99), na.rm=TRUE)
        x_w <- pmin(pmax(x, qs[1]), qs[2])
        x_dm <- x_w - mean(x_w, na.rm=TRUE)
        x_dm / sd(x_dm, na.rm=TRUE)
      }
    ) %>% ungroup()
} else {
  message("Variables 'lev_std' y 'dd_std' ya existen; omitiendo.")
}

if (!all(c("lev_shock","d2d_shock") %in% names(df))) {
  df <- df %>% mutate(
    lev_shock = lev_std * shock,
    d2d_shock = dd_std  * shock
  )
} else {
  message("Variables 'lev_shock' y 'd2d_shock' ya existen; omitiendo.")
}

# --------------------------------------------------------
# 4) Construir variables dinámicas cumFh_int_exp para h = 0…12
# --------------------------------------------------------
dyn_vars <- paste0("cumF", 0:12, "_int_exp")
if (!all(dyn_vars %in% names(df))) {
  df_dyn <- df %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~ {
      tmp <- .x
      for (h in 0:12) {
        tmp[[paste0("cumF",h,"_int_exp")]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$int_exp, .x)), na.rm=TRUE)
      }
      tmp
    }) %>% ungroup()
} else {
  message("Variables dinámicas cumFh_int_exp ya existen; omitiendo.")
  df_dyn <- df
}


# --------------------------------------------------------
# 5) Definir controles firm-level y macro (con detección de faltantes)
# --------------------------------------------------------
controls_firm  <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro <- paste0(c("dlog_gdp", "dlog_cpi", "unemp", "embigl"), "_std")

# Detectar los controles macro efectivamente presentes en df_dyn
present_macro <- intersect(controls_macro, names(df_dyn))
if (length(present_macro) < length(controls_macro)) {
  missing_macro <- setdiff(controls_macro, present_macro)
  warning(
    "Faltan controles macro: ",
    paste(missing_macro, collapse = ", "),
    ". Se omitirán en la estimación."
  )
}

# Construir la cadena de controles sólo con los existentes
all_controls <- paste(c(controls_firm, present_macro), collapse = " + ")

# --------------------------------------------------------
# 6) Estimar dinámicas para h = 0…12 (solo coeficiente × choque puro)
# --------------------------------------------------------
library(purrr)
library(fixest)

res_IR_lev <- map(0:12, function(h) {
  fml <- as.formula(paste0(
    "cumF", h, "_int_exp ~ lev_shock + ", all_controls,
    " | name + Country"
  ))
  feols(fml, data = df_dyn, cluster = ~ name + Country)
})

res_IR_dd <- map(0:12, function(h) {
  fml <- as.formula(paste0(
    "cumF", h, "_int_exp ~ d2d_shock + ", all_controls,
    " | name + Country"
  ))
  feols(fml, data = df_dyn, cluster = ~ name + Country)
})

# --------------------------------------------------------
# 7) Extraer coeficientes y errores estándar
# --------------------------------------------------------
library(dplyr)

IR_lev <- tibble(
  horizon = 0:12,
  beta_lev = map_dbl(res_IR_lev, ~ coef(.x)["lev_shock"]),
  se_lev   = map_dbl(res_IR_lev, ~ sqrt(vcov(.x)["lev_shock", "lev_shock"]))
)

IR_dd <- tibble(
  horizon = 0:12,
  beta_dd  = map_dbl(res_IR_dd, ~ coef(.x)["d2d_shock"]),
  se_dd    = map_dbl(res_IR_dd, ~ sqrt(vcov(.x)["d2d_shock", "d2d_shock"]))
)

# --------------------------------------------------------
# 8) Combinar resultados y presentar tabla con kableExtra
# --------------------------------------------------------
library(kableExtra)
library(knitr)

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



# --------------------------------------------------------
# 9) Crear gráficos con ggplot2 y apilarlos usando patchwork
# --------------------------------------------------------
p_lev <- ggplot(IR_lev, aes(x=horizon, y=beta_lev)) +
  geom_line(size=0.8, color="steelblue") +
  geom_point(size=2, color="steelblue") +
  geom_ribbon(aes(ymin=beta_lev-1.96*se_lev, ymax=beta_lev+1.96*se_lev), alpha=0.2, fill="steelblue") +
  scale_x_continuous(breaks=0:12) +
  labs(title="Dinámica IR: Leverage × Shock", x="Horizonte (trimestres)", y="Δ pagos de interés acumulados") +
  theme_minimal()

p_dd <- ggplot(IR_dd, aes(x=horizon, y=beta_dd)) +
  geom_line(size=0.8, color="steelblue") +
  geom_point(size=2, color="steelblue") +
  geom_ribbon(aes(ymin=beta_dd-1.96*se_dd, ymax=beta_dd+1.96*se_dd), alpha=0.2, fill="steelblue") +
  scale_x_continuous(breaks=0:12) +
  labs(title="Dinámica IR: Distance-to-Default × Shock", x="Horizonte (trimestres)", y="Δ pagos de interés acumulados") +
  theme_minimal()

(p_lev / p_dd) + plot_annotation(title="Figure 6a: Respuestas dinámicas de pagos de interés (horizonte 0–12)")







# --------------------------------------------------------
# Figure 6b: Respuestas dinámicas del proxy de financiamiento externo (horizonte 0–12)
# --------------------------------------------------------

# Crear proxy financ_ext_proxy si no existe
if(!"new_debt"%in%names(df)){
  df<-df%>% group_by(name)%>%arrange(dateq)%>%
    mutate(new_debt=dlcq+dlttq)
} else message("Variable 'new_debt' ya existe; omitiendo.")
if(!"delta_new_debt"%in%names(df)){
  df<-df%>% group_by(name)%>%arrange(dateq)%>%
    mutate(delta_new_debt=new_debt-lag(new_debt))
} else message("Variable 'delta_new_debt' ya existe; omitiendo.")
if(!"financ_ext_proxy"%in%names(df)){
  df<-df%>% mutate(financ_ext_proxy=delta_new_debt/atq)
} else message("Variable 'financ_ext_proxy' ya existe; omitiendo.")

# Dinámicas proxy financiamiento externo
dyn_vars_fe<-paste0("cumF",0:12,"_financ_ext")
if(!all(dyn_vars_fe%in%names(df))){
  df_dyn_fe<-df%>%group_by(name)%>%arrange(dateq)%>%
    group_modify(~{tmp<-.x; for(h in 0:12) tmp[[paste0("cumF",h,"_financ_ext")]]<-
        rowSums(map_dfc(0:h,~lead(tmp$financ_ext_proxy,.x)),na.rm=TRUE); tmp})%>%ungroup()
} else {message("Variables dinámicas cumFh_financ_ext ya existen; omitiendo."); df_dyn_fe<-df}

# Estimar 6b: financiamiento externo
res_FE_lev<-map(0:12,function(h){feols(as.formula(paste0(
  "cumF",h,"_financ_ext~lev_shock+",all_controls,
  "|name+Country")),data=df_dyn_fe,cluster=~name+Country)})
res_FE_dd<-map(0:12,function(h){feols(as.formula(paste0(
  "cumF",h,"_financ_ext~d2d_shock+",all_controls,
  "|name+Country")),data=df_dyn_fe,cluster=~name+Country)})
FE_lev_coefs<-tibble(horizon=0:12,beta_lev=map_dbl(res_FE_lev,~coef(.x)["lev_shock"]),
                     se_lev=map_dbl(res_FE_lev,~sqrt(vcov(.x)["lev_shock","lev_shock"])))
FE_dd_coefs<-tibble(horizon=0:12,beta_dd=map_dbl(res_FE_dd,~coef(.x)["d2d_shock"]),
                    se_dd=map_dbl(res_FE_dd,~sqrt(vcov(.x)["d2d_shock","d2d_shock"])))

# Guardar resultados CSV
write_csv(FE_lev_coefs,file.path(base_dir,"dynamics_FE_lev_0_12.csv"))
write_csv(FE_dd_coefs,file.path(base_dir,"dynamics_FE_dd_0_12.csv"))

# Table for 6b
table_FE<-FE_lev_coefs%>%left_join(FE_dd_coefs,by="horizon")%>%
  mutate(Leverage=sprintf("%.3f (%.3f)",beta_lev,se_lev),
         `Dist-to-Default`=sprintf("%.3f (%.3f)",beta_dd,se_dd))%>%
  select(Horizon=horizon,Leverage,`Dist-to-Default`)

table_FE%>%kable(caption="Figure 6b: Respuestas dinámicas del proxy de financiamiento externo (horizonte 0–12)",
                 align=c("c","c","c"))%>%kable_styling(full_width=FALSE,bootstrap_options=c("striped","hover"))

# Gráficos 6b: financiamiento externo
p_FE_lev<-ggplot(FE_lev_coefs,aes(x=horizon,y=beta_lev))+geom_line(size=0.8,color="steelblue")+
  geom_point(size=2,color="steelblue")+geom_ribbon(aes(ymin=beta_lev-1.96*se_lev,ymax=beta_lev+1.96*se_lev),
                                                   alpha=0.2,fill="steelblue")+scale_x_continuous(breaks=0:12)+
  labs(title="Dinámica FE: Leverage × Shock",x="Horizonte (trimestres)",y="Δ Financiamiento Externo Acumulado")+
  theme_minimal()
p_FE_dd<-ggplot(FE_dd_coefs,aes(x=horizon,y=beta_dd))+geom_line(size=0.8,color="steelblue")+
  geom_point(size=2,color="steelblue")+geom_ribbon(aes(ymin=beta_dd-1.96*se_dd,ymax=beta_dd+1.96*se_dd),
                                                   alpha=0.2,fill="steelblue")+scale_x_continuous(breaks=0:12)+
  labs(title="Dinámica FE: Distance-to-Default × Shock",x="Horizonte (trimestres)",
       y="Δ Financiamiento Externo Acumulado")+theme_minimal()

p_FE_lev/p_FE_dd




# --------------------------------------------------------
# Figure 9: Dynamics, Not Controlling for Differences in Cyclical Sensitivities (horizonte 0–12)
# --------------------------------------------------------

# 2) Generar controles firm-level, estandarizar condiciones financieras y definir choque puro
if (!"rsales_g_std" %in% names(df)) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(rsales_g_std = { x <- log(saleq) - log(lag(saleq));
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }) %>% ungroup()
} else message("Variable 'rsales_g_std' ya existe; omitiendo.")

if (!"size_index" %in% names(df)) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(size_index = { x <- log(atq);
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }) %>% ungroup()
} else message("Variable 'size_index' ya existe; omitiendo.")

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
} else {
  message("Variable 'sh_current_a_std' ya existe; omitiendo.")
}


if (!all(c("lev_std","dd_std") %in% names(df))) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      lev_std = { x <- leverage; (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) },
      dd_std  = { x <- dd;       (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }
    ) %>% ungroup()
} else message("Variables 'lev_std' y 'dd_std' ya existen; omitiendo.")

if (!all(c("lev_shock","d2d_shock") %in% names(df))) {
  df <- df %>% mutate(
    lev_shock = lev_std * shock,
    d2d_shock = dd_std  * shock
  )
} else message("Variables 'lev_shock' y 'd2d_shock' ya existen; omitiendo.")

# 3) Construir cumFh_dlog_capital para h = 0…12
vars_cap <- paste0("cumF", 0:12, "_dlog_capital")
if (!all(vars_cap %in% names(df))) {
  df_dyn_nocy <- df %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~{
      tmp <- .x
      for (h in 0:12) {
        tmp[[paste0("cumF",h,"_dlog_capital")]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm=TRUE)
      }
      tmp
    }) %>% ungroup()
} else {
  message("Variables dinámicas cumFh_dlog_capital ya existen; omitiendo.")
  df_dyn_nocy <- df
}


# --------------------------------------------------------
# 4) Definir controles firm-level y macro (contemporáneos)
# --------------------------------------------------------
controls_firm  <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")

# Detectar cuáles de los controles macro realmente existen en df_dyn_nocy
present_macro <- intersect(controls_macro, names(df_dyn_nocy))
if (length(present_macro) < length(controls_macro)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro, present_macro), collapse = ", "),
    ". Se omitirán."
  )
}

# Construir la cadena de controles sólo con los firm-level y los macro presentes
all_controls_nocy <- paste(c(controls_firm, present_macro), collapse = " + ")

# --------------------------------------------------------
# 5) Estimar dinámicas para h = 0…12 sin componente cíclico
# --------------------------------------------------------
library(purrr)
library(fixest)

# (a) Leverage × shock
res_lev_nocy <- map(0:12, function(h) {
  fml <- as.formula(paste0(
    "cumF", h, "_dlog_capital ~ lev_shock + ", all_controls_nocy,
    " | name + Country"
  ))
  feols(fml,
        data    = df_dyn_nocy,
        cluster = ~ name + Country)
})

# (b) Distance-to-default × shock
res_dd_nocy <- map(0:12, function(h) {
  fml <- as.formula(paste0(
    "cumF", h, "_dlog_capital ~ d2d_shock + ", all_controls_nocy,
    " | name + Country"
  ))
  feols(fml,
        data    = df_dyn_nocy,
        cluster = ~ name + Country)
})



# 6) Extraer coeficientes y errores estándar
if (!exists("dyn_lev_nocy")) {
  dyn_lev_nocy <- tibble(
    horizon = 0:12,
    beta     = map_dbl(res_lev_nocy, ~ coef(.x)["lev_shock"]),
    se       = map_dbl(res_lev_nocy, ~ sqrt(vcov(.x)["lev_shock","lev_shock"]))
  )
} else message("Objeto 'dyn_lev_nocy' ya existe; omitiendo extracción.")

if (!exists("dyn_dd_nocy")) {
  dyn_dd_nocy <- tibble(
    horizon = 0:12,
    beta     = map_dbl(res_dd_nocy, ~ coef(.x)["d2d_shock"]),
    se       = map_dbl(res_dd_nocy, ~ sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
  )
} else message("Objeto 'dyn_dd_nocy' ya existe; omitiendo extracción.")

# 7) Graficar (horizonte 0–12)
p_nocy_lev <- ggplot(dyn_lev_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
              alpha = 0.2, fill = "steelblue") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figure 9(a): Dynamics (No Cyclical Sens.) – Leverage × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) + theme_minimal()

p_nocy_dd <- ggplot(dyn_dd_nocy, aes(x = horizon, y = beta)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
              alpha = 0.2, fill = "steelblue") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figure 9(b): Dynamics (No Cyclical Sens.) – DD × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) + theme_minimal()

(p_nocy_lev / p_nocy_dd) +
  plot_annotation(
    title = "Figure 9: Dynamics (No Cyclical Sensitivities) (horizonte 0–12)"
  )



# ------------------------------------------------------------
# Figura 10: Respuesta Promedio de Inversión al Shock Monetario (horizonte 0–12)
# ------------------------------------------------------------

# (2) Preprocesar y estandarizar variables firm-level y financieras
# ——————————————————————————————————————————————————————————————————

# 2.1) Crecimiento de ventas
if (!"rsales_g" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(rsales_g = log(saleq) - log(lag(saleq))) %>%
    ungroup()
} else message("Variable 'rsales_g' ya existe; omitiendo creación.")

if (!"rsales_g_std" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) /
        sd(rsales_g, na.rm = TRUE)
    ) %>%
    ungroup()
} else message("Variable 'rsales_g_std' ya existe; omitiendo creación.")

# 2.2) Tamaño (log activos)
if (!"size_index" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      size_index = (log(atq) - mean(log(atq), na.rm = TRUE)) /
        sd(log(atq), na.rm = TRUE)
    ) %>%
    ungroup()
} else message("Variable 'size_index' ya existe; omitiendo creación.")

# 2.3) Liquidez corriente
if (!"sh_current_a_std" %in% names(df)) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      sh_current_a_std = (current_ratio - mean(current_ratio, na.rm = TRUE)) /
        sd(current_ratio, na.rm = TRUE)
    ) %>%
    ungroup()
} else message("Variable 'sh_current_a_std' ya existe; omitiendo creación.")

# 2.4) Apalancamiento y distance-to-default
if (!all(c("lev_std","dd_std") %in% names(df))) {
  df <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    mutate(
      lev_std = (leverage - mean(leverage, na.rm = TRUE)) /
        sd(leverage, na.rm = TRUE),
      dd_std  = (dd        - mean(dd,        na.rm = TRUE)) /
        sd(dd,        na.rm = TRUE)
    ) %>%
    ungroup()
} else message("Variables 'lev_std' y 'dd_std' ya existen; omitiendo creación.")

# 2.5) Interacciones con el shock puro
if (!all(c("lev_shock","d2d_shock") %in% names(df))) {
  df <- df %>%
    mutate(
      lev_shock = lev_std * shock,
      d2d_shock = dd_std  * shock
    )
} else message("Variables 'lev_shock' y 'd2d_shock' ya existen; omitiendo creación.")

# (3) Construir variables dinámicas cumFh_dlog_capital para h = 0…12
# ——————————————————————————————————————————————————————————————————
vars_cap <- paste0("cumF", 0:12, "_dlog_capital")
if (!all(vars_cap %in% names(df))) {
  df_dyn <- df %>%
    group_by(name) %>%
    arrange(dateq) %>%
    group_modify(~{
      tmp <- .x
      for (h in 0:12) {
        tmp[[ paste0("cumF",h,"_dlog_capital") ]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  message("Variables dinámicas cumFh_dlog_capital ya existen; omitiendo creación.")
  df_dyn <- df
}

# --------------------------------------------------------
# 4) Definir controles firm-level y macro (contemporáneos)
# --------------------------------------------------------
controls_firm  <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")

# Solo conservamos los macro que realmente existen en df_dyn
present_macro <- intersect(controls_macro, names(df_dyn))
if (length(present_macro) < length(controls_macro)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro, present_macro), collapse = ", "),
    ". Se omitirán en la estimación."
  )
}

# Construimos all_controls con firm-level + los macro presentes
all_controls <- paste(c(controls_firm, present_macro), collapse = " + ")

# --------------------------------------------------------
# 5) Estimar dinámica de respuesta promedio al shock (horizonte 0–12)
# --------------------------------------------------------
library(purrr)
library(fixest)

res_avg <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep_var,
    " ~ shock + lev_shock + d2d_shock + ",
    all_controls,
    " | name + Country"
  ))
  feols(fml,
        data    = df_dyn,
        cluster = ~ name + Country)
})


# (6) Extraer coeficientes y errores estándar para 'shock'
# ——————————————————————————————————————————————————————————————————
if (!exists("avg_coefs")) {
  avg_coefs <- tibble(
    horizon    = 0:12,
    beta_shock = map_dbl(res_avg, ~ coef(.x)["shock"]),
    se_shock   = map_dbl(res_avg, ~ sqrt(vcov(.x)["shock","shock"]))
  )
} else message("Objeto 'avg_coefs' ya existe; omitiendo extracción.")

# (7) Graficar la respuesta promedio al shock
# ——————————————————————————————————————————————————————————————————
p_avg <- ggplot(avg_coefs, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_ribbon(aes(ymin = beta_shock - 1.96 * se_shock,
                  ymax = beta_shock + 1.96 * se_shock),
              alpha = 0.2, fill = "steelblue") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 10: Respuesta Promedio de Inversión al Shock Monetario",
    x     = "Horizonte (trimestres)",
    y     = "Efecto promedio acumulado de inversión"
  ) +
  theme_minimal()

print(p_avg)



# ------------------------------------------------------------
# Figura 10 (ajustada): Respuesta dinámica de inversión al Shock Monetario Heterogénea
# (horizonte 0–12) – sin estimar el efecto promedio
# ------------------------------------------------------------

# 2) Preprocesar y estandarizar variables firm-level y financieras
if (!"rsales_g" %in% names(df)) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(rsales_g = log(saleq) - log(lag(saleq))) %>% ungroup()
} else message("Variable 'rsales_g' ya existe; omitiendo.")
if (!"rsales_g_std" %in% names(df)) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(rsales_g_std = (rsales_g - mean(rsales_g, na.rm=TRUE)) / sd(rsales_g, na.rm=TRUE)) %>% ungroup()
} else message("Variable 'rsales_g_std' ya existe; omitiendo.")
if (!"size_index" %in% names(df)) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(size_index = (log(atq) - mean(log(atq), na.rm=TRUE)) / sd(log(atq), na.rm=TRUE)) %>% ungroup()
} else message("Variable 'size_index' ya existe; omitiendo.")
if (!"sh_current_a_std" %in% names(df)) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(sh_current_a_std = (current_ratio - mean(current_ratio, na.rm=TRUE)) / sd(current_ratio, na.rm=TRUE)) %>% ungroup()
} else message("Variable 'sh_current_a_std' ya existe; omitiendo.")
if (!all(c("lev_std","dd_std") %in% names(df))) {
  df <- df %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      lev_std = (leverage - mean(leverage, na.rm=TRUE)) / sd(leverage, na.rm=TRUE),
      dd_std  = (dd        - mean(dd,        na.rm=TRUE)) / sd(dd,        na.rm=TRUE)
    ) %>% ungroup()
} else message("Variables 'lev_std' y 'dd_std' ya existen; omitiendo.")
if (!all(c("lev_shock","d2d_shock") %in% names(df))) {
  df <- df %>% mutate(
    lev_shock = lev_std * shock,
    d2d_shock = dd_std  * shock
  )
} else message("Variables 'lev_shock' y 'd2d_shock' ya existen; omitiendo.")

# 3) Construir dinámicas cumFh_dlog_capital para h = 0…12
vars_cap10 <- paste0("cumF",0:12,"_dlog_capital")
if (!all(vars_cap10 %in% names(df))) {
  df_dyn10 <- df %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~ { tmp <- .x; for(h in 0:12) tmp[[paste0("cumF",h,"_dlog_capital")]] <-
        rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm=TRUE); tmp }) %>% ungroup()
} else { message("Variables dinámicas cumFh_dlog_capital ya existen; omitiendo."); df_dyn10 <- df }


# ------------------------------------------------------------
# 4) Definir controles firm-level y macro (contemporáneos)
# ------------------------------------------------------------
controls_firm10  <- c("rsales_g_std", "size_index", "sh_current_a_std")
controls_macro10 <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")

# Solo conservamos los controles macro que realmente existen en df_dyn10
present_macro10 <- intersect(controls_macro10, names(df_dyn10))
if (length(present_macro10) < length(controls_macro10)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro10, present_macro10), collapse = ", "),
    ". Se omitirán."
  )
}

# Construimos la cadena de controles con firm-level + los macro presentes
all_controls10 <- paste(c(controls_firm10, present_macro10), collapse = " + ")

# ------------------------------------------------------------
# 5) Estimar dinámicas sin componente promedio (horizonte 0–12)
# ------------------------------------------------------------
library(purrr)
library(fixest)

# (a) Leverage × shock
res_lev10 <- map(0:12, function(h) {
  dep <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep, " ~ lev_shock + ", all_controls10,
    " | name + Country"
  ))
  feols(fml,
        data    = df_dyn10,
        cluster = ~ name + Country)
})

# (b) Distance‐to‐default × shock
res_dd10 <- map(0:12, function(h) {
  dep <- paste0("cumF", h, "_dlog_capital")
  fml <- as.formula(paste0(
    dep, " ~ d2d_shock + ", all_controls10,
    " | name + Country"
  ))
  feols(fml,
        data    = df_dyn10,
        cluster = ~ name + Country)
})



# 6) Extraer coeficientes y errores
if (!exists("coef_lev10")) {
  coef_lev10 <- tibble(
    horizon = 0:12,
    beta_lev = map_dbl(res_lev10, ~ coef(.x)["lev_shock"]),
    se_lev   = map_dbl(res_lev10, ~ sqrt(vcov(.x)["lev_shock","lev_shock"]))
  )
} else message("Objeto 'coef_lev10' ya existe; omitiendo.")
if (!exists("coef_dd10")) {
  coef_dd10 <- tibble(
    horizon = 0:12,
    beta_dd  = map_dbl(res_dd10, ~ coef(.x)["d2d_shock"]),
    se_dd    = map_dbl(res_dd10, ~ sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
  )
} else message("Objeto 'coef_dd10' ya existe; omitiendo.")

# 7) Graficar dinámicas (horizonte 0–12)
p10_lev <- ggplot(coef_lev10, aes(x=horizon, y=beta_lev)) +
  geom_line(size=1, color="steelblue") + geom_point(size=2, color="steelblue") +
  geom_ribbon(aes(ymin=beta_lev-1.96*se_lev, ymax=beta_lev+1.96*se_lev), alpha=0.2, fill="steelblue") +
  scale_x_continuous(breaks=0:12) +
  labs(title="Dinámica: Leverage × Shock (sin promedio)", x="Horizonte", y="β") + theme_minimal()
p10_dd <- ggplot(coef_dd10, aes(x=horizon, y=beta_dd)) +
  geom_line(size=1, color="firebrick") + geom_point(size=2, color="firebrick") +
  geom_ribbon(aes(ymin=beta_dd-1.96*se_dd, ymax=beta_dd+1.96*se_dd), alpha=0.2, fill="firebrick") +
  scale_x_continuous(breaks=0:12) +
  labs(title="Dinámica: DD × Shock (sin promedio)", x="Horizonte", y="β") + theme_minimal()

(p10_lev / p10_dd) +
  plot_annotation(title = "Figura 10: Dinámica de Inversión sin Componente Promedio (0–12)")




# ------------------------------------------------------------
# Figura 12 (ajustada): Dinámica de Respuestas al Shock Monetario por Tamaño (toda la muestra)
# ------------------------------------------------------------

library(dplyr)
library(purrr)
library(fixest)
library(zoo)
library(ggplot2)

# 1) Partimos de df original
df_cntl <- df

# 2) Calcular ventas, tamaño, liquidez y apalancamiento estandarizado
df_cntl <- df_cntl %>%
  group_by(name) %>%
  arrange(dateq, .by_group = TRUE) %>%
  # 2.1) crecimiento de ventas
  mutate(
    rsales_g = if (!"rsales_g"   %in% names(.))   log(saleq) - log(lag(saleq)) else rsales_g,
    rsales_g_std = if (!"rsales_g_std" %in% names(.))
      (rsales_g - mean(rsales_g, na.rm = TRUE)) / sd(rsales_g, na.rm = TRUE)
    else rsales_g_std,
    # 2.2) size_index
    size_index = if (!"size_index" %in% names(.))
      (log(atq) - mean(log(atq), na.rm = TRUE)) / sd(log(atq), na.rm = TRUE)
    else size_index,
    # 2.3) liquidez corriente
    sh_current_a_std = if (!"sh_current_a_std" %in% names(.))
      (current_ratio - mean(current_ratio, na.rm = TRUE)) / sd(current_ratio, na.rm = TRUE)
    else sh_current_a_std,
    # 2.4) apalancamiento y distance-to-default
    lev_std = if (!"lev_std" %in% names(.))
      (leverage - mean(leverage, na.rm = TRUE)) / sd(leverage, na.rm = TRUE)
    else lev_std,
    dd_std  = if (!"dd_std"  %in% names(.))
      (dd        - mean(dd,        na.rm = TRUE)) / sd(dd,        na.rm = TRUE)
    else dd_std
  ) %>%
  ungroup()

# 3) Crear interacciones size × GDP y size × shock
# — solo si existen dlog_gdp y shock en df_cntl
if ("dlog_gdp" %in% names(df_cntl)) {
  df_cntl <- df_cntl %>%
    mutate(
      size_gdp   = if (!"size_gdp"   %in% names(.)) size_index * dlog_gdp else size_gdp
    )
} else {
  warning("No existe 'dlog_gdp'; se omite creación de 'size_gdp'.")
}

if ("shock" %in% names(df_cntl)) {
  df_cntl <- df_cntl %>%
    mutate(
      size_shock = if (!"size_shock" %in% names(.)) size_index * shock else size_shock
    )
} else {
  stop("Falta la variable 'shock' en el data frame.")
}

# 4) Generar dinámicas acumuladas cumFh_dlog_capital para h = 0…12
vars_cum <- paste0("cumF", 0:12, "_dlog_capital")
if (!all(vars_cum %in% names(df_cntl))) {
  df_dyn12 <- df_cntl %>%
    group_by(name) %>%
    arrange(dateq, .by_group = TRUE) %>%
    group_modify(~{
      tmp <- .x
      for (h in 0:12) {
        tmp[[ paste0("cumF", h, "_dlog_capital") ]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  message("Variables dinámicas cumFh_dlog_capital ya existen; omitiendo.")
  df_dyn12 <- df_cntl
}

# 5) Definir controles firm-level y macro (solo los macro presentes)
controls_firm  <- c("rsales_g_std", "sh_current_a_std")
controls_macro <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")
present_macro  <- intersect(controls_macro, names(df_dyn12))
if (length(present_macro) < length(controls_macro)) {
  warning("Faltan controles macro: ",
          paste(setdiff(controls_macro, present_macro), collapse = ", "),
          ". Se omitirán.")
}
all_controls <- c(controls_firm, present_macro) %>% paste(collapse = " + ")


# ------------------------------------------------------------
# 6) Estimar efecto heterogéneo size × shock para h = 0…12 (ajustado)
# ------------------------------------------------------------
library(purrr)
library(fixest)

# Flags de presencia
has_size_gdp   <- "size_gdp"   %in% names(df_dyn12)
present_macro  <- intersect(controls_macro, names(df_dyn12))

res_size <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  
  # Construir dinámicamente los términos RHS
  rhs_terms <- c(
    if (has_size_gdp) "size_gdp",
    "size_shock",
    controls_firm,
    present_macro
  )
  
  fml_str <- paste0(
    dep_var, " ~ ",
    paste(rhs_terms, collapse = " + "),
    " | name + Country"
  )
  
  feols(
    as.formula(fml_str),
    data    = df_dyn12,
    cluster = ~ name + Country
  )
})


# 7) Extraer coeficiente y error estándar de size_shock
if (!exists("coef_size")) {
  coef_size <- tibble(
    horizon    = 0:12,
    beta_shock = map_dbl(res_size, ~ coef(.x)["size_shock"]),
    se_shock   = map_dbl(res_size, ~ sqrt(vcov(.x)["size_shock", "size_shock"]))
  )
} else {
  message("Objeto 'coef_size' ya existe; omitiendo extracción.")
}

# 8) Graficar la dinámica de size_shock
ggplot(coef_size, aes(x = horizon, y = beta_shock)) +
  geom_line(size = 1, color = "darkgreen") +
  geom_point(size = 2, color = "darkgreen") +
  geom_ribbon(aes(
    ymin = beta_shock - 1.96 * se_shock,
    ymax = beta_shock + 1.96 * se_shock
  ), alpha = 0.2, fill = "darkgreen") +
  scale_x_continuous(breaks = 0:12) +
  labs(
    title = "Figura 12: Dinámica – Size × Shock (toda la muestra)",
    x     = "Horizonte (trimestres)",
    y     = "Coeficiente β_h"
  ) +
  theme_minimal()




# ------------------------------------------------------------
# Figura 13 (ajustada): Dinámica Conjunta de Posición Financiera y Tamaño
# (horizonte 0–12 para toda la muestra)
# ------------------------------------------------------------

# 1) Inicializar df_cntl13 a partir de df (debe tener dateq, name, Country)
# ------------------------------------------------------------------
df_cntl13 <- df

# 2) Crear/Verificar series macro: dlog_gdp y dlog_cpi (ajustado)
# ------------------------------------------------------------------
# dlog_gdp
if ("gdp" %in% names(df_cntl13)) {
  if (!"dlog_gdp" %in% names(df_cntl13)) {
    df_cntl13 <- df_cntl13 %>%
      group_by(Country) %>%
      arrange(dateq, .by_group = TRUE) %>%
      mutate(
        dlog_gdp = if_else(
          !is.na(gdp) & gdp > 0 &
            !is.na(lag(gdp)) & lag(gdp) > 0,
          log(gdp) - log(lag(gdp)),
          NA_real_
        )
      ) %>%
      ungroup()
  } else {
    message("Variable 'dlog_gdp' ya existe; omitiendo.")
  }
} else {
  warning("No se encontró 'gdp' en df_cntl13; se omite creación de 'dlog_gdp'.")
}

# dlog_cpi
if ("ipc" %in% names(df_cntl13)) {
  if (!"dlog_cpi" %in% names(df_cntl13)) {
    df_cntl13 <- df_cntl13 %>%
      group_by(Country) %>%
      arrange(dateq, .by_group = TRUE) %>%
      mutate(
        dlog_cpi = if_else(
          !is.na(ipc) & ipc > 0 &
            !is.na(lag(ipc)) & lag(ipc) > 0,
          log(ipc) - log(lag(ipc)),
          NA_real_
        )
      ) %>%
      ungroup()
  } else {
    message("Variable 'dlog_cpi' ya existe; omitiendo.")
  }
} else {
  warning("No se encontró 'ipc' en df_cntl13; se omite creación de 'dlog_cpi'.")
}


# 3) Crear proxy size_index_std y condiciones financieras
# ------------------------------------------------------------------
if (!"size_index_std" %in% names(df_cntl13)) {
  df_cntl13 <- df_cntl13 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      log_atq        = log(atq),
      log_atq_std    = (log_atq - mean(log_atq, na.rm=TRUE))/sd(log_atq, na.rm=TRUE),
      saleq_std      = (saleq   - mean(saleq,   na.rm=TRUE))/sd(saleq,   na.rm=TRUE),
      size_index     = (log_atq_std + saleq_std)/2,
      size_index_std = (size_index - mean(size_index, na.rm=TRUE))/sd(size_index, na.rm=TRUE),
      lev_std        = (leverage - mean(leverage, na.rm=TRUE))/sd(leverage, na.rm=TRUE),
      dd_std         = (dd        - mean(dd,        na.rm=TRUE))/sd(dd,        na.rm=TRUE)
    ) %>%
    ungroup() %>%
    select(-log_atq,-log_atq_std,-saleq_std)
} else message("Variables de tamaño y condiciones financieras ya existen; omitiendo.")

# 4) Crear interacciones
# ------------------------------------------------------------------
if (!all(c("size_gdp","size_shock","lev_shock_gdp","lev_shock","d2d_shock_gdp","d2d_shock") %in% names(df_cntl13))) {
  df_cntl13 <- df_cntl13 %>% mutate(
    size_gdp      = size_index_std * dlog_gdp,
    size_shock    = size_index_std * shock,
    lev_shock_gdp = lev_std          * dlog_gdp,
    lev_shock     = lev_std          * shock,
    d2d_shock_gdp = dd_std           * dlog_gdp,
    d2d_shock     = dd_std           * shock
  )
} else message("Interacciones ya existen; omitiendo.")

# 5) Crear dlog_capital y controles adicionales
# ------------------------------------------------------------------
if (!"dlog_capital" %in% names(df_cntl13)) {
  df_cntl13 <- df_cntl13 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      dlog_capital     = log(capital) - log(lag(capital)),
      rsales_g_std     = {(x<-log(saleq)-log(lag(saleq)));(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)},
      sh_current_a_std = {(x<-current_ratio);(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
    ) %>% ungroup()
} else message("Control de inversión o liquidez ya existe; omitiendo.")

# 6) Construir dinámicas cumFh_dlog_capital para h=0…12
# ------------------------------------------------------------------
vars13 <- paste0("cumF",0:12,"_dlog_capital")
if (!all(vars13 %in% names(df_cntl13))) {
  df_dyn13 <- df_cntl13 %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~{tmp<-.x;for(h in 0:12) tmp[[paste0("cumF",h,"_dlog_capital")]]<-rowSums(map_dfc(0:h,~lead(tmp$dlog_capital,.x)),na.rm=TRUE);tmp})%>%ungroup()
} else { message("Variables dinámicas ya existen; omitiendo."); df_dyn13<-df_cntl13 }



# ------------------------------------------------------------
# 7.5) Definir controles firm-level y macro (ajustado)
# ------------------------------------------------------------
controls_firm13  <- c("rsales_g_std", "sh_current_a_std", "size_index_std")
controls_macro13 <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")

present_macro13 <- intersect(controls_macro13, names(df_dyn13))
if (length(present_macro13) < length(controls_macro13)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro13, present_macro13), collapse = ", "),
    ". Se omitirán."
  )
}

all_controls13 <- paste(c(controls_firm13, present_macro13), collapse = " + ")

# ------------------------------------------------------------
# 8a) Estimar dinámica Leverage + Size (ajustado)
# ------------------------------------------------------------
res_lev_size <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var, " ~ size_shock + lev_shock + lev_shock_gdp + ",
    all_controls13,
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn13,
    cluster = ~ name + Country
  )
})



# 8b) Estimar dinámica DD + Size
# ------------------------------------------------------------------
res_dd_size <- map(0:12, ~feols(
  as.formula(paste0(
    "cumF",.x,"_dlog_capital ~ size_shock + d2d_shock + d2d_shock_gdp + ",
    all_controls13," | name + Country"
  )),data=df_dyn13,cluster=~name+dateq))

# 9) Extraer coeficientes y errores
# ------------------------------------------------------------------
lev_size_coefs <- tibble(
  horizon   = 0:12,
  beta_size = map_dbl(res_lev_size,~coef(.x)["size_shock"]),
  se_size   = map_dbl(res_lev_size,~sqrt(vcov(.x)["size_shock","size_shock"])),
  beta_lev  = map_dbl(res_lev_size,~coef(.x)["lev_shock"]),
  se_lev    = map_dbl(res_lev_size,~sqrt(vcov(.x)["lev_shock","lev_shock"]))
)
dd_size_coefs <- tibble(
  horizon   = 0:12,
  beta_size = map_dbl(res_dd_size,~coef(.x)["size_shock"]),
  se_size   = map_dbl(res_dd_size,~sqrt(vcov(.x)["size_shock","size_shock"])),
  beta_dd   = map_dbl(res_dd_size,~coef(.x)["d2d_shock"]),
  se_dd     = map_dbl(res_dd_size,~sqrt(vcov(.x)["d2d_shock","d2d_shock"]))
)

# 10) Guardar resultados opcional
# ------------------------------------------------------------------
write_csv(lev_size_coefs, file.path(base_dir,"dynamics_lev_size.csv"))
write_csv(dd_size_coefs,  file.path(base_dir,"dynamics_dd_size.csv"))

# 11a) Gráfico Size × Shock
p13a <- ggplot(lev_size_coefs, aes(x = horizon)) +
  geom_line(aes(y = beta_size), size = 1, color = "steelblue") +
  geom_ribbon(aes(ymin = beta_size - 1.96 * se_size,
                  ymax = beta_size + 1.96 * se_size),
              alpha = 0.2, fill = "steelblue") +
  labs(
    title = "Figura 13(a): Size × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Coeficiente β_h"
  ) +
  scale_x_continuous(breaks = 0:12) +
  theme_minimal()

# 11b) Gráfico Leverage × Shock
p13b <- ggplot(lev_size_coefs, aes(x = horizon)) +
  geom_line(aes(y = beta_lev), size = 1, color = "firebrick") +
  geom_ribbon(aes(ymin = beta_lev - 1.96 * se_lev,
                  ymax = beta_lev + 1.96 * se_lev),
              alpha = 0.2, fill = "firebrick") +
  labs(
    title = "Figura 13(b): Leverage × Shock",
    x     = "Horizonte (trimestres)",
    y     = "Coeficiente β_h"
  ) +
  scale_x_continuous(breaks = 0:12) +
  theme_minimal()

# 12) Mostrar ambas gráficas uno debajo de la otra
(p13a / p13b) +
  plot_annotation(
    title = "Figura 13: Dinámica Conjunta de Posición Financiera y Tamaño"
  )




# ------------------------------------------------------------
# Figura 14 (ajustada): Dinámica Conjunta de Posición Financiera (Leverage o DD) y Edad
# (toda la muestra)
# ------------------------------------------------------------

# Partir de df con dateq ya trimestral y name, Country
# ------------------------------------------------------------------
df_cntl14 <- df

# 1) Crear dlog_capital y Ldl_capital si faltan
if (!all(c("dlog_capital","Ldl_capital") %in% names(df_cntl14))) {
  df_cntl14 <- df_cntl14 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      capital2     = if_else(!is.na(capital)&capital>0, capital, NA_real_),
      log_capital  = log(capital2),
      lag_log_cap  = lag(log_capital),
      dlog_capital = log_capital - lag_log_cap,
      Ldl_capital  = lag(dlog_capital,1)
    ) %>% ungroup() %>% select(-capital2, -log_capital, -lag_log_cap)
} else message("Variables dlog_capital y Ldl_capital ya existen; omitiendo.")

# 2) Crear controles firm-level si faltan
if (!"rsales_g_std" %in% names(df_cntl14)) {
  df_cntl14 <- df_cntl14 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      saleq2       = if_else(!is.na(saleq)&saleq>0, saleq, NA_real_),
      rsales_g     = log(saleq2) - lag(log(saleq2)),
      rsales_g_std = (rsales_g - mean(rsales_g, na.rm=TRUE)) / sd(rsales_g, na.rm=TRUE)
    ) %>% ungroup() %>% select(-saleq2, -rsales_g)
} else message("Variable rsales_g_std ya existe; omitiendo.")

if (!all(c("lev_std","dd_std") %in% names(df_cntl14))) {
  df_cntl14 <- df_cntl14 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      lev_std = (leverage - mean(leverage, na.rm=TRUE)) / sd(leverage, na.rm=TRUE),
      dd_std  = (dd        - mean(dd,        na.rm=TRUE)) / sd(dd,        na.rm=TRUE)
    ) %>% ungroup()
} else message("Variables lev_std y dd_std ya existen; omitiendo.")

# 3) Crear dummies de edad (joven, middle, old) si faltan
if (!all(c("AgeYoung","AgeMiddle","AgeOld") %in% names(df_cntl14))) {
  df_cntl14 <- df_cntl14 %>%
    group_by(dateq) %>%
    mutate(
      p25 = quantile(age_sample, 0.25, na.rm=TRUE),
      p75 = quantile(age_sample, 0.75, na.rm=TRUE)
    ) %>%
    ungroup() %>%
    group_by(name) %>% arrange(dateq) %>%
    mutate(
      AgeYoung  = as.numeric(age_sample <  p25),
      AgeMiddle = as.numeric(age_sample >= p25 & age_sample < p75),
      AgeOld    = as.numeric(age_sample >= p75)
    ) %>%
    ungroup() %>%
    select(-p25, -p75)
} else message("Dummies de edad ya existen; omitiendo.")


# 4) Crear interacciones con GDP y Shock (sólo si existen las variables)
# -----------------------------------------------------------------------
# (a) Interacciones con dlog_gdp
if ("dlog_gdp" %in% names(df_cntl14)) {
  df_cntl14 <- df_cntl14 %>%
    mutate(
      AgeMiddle_gdp = AgeMiddle * dlog_gdp,
      AgeOld_gdp    = AgeOld    * dlog_gdp,
      lev_gdp       = lev_std   * dlog_gdp,
      dd_gdp        = dd_std    * dlog_gdp
    )
} else {
  warning("La variable 'dlog_gdp' no existe en df_cntl14; se omiten interacciones con GDP.")
}

# (b) Interacciones con shock
if ("shock" %in% names(df_cntl14)) {
  df_cntl14 <- df_cntl14 %>%
    mutate(
      AgeMiddle_shock = AgeMiddle * shock,
      AgeOld_shock    = AgeOld    * shock,
      lev_shock       = lev_std   * shock,
      dd_shock        = dd_std    * shock
    )
} else {
  stop("Falta la variable 'shock' en df_cntl14; no se pueden crear interacciones con shock.")
}

# 5) Construir dinámicas acumuladas cumFh_dlog_capital h=0…12
vars14 <- paste0("cumF",0:12,"_dlog_capital")
if (!all(vars14 %in% names(df_cntl14))) {
  df_dyn14 <- df_cntl14 %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~{tmp<-.x;for(h in 0:12) tmp[[paste0("cumF",h,"_dlog_capital")]]<-rowSums(map_dfc(0:h,~lead(tmp$dlog_capital,.x)),na.rm=TRUE);tmp})%>%ungroup()
} else { message("Variables dinámicas ya existen; omitiendo."); df_dyn14<-df_cntl14 }

# 6) Definir controles firm-level y macro
controls_firm14  <- c("rsales_g_std","Ldl_capital")
controls_macro14 <- c("dlog_gdp","dlog_cpi","unemp","embigl")
all_controls14   <- paste(c(controls_firm14,controls_macro14), collapse=" + ")



library(purrr)
library(fixest)

# 7a) Estimar Joint Dynamics – Leverage y Edad (ajustado)
# --------------------------------------------------------

# Detectar qué interacciones y controles están presentes
age_gdp_terms    <- intersect(c("AgeMiddle_gdp", "AgeOld_gdp"), names(df_dyn14))
age_shock_terms  <- intersect(c("AgeMiddle_shock", "AgeOld_shock"), names(df_dyn14))
lev_gdp_term     <- if ("lev_gdp"   %in% names(df_dyn14)) "lev_gdp"   else NULL
lev_shock_term   <- if ("lev_shock" %in% names(df_dyn14)) "lev_shock" else NULL

firm_terms       <- intersect(controls_firm14, names(df_dyn14))
macro_terms      <- intersect(controls_macro14, names(df_dyn14))

# Combinar todos los RHS terms en orden lógico
rhs_base <- c(age_gdp_terms,
              age_shock_terms,
              lev_gdp_term,
              lev_shock_term,
              firm_terms,
              macro_terms)

res_lev_age <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var, " ~ ",
    paste(rhs_base, collapse = " + "),
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn14,
    cluster = ~ name + Country
  )
})

# ------------------------------------------------------------
# 7b) Estimar Joint Dynamics – DD y Edad (ajustado)
# ------------------------------------------------------------
# Detectar qué términos de DD & Edad están presentes
age_gdp_terms    <- intersect(c("AgeMiddle_gdp", "AgeOld_gdp"), names(df_dyn14))
age_shock_terms  <- intersect(c("AgeMiddle_shock", "AgeOld_shock"), names(df_dyn14))
dd_gdp_term      <- if ("dd_gdp"    %in% names(df_dyn14)) "dd_gdp"    else NULL
dd_shock_term    <- if ("dd_shock"  %in% names(df_dyn14)) "dd_shock"  else NULL

# Ya tienes controls_firm14 y controls_macro14 de antes
firm_terms       <- intersect(controls_firm14, names(df_dyn14))
macro_terms      <- intersect(controls_macro14, names(df_dyn14))

rhs_base_dd <- c(
  age_gdp_terms,
  age_shock_terms,
  dd_gdp_term,
  dd_shock_term,
  firm_terms,
  macro_terms
)

res_dd_age <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var, " ~ ",
    paste(rhs_base_dd, collapse = " + "),
    " | name + Country"
  )
  feols(as.formula(fml_str),
        data    = df_dyn14,
        cluster = ~ name + Country)
})


# ------------------------------------------------------------
# 8) Extraer coeficientes y errores para DD & Edad
# ------------------------------------------------------------
dd_age_coefs <- tibble(
  horizon      = 0:12,
  beta_mid_age = map_dbl(res_dd_age, ~ coef(.x)["AgeMiddle_shock"]),
  se_mid_age   = map_dbl(res_dd_age, ~ sqrt(vcov(.x)["AgeMiddle_shock", "AgeMiddle_shock"])),
  beta_old_age = map_dbl(res_dd_age, ~ coef(.x)["AgeOld_shock"]),
  se_old_age   = map_dbl(res_dd_age, ~ sqrt(vcov(.x)["AgeOld_shock", "AgeOld_shock"])),
  beta_dd_shock= map_dbl(res_dd_age, ~ coef(.x)["dd_shock"]),
  se_dd_shock  = map_dbl(res_dd_age, ~ sqrt(vcov(.x)["dd_shock", "dd_shock"]))
)


# 9) Guardar resultados opcional
write_csv(lev_age_coefs, file.path(base_dir,"dynamics_lev_age.csv"))
write_csv(dd_age_coefs,  file.path(base_dir,"dynamics_dd_age.csv"))

# 10a) Gráfico Joint Dynamics – Leverage & Edad
# Preparar datos en formato largo para leyenda
lev_plot <- lev_age_coefs %>%
  select(horizon,
         AgeMiddle = beta_mid_age, AgeOld = beta_old_age, Leverage = beta_lev_shock) %>%
  pivot_longer(-horizon, names_to = "series", values_to = "beta")
se_plot_lev <- lev_age_coefs %>%
  select(horizon,
         AgeMiddle = se_mid_age, AgeOld = se_old_age, Leverage = se_lev_shock) %>%
  pivot_longer(-horizon, names_to = "series", values_to = "se")
lev_plot <- left_join(lev_plot, se_plot_lev, by = c("horizon", "series"))

p14a <- ggplot(lev_plot, aes(x = horizon, y = beta, color = series, fill = series)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), alpha = 0.2) +
  scale_color_manual(
    values = c("AgeMiddle" = "blue", "AgeOld" = "orange", "Leverage" = "red"),
    labels = c("Edad media", "Edad alta", "Leverage × Shock"),
    name = "Serie"
  ) +
  scale_fill_manual(
    values = c("AgeMiddle" = "blue", "AgeOld" = "orange", "Leverage" = "red"),
    labels = c("Edad media", "Edad alta", "Leverage × Shock"),
    name = "Serie"
  ) +
  labs(
    title = "Figura 14(a): Leverage & Edad",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  scale_x_continuous(breaks = 0:12) +
  theme_minimal()

# 10b) Gráfico Joint Dynamics – DD & Edad
# Preparar datos en formato largo
dd_plot <- dd_age_coefs %>%
  select(horizon,
         AgeMiddle = beta_mid_age, AgeOld = beta_old_age, DD = beta_dd_shock) %>%
  pivot_longer(-horizon, names_to = "series", values_to = "beta")
se_plot_dd <- dd_age_coefs %>%
  select(horizon,
         AgeMiddle = se_mid_age, AgeOld = se_old_age, DD = se_dd_shock) %>%
  pivot_longer(-horizon, names_to = "series", values_to = "se")
dd_plot <- left_join(dd_plot, se_plot_dd, by = c("horizon", "series"))

p14b <- ggplot(dd_plot, aes(x = horizon, y = beta, color = series, fill = series)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), alpha = 0.2) +
  scale_color_manual(
    values = c("AgeMiddle" = "blue", "AgeOld" = "orange", "DD" = "darkgreen"),
    labels = c("Edad media", "Edad alta", "DD × Shock"),
    name = "Serie"
  ) +
  scale_fill_manual(
    values = c("AgeMiddle" = "blue", "AgeOld" = "orange", "DD" = "darkgreen"),
    labels = c("Edad media", "Edad alta", "DD × Shock"),
    name = "Serie"
  ) +
  labs(
    title = "Figura 14(b): DD & Edad",
    x     = "Horizonte (trimestres)",
    y     = "Coeficiente β_h"
  ) +
  scale_x_continuous(breaks = 0:12) +
  theme_minimal()

# 11) Mostrar ambas gráficas apiladas
(p14a / p14b) +
  plot_annotation(title = "Figura 14: Joint Dynamics Posición Financiera y Edad")
p14b <- ggplot(dd_age_coefs, aes(x=horizon)) +
  geom_line(aes(y=beta_mid_age), size=1, color="blue") +
  geom_ribbon(aes(ymin=beta_mid_age-1.96*se_mid_age,ymax=beta_mid_age+1.96*se_mid_age), fill="blue", alpha=0.2) +
  geom_line(aes(y=beta_old_age), size=1, color="orange") +
  geom_ribbon(aes(ymin=beta_old_age-1.96*se_old_age,ymax=beta_old_age+1.96*se_old_age), fill="orange", alpha=0.2) +
  geom_line(aes(y=beta_dd_shock), size=1, color="darkgreen") +
  geom_ribbon(aes(ymin=beta_dd_shock-1.96*se_dd_shock,ymax=beta_dd_shock+1.96*se_dd_shock), fill="darkgreen", alpha=0.2) +
  labs(
    title = "Figura 14(b): DD & Edad",
    x     = "Horizonte (trimestres)",
    y     = "Efecto acumulado de inversión"
  ) +
  scale_x_continuous(breaks = 0:12) +
  theme_minimal()

# 11) Mostrar ambas gráficas apiladas
(p14a / p14b) +
  plot_annotation(title = "Figura 14: Joint Dynamics Posición Financiera y Edad")




# ------------------------------------------------------------
# Figura 16 (ajustada): Dinámica Conjunta de Posición Financiera (Leverage o DD) y Liquidez
# (toda la muestra)
# ------------------------------------------------------------

df_cntl16 <- df

# 1) dlog_capital y Ldl_capital
if (!all(c("dlog_capital","Ldl_capital") %in% names(df_cntl16))) {
  df_cntl16 <- df_cntl16 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      capital2     = if_else(!is.na(capital)&capital>0, capital, NA_real_),
      log_capital  = log(capital2),
      dlog_capital = log_capital - lag(log_capital),
      Ldl_capital  = lag(dlog_capital)
    ) %>% ungroup() %>% select(-capital2, -log_capital)
} else message("dlog_capital y Ldl_capital ya existen; omitiendo.")

# 2) Controles firm-level
if (!"rsales_g_std" %in% names(df_cntl16)) {
  df_cntl16 <- df_cntl16 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      rsales_g     = log(saleq) - log(lag(saleq)),
      rsales_g_std = (rsales_g - mean(rsales_g,na.rm=TRUE))/sd(rsales_g,na.rm=TRUE)
    ) %>% ungroup() %>% select(-rsales_g)
} else message("rsales_g_std ya existe; omitiendo.")

if (!all(c("lev_std","dd_std","liq_std") %in% names(df_cntl16))) {
  df_cntl16 <- df_cntl16 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      lev_std  = (leverage - mean(leverage,na.rm=TRUE))/sd(leverage,na.rm=TRUE),
      dd_std   = (dd - mean(dd,na.rm=TRUE))/sd(dd,na.rm=TRUE),
      liq_std  = (current_ratio - mean(current_ratio,na.rm=TRUE))/sd(current_ratio,na.rm=TRUE)
    ) %>% ungroup()
} else message("lev_std, dd_std o liq_std ya existen; omitiendo.")

# 3) Interacciones
ints16 <- c("lev_shock","dd_shock","liq_shock")
if (!all(ints16 %in% names(df_cntl16))) {
  df_cntl16 <- df_cntl16 %>% mutate(
    lev_shock = lev_std * shock,
    dd_shock  = dd_std  * shock,
    liq_shock = liq_std * shock
  )
} else message("Interacciones de choque ya existen; omitiendo.")

# 4) Dinámicas acumuladas
vars16 <- paste0("cumF",0:12,"_dlog_capital")
if (!all(vars16 %in% names(df_cntl16))) {
  df_dyn16 <- df_cntl16 %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~{tmp<-.x; for(h in 0:12) tmp[[paste0("cumF",h,"_dlog_capital")]]<-rowSums(map_dfc(0:h,~lead(tmp$dlog_capital,.x)),na.rm=TRUE); tmp}) %>% ungroup()
} else { message("Variables dinámicas ya existen; omitiendo."); df_dyn16 <- df_cntl16 }

# ------------------------------------------------------------
# 5) Definir controles firm-level y macro (solo los macro presentes)
# ------------------------------------------------------------
controls_firm  <- c("rsales_g_std", "Ldl_capital")
controls_macro <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")

# Detectar qué controles macro existen en df_dyn16
present_macro <- intersect(controls_macro, names(df_dyn16))
if (length(present_macro) < length(controls_macro)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro, present_macro), collapse = ", "),
    ". Se omitirán."
  )
}

# Construir vector de todos los RHS
rhs_base <- c("lev_shock", "liq_shock", controls_firm, present_macro)

# ------------------------------------------------------------
# 6a) Estimar Leverage & Liquidez
# ------------------------------------------------------------
res_lev_liq <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var, " ~ ",
    paste(rhs_base, collapse = " + "),
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn16,
    cluster = ~ name + Country
  )
})

# 6b) Estimar DD & Liquidez
res_dd_liq <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var, " ~ dd_shock + liq_shock + ",
    paste(rhs_base[-1], collapse = " + "),  # omitimos lev_shock
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn16,
    cluster = ~ name + Country
  )
})

# 7) Extraer coeficientes
tbl_lev_liq <- tibble(
  horizon      = 0:12,
  lev_beta     = map_dbl(res_lev_liq,~coef(.x)["lev_shock"]),
  lev_se       = map_dbl(res_lev_liq,~sqrt(vcov(.x)["lev_shock","lev_shock"])),
  liq_beta     = map_dbl(res_lev_liq,~coef(.x)["liq_shock"]),
  liq_se       = map_dbl(res_lev_liq,~sqrt(vcov(.x)["liq_shock","liq_shock"]))
)
tbl_dd_liq <- tibble(
  horizon      = 0:12,
  dd_beta      = map_dbl(res_dd_liq,~coef(.x)["dd_shock"]),
  dd_se        = map_dbl(res_dd_liq,~sqrt(vcov(.x)["dd_shock","dd_shock"])),
  liq_beta     = map_dbl(res_dd_liq,~coef(.x)["liq_shock"]),
  liq_se       = map_dbl(res_dd_liq,~sqrt(vcov(.x)["liq_shock","liq_shock"]))
)

# 8) Gráficos con leyenda
data16a <- tbl_lev_liq %>% select(horizon, Leverage=lev_beta, Liquidez=liq_beta) %>% pivot_longer(-horizon,names_to="serie",values_to="beta")
data16a_se <- tbl_lev_liq %>% select(horizon, Leverage=lev_se, Liquidez=liq_se) %>% pivot_longer(-horizon,names_to="serie",values_to="se")
data16a <- left_join(data16a,data16a_se,by=c("horizon","serie"))

p16a <- ggplot(data16a,aes(x=horizon,y=beta,color=serie,fill=serie))+
  geom_line(size=1)+geom_ribbon(aes(ymin=beta-1.96*se,ymax=beta+1.96*se),alpha=0.2)+
  scale_color_manual(values=c("Leverage"="red","Liquidez"="blue"),name="Serie")+
  scale_fill_manual(values=c("Leverage"="red","Liquidez"="blue"),name="Serie")+
  labs(title="Figura 16(a): Leverage & Liquidez",x="Horizonte",y="Efecto acumulado de inversión")+
  scale_x_continuous(breaks=0:12)+theme_minimal()

# 9) DD & Liquidez
data16b <- tbl_dd_liq %>% select(horizon, DD=dd_beta, Liquidez=liq_beta) %>% pivot_longer(-horizon,names_to="serie",values_to="beta")
data16b_se <- tbl_dd_liq %>% select(horizon, DD=dd_se, Liquidez=liq_se) %>% pivot_longer(-horizon,names_to="serie",values_to="se")
data16b <- left_join(data16b,data16b_se,by=c("horizon","serie"))

p16b <- ggplot(data16b,aes(x=horizon,y=beta,color=serie,fill=serie))+
  geom_line(size=1)+geom_ribbon(aes(ymin=beta-1.96*se,ymax=beta+1.96*se),alpha=0.2)+
  scale_color_manual(values=c("DD"="darkgreen","Liquidez"="blue"),name="Serie")+
  scale_fill_manual(values=c("DD"="darkgreen","Liquidez"="blue"),name="Serie")+
  labs(title="Figura 16(b): DD & Liquidez",x="Horizonte",y="Efecto acumulado de inversión")+
  scale_x_continuous(breaks=0:12)+theme_minimal()

# 10) Mostrar
graph <- p16a / p16b + plot_annotation(title="Figura 16: Joint Dynamics Liquidez y Estado Financiero")
graph



# ------------------------------------------------------------
# Figura 21 (ajustada): Dinámica de Respuesta Diferencial a Choques Monetarios
# por Volatilidad de Ventas 5‑Year (toda la muestra)
# ------------------------------------------------------------

df_cntl21 <- df

# 1) dlog_capital y Ldl_capital
if (!all(c("dlog_capital","Ldl_capital") %in% names(df_cntl21))) {
  df_cntl21 <- df_cntl21 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      capital2     = if_else(!is.na(capital)&capital>0,capital,NA_real_),
      log_capital  = log(capital2),
      dlog_capital = log_capital - lag(log_capital),
      Ldl_capital  = lag(dlog_capital)
    ) %>% ungroup() %>% select(-capital2,-log_capital)
} else message("dlog_capital y Ldl_capital ya existen; omitiendo.")

# # 2) Controles firm-level: ventas y liquidez
if (!"rsales_g_std" %in% names(df_cntl21)) {
  df_cntl21 <- df_cntl21 %>%
    group_by(name) %>%
    arrange(dateq, .by_group = TRUE) %>%
    mutate(
      rsales_g     = log(saleq) - log(lag(saleq)),
      rsales_g_std = (rsales_g - mean(rsales_g, na.rm = TRUE)) /
        sd(rsales_g, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(-rsales_g)
} else {
  message("rsales_g_std ya existe; omitiendo.")
}

if (!"liq_std" %in% names(df_cntl21)) {
  df_cntl21 <- df_cntl21 %>%
    group_by(name) %>%
    arrange(dateq, .by_group = TRUE) %>%
    mutate(
      liq_std = (current_ratio - mean(current_ratio, na.rm = TRUE)) /
        sd(current_ratio, na.rm = TRUE)
    ) %>%
    ungroup()
} else {
  message("liq_std ya existe; omitiendo.")
}

# 3) Volatilidad 5-Year
if (!"vol_5yr_std" %in% names(df_cntl21)) {
  df_cntl21 <- df_cntl21 %>%
    group_by(name) %>%
    arrange(dateq, .by_group = TRUE) %>%
    mutate(
      yoy_sales    = log(saleq) - lag(log(saleq), 4),
      vol_5yr      = zoo::rollapply(yoy_sales, 20, sd, na.rm = TRUE, fill = NA, align = "right"),
      vol_5yr_std  = (vol_5yr - mean(vol_5yr, na.rm = TRUE)) / sd(vol_5yr, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(-yoy_sales, -vol_5yr)
} else {
  message("vol_5yr_std ya existe; omitiendo.")
}


# 4) Interacciones con Volatilidad 5-Year
if (!all(c("vol_5yr_std_gdp", "vol_5yr_std_shock") %in% names(df_cntl21))) {
  
  # (a) Sólo si existe dlog_gdp
  if ("dlog_gdp" %in% names(df_cntl21)) {
    df_cntl21 <- df_cntl21 %>%
      mutate(
        vol_5yr_std_gdp   = vol_5yr_std * dlog_gdp
      )
  } else {
    warning("No existe 'dlog_gdp'; se omiten interacciones vol_5yr_std_gdp.")
  }
  
  # (b) Sólo si existe shock
  if ("shock" %in% names(df_cntl21)) {
    df_cntl21 <- df_cntl21 %>%
      mutate(
        vol_5yr_std_shock = vol_5yr_std * shock
      )
  } else {
    stop("Falta la variable 'shock'; no se pueden crear interacciones vol_5yr_std_shock.")
  }
  
} else {
  message("Interacciones de volatilidad ya existen; omitiendo.")
}

# 5) Dinámicas acumuladas
vars21 <- paste0("cumF", 0:12, "_dlog_capital")

if (!all(vars21 %in% names(df_cntl21))) {
  df_dyn21 <- df_cntl21 %>%
    group_by(name) %>%
    arrange(dateq, .by_group = TRUE) %>%
    group_modify(~{
      tmp <- .x
      for (h in 0:12) {
        tmp[[ paste0("cumF", h, "_dlog_capital") ]] <-
          rowSums(map_dfc(0:h, ~ lead(tmp$dlog_capital, .x)), na.rm = TRUE)
      }
      tmp
    }) %>%
    ungroup()
} else {
  message("Variables dinámicas ya existen; omitiendo.")
  df_dyn21 <- df_cntl21
}


library(dplyr)
library(purrr)
library(fixest)

# 6) Definir controles firm-level y macro (solo los macro presentes)
controls_firm21  <- c("rsales_g_std", "liq_std")
controls_macro21 <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")
present_macro21  <- intersect(controls_macro21, names(df_dyn21))
if (length(present_macro21) < length(controls_macro21)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro21, present_macro21), collapse = ", "),
    ". Se omitirán."
  )
}
controls21_terms <- c(controls_firm21, present_macro21)

# 7) Modelos vol_5yr_std_shock (ajustado)
has_vol_gdp <- "vol_5yr_std_gdp" %in% names(df_dyn21)

res21 <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  rhs_terms <- c(
    if (has_vol_gdp) "vol_5yr_std_gdp",
    "vol_5yr_std_shock",
    controls21_terms
  )
  fml_str <- paste0(
    dep_var, " ~ ",
    paste(rhs_terms, collapse = " + "),
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn21,
    cluster = ~ name + Country
  )
})

# 8) Extraer coeficientes y errores estándar
tbl21 <- tibble(
  horizon        = 0:12,
  beta_vol_shock = map_dbl(res21, ~ coef(.x)["vol_5yr_std_shock"]),
  se_vol_shock   = map_dbl(res21, ~ sqrt(vcov(.x)["vol_5yr_std_shock","vol_5yr_std_shock"]))
)

# Ahora puedes graficar tbl21 sin errores


# 9) Gráfico con leyenda
plot21 <- tbl21 %>%
  mutate(Serie="Volatilidad 5y × Shock") %>%
  ggplot(aes(x=horizon,y=beta_vol_shock,color=Serie,fill=Serie)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=beta_vol_shock-1.96*se_vol_shock,ymax=beta_vol_shock+1.96*se_vol_shock),alpha=0.2) +
  scale_color_manual(values=c("Volatilidad 5y × Shock"="purple"), name="Serie") +
  scale_fill_manual(values=c("Volatilidad 5y × Shock"="purple"), name="Serie") +
  labs(title="Figura 21: Volatilidad 5‑Year × Shock", x="Horizonte", y="Efecto acumulado de inversión") +
  scale_x_continuous(breaks=0:12) + theme_minimal()

print(plot21)



# ------------------------------------------------------------
# Figura 22 (adaptada): Dinámica Conjunta de Posición Financiera y Volatilidad 5-Year
# (toda la muestra, h = 0…12)
# ------------------------------------------------------------

df_cntl22 <- df

# 1) dlog_capital y Ldl_capital
if (!all(c("dlog_capital","Ldl_capital") %in% names(df_cntl22))) {
  df_cntl22 <- df_cntl22 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      capital2     = if_else(!is.na(capital)&capital>0, capital, NA_real_),
      log_capital  = log(capital2),
      dlog_capital = log_capital - lag(log_capital),
      Ldl_capital  = lag(dlog_capital)
    ) %>% ungroup() %>% select(-capital2,-log_capital)
} else message("dlog_capital y Ldl_capital ya existen; omitiendo.")

# 2) Controles firm-level
if (!"rsales_g_std" %in% names(df_cntl22)) {
  df_cntl22 <- df_cntl22 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      rsales_g     = log(saleq) - log(lag(saleq)),
      rsales_g_std = (rsales_g - mean(rsales_g,na.rm=TRUE))/sd(rsales_g,na.rm=TRUE)
    ) %>% ungroup() %>% select(-rsales_g)
} else message("rsales_g_std ya existe; omitiendo.")
if (!"liq_std" %in% names(df_cntl22)) {
  df_cntl22 <- df_cntl22 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(liq_std = (current_ratio - mean(current_ratio,na.rm=TRUE))/sd(current_ratio,na.rm=TRUE)) %>% ungroup()
} else message("liq_std ya existe; omitiendo.")
if (!all(c("lev_std","dd_std") %in% names(df_cntl22))) {
  df_cntl22 <- df_cntl22 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      lev_std = (leverage - mean(leverage,na.rm=TRUE))/sd(leverage,na.rm=TRUE),
      dd_std  = (dd        - mean(dd,       na.rm=TRUE))/sd(dd,       na.rm=TRUE)
    ) %>% ungroup()
} else message("lev_std y dd_std ya existen; omitiendo.")

# 3) Volatilidad 5-Year
if (!"vol_5yr_std" %in% names(df_cntl22)) {
  df_cntl22 <- df_cntl22 %>% group_by(name) %>% arrange(dateq) %>%
    mutate(
      yoy_sales    = log(saleq) - lag(log(saleq),4),
      vol_5yr      = zoo::rollapply(yoy_sales,20,sd,na.rm=TRUE,fill=NA,align="right"),
      vol_5yr_std  = (vol_5yr - mean(vol_5yr,na.rm=TRUE))/sd(vol_5yr,na.rm=TRUE)
    ) %>% ungroup() %>% select(-yoy_sales,-vol_5yr)
} else message("vol_5yr_std ya existe; omitiendo.")


# 4) Interacciones necesarias (ajustado)
# ---------------------------------------
# Solo creamos cada interacción si la(s) variable(s) base existen

# (a) Interacciones con dlog_gdp
if ("dlog_gdp" %in% names(df_cntl22)) {
  df_cntl22 <- df_cntl22 %>%
    mutate(
      vol_5yr_std_gdp = if (!"vol_5yr_std_gdp" %in% names(.)) vol_5yr_std * dlog_gdp else vol_5yr_std_gdp
    )
} else {
  warning("No existe 'dlog_gdp'; se omiten interacciones vol_5yr_std_gdp.")
}

# (b) Interacciones con shock
if ("shock" %in% names(df_cntl22)) {
  df_cntl22 <- df_cntl22 %>%
    mutate(
      vol_5yr_std_shock  = if (!"vol_5yr_std_shock"  %in% names(.)) vol_5yr_std * shock  else vol_5yr_std_shock,
      lev_wins_std_shock = if (!"lev_wins_std_shock" %in% names(.)) lev_std   * shock  else lev_wins_std_shock,
      dd_wins_std_shock  = if (!"dd_wins_std_shock"  %in% names(.)) dd_std    * shock  else dd_wins_std_shock
    )
} else {
  stop("Falta la variable 'shock'; no se pueden crear interacciones con shock.")
}



# 5) Dinámicas acumuladas
vars22 <- paste0("cumF",0:12,"_dlog_capital")
if (!all(vars22 %in% names(df_cntl22))) {
  df_dyn22 <- df_cntl22 %>% group_by(name) %>% arrange(dateq) %>%
    group_modify(~{tmp<-.x;for(h in 0:12) tmp[[paste0("cumF",h,"_dlog_capital")]]<-rowSums(map_dfc(0:h,~lead(tmp$dlog_capital,.x)),na.rm=TRUE);tmp})%>%ungroup()
} else { message("Variables dinámicas ya existen; omitiendo."); df_dyn22<-df_cntl22 }

# 6) Definir controles firm-level y macro (solo los macro presentes)
controls_firm22  <- c("rsales_g_std", "liq_std", "Ldl_capital")
controls_macro22 <- c("dlog_gdp", "dlog_cpi", "unemp", "embigl")

present_macro22 <- intersect(controls_macro22, names(df_dyn22))
if (length(present_macro22) < length(controls_macro22)) {
  warning(
    "Faltan controles macro: ",
    paste(setdiff(controls_macro22, present_macro22), collapse = ", "),
    ". Se omitirán."
  )
}

auto_controls22 <- paste(c(controls_firm22, present_macro22), collapse = " + ")

# 7a) Modelos DD & Vol 5y (ajustado)
res_dd22 <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var,
    " ~ dd_wins_std_shock + vol_5yr_std_shock + ",
    auto_controls22,
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn22,
    cluster = ~ name + Country
  )
})

# 7b) Modelos Leverage & Vol 5y (ajustado)
res_lev22 <- map(0:12, function(h) {
  dep_var <- paste0("cumF", h, "_dlog_capital")
  fml_str <- paste0(
    dep_var,
    " ~ lev_wins_std_shock + vol_5yr_std_shock + ",
    auto_controls22,
    " | name + Country"
  )
  feols(
    as.formula(fml_str),
    data    = df_dyn22,
    cluster = ~ name + Country
  )
})


# 8) Extraer coefs
tbl_dd22 <- tibble(
  horizon   = 0:12,
  beta_dd   = map_dbl(res_dd22, ~coef(.x)["dd_wins_std_shock"]),
  se_dd     = map_dbl(res_dd22, ~sqrt(vcov(.x)["dd_wins_std_shock","dd_wins_std_shock"])),
  beta_vol  = map_dbl(res_dd22, ~coef(.x)["vol_5yr_std_shock"]),
  se_vol    = map_dbl(res_dd22, ~sqrt(vcov(.x)["vol_5yr_std_shock","vol_5yr_std_shock"]))
)
tbl_lev22 <- tibble(
  horizon     = 0:12,
  beta_lev    = map_dbl(res_lev22, ~coef(.x)["lev_wins_std_shock"]),
  se_lev      = map_dbl(res_lev22, ~sqrt(vcov(.x)["lev_wins_std_shock","lev_wins_std_shock"])),
  beta_vol    = map_dbl(res_lev22, ~coef(.x)["vol_5yr_std_shock"]),
  se_vol      = map_dbl(res_lev22, ~sqrt(vcov(.x)["vol_5yr_std_shock","vol_5yr_std_shock"]))
)

# 9a) Gráfico DD & Volatilidad
long_dd22 <- tbl_dd22 %>% select(horizon, DD=beta_dd, Vol=beta_vol) %>% pivot_longer(-horizon,names_to="serie",values_to="beta")
se_dd22   <- tbl_dd22 %>% select(horizon, DD=se_dd, Vol=se_vol)     %>% pivot_longer(-horizon,names_to="serie",values_to="se")
dd22_plot <- left_join(long_dd22,se_dd22,by=c("horizon","serie"))

p22a <- ggplot(dd22_plot,aes(x=horizon,y=beta,color=serie,fill=serie))+
  geom_line(size=1)+geom_ribbon(aes(ymin=beta-1.96*se,ymax=beta+1.96*se),alpha=0.2)+
  scale_color_manual(values=c("DD"="red","Vol"="blue"),name="Serie")+
  scale_fill_manual(values=c("DD"="red","Vol"="blue"),name="Serie")+
  labs(title="Figura 22(a): DD vs Volatilidad 5y",x="Horizonte",y="Efecto acumulado de inversión")+
  theme_minimal()

# 9b) Gráfico Leverage & Volatilidad
long_lev22 <- tbl_lev22 %>% select(horizon, Lev=beta_lev, Vol=beta_vol) %>% pivot_longer(-horizon,names_to="serie",values_to="beta")
se_lev22   <- tbl_lev22 %>% select(horizon, Lev=se_lev, Vol=se_vol)     %>% pivot_longer(-horizon,names_to="serie",values_to="se")
lev22_plot <- left_join(long_lev22,se_lev22,by=c("horizon","serie"))

p22b <- ggplot(lev22_plot,aes(x=horizon,y=beta,color=serie,fill=serie))+
  geom_line(size=1)+geom_ribbon(aes(ymin=beta-1.96*se,ymax=beta+1.96*se),alpha=0.2)+
  scale_color_manual(values=c("Lev"="darkgreen","Vol"="blue"),name="Serie")+
  scale_fill_manual(values=c("Lev"="darkgreen","Vol"="blue"),name="Serie")+
  labs(title="Figura 22(b): Leverage vs Volatilidad 5y",x="Horizonte",y="Efecto acumulado de inversión")+
  theme_minimal()

# 10) Mostrar ambos
graph22 <- p22a / p22b + plot_annotation(title="Figura 22: Posición Financiera y Volatilidad 5-Year")
print(graph22)



  