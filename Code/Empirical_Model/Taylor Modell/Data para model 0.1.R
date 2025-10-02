# Cargar librerías necesarias
library(readxl)
library(tempdisagg)
library(writexl)
library(dplyr)

# Definir rutas
ruta_entrada <- "C:/Users/joser/Desktop/GDP.xlsx"
ruta_salida_folder <- "C:/Users/joser/Desktop/Trimestralizado 2/"

# Crear carpeta de salida si no existe
if (!dir.exists(ruta_salida_folder)) {
  dir.create(ruta_salida_folder, recursive = TRUE)
}

# Obtener nombres de las hojas del archivo
hojas <- excel_sheets(ruta_entrada)

# Recorrer cada hoja
for (hoja in hojas) {
  # Leer datos de la hoja actual
  datos <- read_excel(ruta_entrada, sheet = hoja)
  
  # Verificar que haya exactamente dos columnas
  if (ncol(datos) != 2) {
    stop(paste0("La hoja '", hoja, "' debe tener exactamente 2 columnas: Año y GDP."))
  }
  
  # Renombrar columnas (remover espacios en nombres)
  nombres_actuales <- names(datos)
  nombres_limpios <- trimws(nombres_actuales)
  names(datos) <- nombres_limpios
  
  # Asumir que la primera columna es el año numérico y la segunda es GDP numérico
  nombre_anio <- nombres_limpios[1]
  nombre_gdp  <- nombres_limpios[2]
  
  # Verificar que la primera columna sea numérica
  if (!is.numeric(datos[[nombre_anio]])) {
    stop(paste0("En la hoja '", hoja, "', la primera columna ('", nombre_anio,
                "') debe contener los años en formato numérico."))
  }
  # Verificar que la segunda columna sea numérica
  if (!is.numeric(datos[[nombre_gdp]])) {
    stop(paste0("En la hoja '", hoja, "', la segunda columna ('", nombre_gdp,
                "') debe contener valores numéricos de GDP."))
  }
  
  # Extraer año inicial
  anio_inicio <- datos[[nombre_anio]][1]
  
  # Convertir la columna de GDP a serie de tiempo anual
  valores_gdp <- datos[[nombre_gdp]]
  serie_gdp_ts <- ts(as.numeric(valores_gdp),
                     start = c(anio_inicio, 1),
                     frequency = 1)
  
  # Aplicar Denton–Cholette para trimestralizar la serie de GDP
  modelo_td <- td(
    formula = serie_gdp_ts ~ 1,
    to      = "quarterly",
    method  = "denton-cholette"
  )
  serie_gdp_trimestral <- predict(modelo_td)
  
  # Determinar número de periodos trimestrales y construir vector de fechas
  n_periodos <- length(serie_gdp_trimestral)
  fechas_trimestrales <- seq(
    from = as.Date(paste0(anio_inicio, "-01-01")),
    by   = "quarter",
    length.out = n_periodos
  )
  
  # Armado del data.frame resultante
  data_final <- data.frame(
    Periodo = fechas_trimestrales,
    GDP     = as.numeric(serie_gdp_trimestral)
  )
  
  # Definir ruta de salida para esta hoja (un archivo Excel por hoja)
  nombre_archivo_salida <- paste0(ruta_salida_folder, hoja, ".xlsx")
  
  # Guardar el data.frame trimestralizado en un archivo Excel
  write_xlsx(data_final, path = nombre_archivo_salida)
  
  cat("✅ Hoja '", hoja, "' trimestralizada y guardada en:\n   ", 
      nombre_archivo_salida, "\n", sep = "")
}

cat("\nTodas las hojas han sido procesadas.\n")



#----------------------------------------------------------------------------

# Cargar librerías necesarias
library(readxl)
library(tempdisagg)
library(writexl)
library(dplyr)

# Definir rutas
ruta_entrada     <- "C:/Users/joser/Desktop/Inflacion.xlsx"
ruta_salida_base <- "C:/Users/joser/Desktop/Trimestralizado 2/Inflacion/"

# Crear carpeta de salida si no existe
if (!dir.exists(ruta_salida_base)) {
  dir.create(ruta_salida_base, recursive = TRUE)
}

# Obtener nombres de las hojas del archivo
hojas <- excel_sheets(ruta_entrada)

# Recorrer cada hoja (país)
for (hoja in hojas) {
  # Leer datos de la hoja actual
  datos <- read_excel(ruta_entrada, sheet = hoja)
  
  # Verificar que haya exactamente 2 columnas: Año y Variación (%) de inflación
  if (ncol(datos) != 2) {
    stop(paste0("La hoja '", hoja, "' debe tener exactamente 2 columnas: Año e inflación (%)."))
  }
  
  # Limpiar nombres de columnas (quitar espacios en blanco al inicio/final)
  names(datos) <- trimws(names(datos))
  nombre_anio    <- names(datos)[1]
  nombre_infl_pct <- names(datos)[2]
  
  # Verificar tipos de datos
  if (!is.numeric(datos[[nombre_anio]])) {
    stop(paste0("En '", hoja, "', la columna '", nombre_anio, "' debe ser numérica (años)."))
  }
  if (!is.numeric(datos[[nombre_infl_pct]])) {
    stop(paste0("En '", hoja, "', la columna '", nombre_infl_pct, "' debe ser numérica (inflación en %)."))
  }
  
  # Extraer vector de años y de variaciones porcentuales anuales
  vector_anios     <- datos[[nombre_anio]]
  infl_pct_anual   <- datos[[nombre_infl_pct]]
  
  # 1) Convertir variaciones porcentuales anuales en índice de precios
  #    Suponer índice inicial = 100 en el primer año
  n_anios <- length(vector_anios)
  indice_anual <- numeric(n_anios)
  indice_anual[1] <- 100
  for (i in 2:n_anios) {
    indice_anual[i] <- indice_anual[i - 1] * (1 + infl_pct_anual[i] / 100)
  }
  # Ahora 'indice_anual' es un vector tipo índice (base 100 en primer año)
  
  # 2) Convertir ese índice anual a serie de tiempo (freq = 1)
  anio_inicio  <- vector_anios[1]
  ts_indice_anual <- ts(indice_anual, start = c(anio_inicio, 1), frequency = 1)
  
  # 3) Trimestralizar el índice anual usando Denton–Cholette
  modelo_td <- td(
    formula = ts_indice_anual ~ 1,
    to      = "quarterly",
    method  = "denton-cholette"
  )
  ts_indice_trimestral <- predict(modelo_td)
  # 'ts_indice_trimestral' ahora es un objeto ts de frecuencia 4
  
  # 4) Convertir el índice trimestral de vuelta a variaciones porcentuales (inflación trimestral)
  #    Inflación trimestral t = (Índice_t / Índice_{t-1} - 1) * 100
  valores_indice_q <- as.numeric(ts_indice_trimestral)
  n_periodos_q     <- length(valores_indice_q)
  infl_pct_trimestral <- numeric(n_periodos_q)
  infl_pct_trimestral[1] <- NA  # primer trimestre no tiene variación previa
  for (j in 2:n_periodos_q) {
    infl_pct_trimestral[j] <- (valores_indice_q[j] / valores_indice_q[j - 1] - 1) * 100
  }
  
  # 5) Construir vector de fechas trimestrales (al primer día de cada trimestre)
  fechas_trimestrales <- seq(
    from       = as.Date(paste0(anio_inicio, "-01-01")),
    by         = "quarter",
    length.out = n_periodos_q
  )
  
  # 6) Armar data.frame final con Periodo (fecha) e inflación trimestral (%)
  data_final <- data.frame(
    Periodo = fechas_trimestrales,
    Inflacion_pct = round(infl_pct_trimestral, 4)
  )
  
  # 7) Guardar el resultado en un archivo Excel (una hoja por país)
  nombre_archivo_salida <- paste0(ruta_salida_base, hoja, ".xlsx")
  write_xlsx(data_final, path = nombre_archivo_salida)
  
  cat("✅ '", hoja, "' procesada: índice anual → trimestralización → inflación trimestral guardada en:\n  ",
      nombre_archivo_salida, "\n", sep = "")
}

cat("\nTodas las hojas se han procesado correctamente.\n")
















