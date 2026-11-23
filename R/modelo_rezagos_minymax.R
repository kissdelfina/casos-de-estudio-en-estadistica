# ==============================================================================
# MODELO CON REZAGOS SEPARADOS DE LA TMAX EXTERIOR
# ==============================================================================
# En vez del promedio de 3 días (TMAX_prom3d) de modelo_temperatura.R, cada día
# anterior entra como predictora separada:
#
#   temp_max ~ TMAX_l0 + TMAX_l1 + ... + TMAX_lK + sin_year + cos_year,  AR(1)
#
# donde TMAX_lk es la máxima exterior de hace k días. Para cada serie:
#   - se prueban K = 0, 1, ..., max_rezago y se elige el K con menor AIC;
#   - se compara contra el modelo del promedio de 3 días y contra un modelo
#     sin exterior (para medir el aporte real del exterior);
#   - se testea si los pesos iguales del promedio son razonables: el promedio
#     de 3 días es el modelo de rezagos 0 a 2 con los 3 coeficientes iguales
#     (modelos anidados -> test de razón de verosimilitud);
#   - se informa el efecto acumulado (suma de los coeficientes de los rezagos)
#     con su error estándar. Los rezagos están muy correlacionados entre sí, así
#     que cada coeficiente por separado es inestable, pero la suma no.
#
# Todos los modelos se ajustan por ML (no REML) para que los AIC sean comparables,
# y todos usan exactamente las mismas filas dentro de cada serie.
#
# Diferencias con modelo_temperatura.R:
#   - Los rezagos se calculan sobre la serie completa del SMN, así que no se
#     pierden los primeros días de cada serie.
#   - Se descartan el primer y el último día de cada archivo (están incompletos).
#   - Los .txt se leen con read.csv2 (coma decimal) y ruta relativa al repo.
#   - Si una serie falla, se avisa en vez de descartarla en silencio.
# ==============================================================================

library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)

max_rezago <- 14  # cantidad máxima de días hacia atrás que se prueban

# ------------------------------------------------------------------------------
# 1. SMN: TMAX y TMIN diaria de la estación 87585
# ------------------------------------------------------------------------------

# Cachés del xlsx 1991-2020 (72 MB, lento y pesado de leer): la TMAX es la
# misma que usan los otros scripts; la TMIN va en un caché aparte para no
# modificar el que leen ellos. Si falta alguno, se lee el xlsx una sola vez.
ruta_cache_tmax <- "datos/procesados/smn_87585_1991_2020.rds"
ruta_cache_tmin <- "datos/procesados/smn_87585_tmin_1991_2020.rds"

if (!file.exists(ruta_cache_tmax) || !file.exists(ruta_cache_tmin)) {
  smn_xlsx <- read_xlsx(
    "datos/smn/smn_datos_meteorologicos_1991_2020.xlsx"
  ) %>%
    filter(NRO_OMM == 87585) %>%
    transmute(
      Date = as.Date(FECHA),
      TMAX = suppressWarnings(as.numeric(gsub(",", ".", TMAX))),
      TMIN = suppressWarnings(as.numeric(gsub(",", ".", TMIN)))
    )
  dir.create(dirname(ruta_cache_tmax), recursive = TRUE,
             showWarnings = FALSE)
  if (!file.exists(ruta_cache_tmax)) saveRDS(select(smn_xlsx, Date, TMAX), ruta_cache_tmax)
  if (!file.exists(ruta_cache_tmin)) saveRDS(select(smn_xlsx, Date, TMIN), ruta_cache_tmin)
  rm(smn_xlsx)
  gc()
}

smn_1991 <- readRDS(ruta_cache_tmax) %>%
  left_join(readRDS(ruta_cache_tmin), by = "Date")

smn_2021 <- read.delim(
  "datos/smn/smn_datos_meteorologicos_desde_2021.lst"
) %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(
    Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
    TMAX = suppressWarnings(as.numeric(trimws(TMAX))),
    TMIN = suppressWarnings(as.numeric(trimws(TMIN)))
  )

# Grilla diaria completa
smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date)

# Rezagos de ambas variables
for (k in 0:max_rezago) {
  smn[[paste0("TMAX_l", k)]] <- dplyr::lag(smn$TMAX, k)
  smn[[paste0("TMIN_l", k)]] <- dplyr::lag(smn$TMIN, k)
}

# Promedios móviles de tres días
smn <- smn %>%
  mutate(
    TMAX_prom3d = (TMAX_l0 + TMAX_l1 + TMAX_l2) / 3,
    TMIN_prom3d = (TMIN_l0 + TMIN_l1 + TMIN_l2) / 3
  )
# ------------------------------------------------------------------------------
# 2. HOBO: temperatura máxima diaria por archivo
# ------------------------------------------------------------------------------
# Misma lógica de etiquetas que 03_modelo_rezagos.R, para que las series sean
# las mismas y los resultados se puedan comparar.
extraer_info_sala <- function(filepath) {
  fname <- basename(filepath)
  fname <- gsub("[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}", "", fname)
  fname <- gsub("\\.txt$", "", fname)
  
  subsector <- "General / Normal"
  if (grepl("Anexo", fname, ignore.case = TRUE)) subsector <- "Anexo"
  if (grepl("Ameghino", fname, ignore.case = TRUE)) subsector <- "Ameghino"
  if (grepl("coleccion", fname, ignore.case = TRUE)) subsector <- "Colección"
  if (grepl("laboratorio", fname, ignore.case = TRUE)) subsector <- "Laboratorio"
  if (grepl("abajo", fname, ignore.case = TRUE)) subsector <- "Abajo"
  if (grepl("arriba", fname, ignore.case = TRUE)) subsector <- "Arriba"
  
  # MACN-Or: son dos dataloggers con fechas superpuestas (colección antigua y
  # Azara), igual que en 03_modelo_rezagos.R
  if (grepl("azara", fname, ignore.case = TRUE) || grepl("04-01_2021", basename(filepath))) subsector <- "Azara"
  if (grepl("MACN-Or-03_02_2020", basename(filepath))) subsector <- "Antigua"
  
  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba|azara|antigua)", "", fname, ignore.case = TRUE)
  sala_base <- gsub("-", "_", sala_base)
  
  paste(sala_base, subsector, sep = " - ")
}

leer_hobo <- function(path) {
  d <- read.csv2(path, fileEncoding = "UTF-8-BOM")
  
  tibble(
    Date = dmy(d[[1]]),
    temp = as.numeric(d[[3]])
  ) %>%
    filter(!is.na(Date), !is.na(temp)) %>%
    group_by(Date) %>%
    summarise(
      temp_max = max(temp),
      temp_min = min(temp),
      .groups = "drop"
    ) %>%
    filter(Date > min(Date), Date < max(Date)) %>%
    mutate(
      Archivo = basename(path),
      Etiqueta = extraer_info_sala(path)
    )
}
archivos <- list.files("datos/hobo", pattern = "\\.txt$", full.names = TRUE)
datos_hobo <- map_dfr(archivos, leer_hobo)

# ------------------------------------------------------------------------------
# 3. Ajuste por serie y variable de temperatura
# ------------------------------------------------------------------------------

nombres_rezagos <- function(variable, k) {
  paste0(variable, "_l", 0:k)
}

# Configuraciones de los dos análisis
configuraciones <- list(
  list(
    Respuesta = "temp_max",
    Exterior = "TMAX"
  ),
  list(
    Respuesta = "temp_min",
    Exterior = "TMIN"
  )
)

ajustar_serie <- function(df_serie, respuesta, exterior) {
  
  etiqueta <- df_serie$Etiqueta[1]
  
  # Con fechas repetidas (dos dataloggers en la misma serie)
  if (anyDuplicated(df_serie$Date) > 0) {
    message(
      "AVISO: '", etiqueta,
      "' tiene fechas repetidas. Se saltea."
    )
    return(NULL)
  }
  
  # Nombres de las variables exteriores
  rezagos_max <- nombres_rezagos(exterior, max_rezago)
  prom3d <- paste0(exterior, "_prom3d")
  
  # Unir HOBO y SMN
  df <- df_serie %>%
    left_join(smn, by = "Date") %>%
    mutate(
      yday = yday(Date),
      sin_year = sin(2 * pi * yday / 365.25),
      cos_year = cos(2 * pi * yday / 365.25)
    ) %>%
    drop_na(
      all_of(c(respuesta, rezagos_max, prom3d))
    ) %>%
    arrange(Date) %>%
    mutate(t = as.integer(Date - min(Date)))
  
  if (nrow(df) < 30) {
    message(
      "AVISO: '", etiqueta, "' (", respuesta,
      ") tiene solo ", nrow(df), " días completos. Se saltea."
    )
    return(NULL)
  }
  
  # Todos los modelos utilizan exactamente las mismas filas
  formulas <- list(
    sin_exterior = reformulate(
      c("sin_year", "cos_year"),
      response = respuesta
    ),
    prom3d = reformulate(
      c(prom3d, "sin_year", "cos_year"),
      response = respuesta
    )
  )
  
  # Modelos con rezagos separados
  for (k in 0:max_rezago) {
    formulas[[paste0("rezagos_0a", k)]] <-
      reformulate(
        c(nombres_rezagos(exterior, k),
          "sin_year", "cos_year"),
        response = respuesta
      )
  }
  
  # Ajuste de todos los modelos mediante ML y AR(1)
  modelos <- imap(formulas, function(f, nombre) {
    tryCatch(
      gls(
        f,
        data = df,
        correlation = corAR1(form = ~ t),
        method = "ML"
      ),
      error = function(e) {
        message(
          "AVISO: '", etiqueta, "', ", respuesta,
          ", modelo ", nombre, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  })
  
  if (any(map_lgl(modelos, is.null))) {
    message(
      "AVISO: '", etiqueta, "' (", respuesta,
      ") se saltea porque falló algún modelo."
    )
    return(NULL)
  }
  
  # --------------------------------------------------------------------------
  # Tabla de AIC
  # --------------------------------------------------------------------------
  
  aic <- tibble(
    Etiqueta = etiqueta,
    Respuesta = respuesta,
    Exterior = exterior,
    Modelo = names(modelos),
    AIC = map_dbl(modelos, AIC),
    phi = map_dbl(
      modelos,
      ~ coef(.x$modelStruct$corStruct, unconstrained = FALSE)
    )
  )
  
  # --------------------------------------------------------------------------
  # Mejor modelo según AIC
  # --------------------------------------------------------------------------
  
  aic_rez <- aic %>%
    filter(startsWith(Modelo, "rezagos"))
  
  mejor_nm <- aic_rez$Modelo[which.min(aic_rez$AIC)]
  k_mejor <- as.integer(sub("rezagos_0a", "", mejor_nm))
  
  mejor <- modelos[[mejor_nm]]
  
  # Coeficientes de los rezagos del mejor modelo
  b <- coef(mejor)
  V <- vcov(mejor)
  
  nm <- nombres_rezagos(exterior, k_mejor)
  
  coefs <- tibble(
    Etiqueta = etiqueta,
    Respuesta = respuesta,
    Exterior = exterior,
    Rezago = 0:k_mejor,
    Coef = unname(b[nm]),
    EE = unname(sqrt(diag(V)[nm]))
  )
  
  # Efecto acumulado
  efecto_acum <- sum(b[nm])
  ee_acum <- sqrt(sum(V[nm, nm]))
  
  # --------------------------------------------------------------------------
  # Test de pesos iguales
  # --------------------------------------------------------------------------
  
  # El promedio de 3 días está anidado en el modelo 0a2
  p_pesos_iguales <- anova(
    modelos$prom3d,
    modelos$rezagos_0a2
  )$`p-value`[2]
  
  # --------------------------------------------------------------------------
  # Resumen
  # --------------------------------------------------------------------------
  
  resumen <- tibble(
    Etiqueta = etiqueta,
    Respuesta = respuesta,
    Exterior = exterior,
    N_dias = nrow(df),
    
    Rezagos_mejor = paste0("0 a ", k_mejor),
    Efecto_acumulado = efecto_acum,
    EE_acumulado = ee_acum,
    
    Pendiente_prom3d = unname(coef(modelos$prom3d)[prom3d]),
    
    dAIC_vs_prom3d = AIC(mejor) - AIC(modelos$prom3d),
    dAIC_vs_sin_ext = AIC(mejor) - AIC(modelos$sin_exterior),
    
    p_pesos_iguales = p_pesos_iguales,
    
    cor_l0_l1 = cor(
      df[[paste0(exterior, "_l0")]],
      df[[paste0(exterior, "_l1")]]
    )
  )
  
  list(
    resumen = resumen,
    coefs = coefs,
    aic = aic
  )
}

# ------------------------------------------------------------------------------
# Ejecutar ambos análisis para todas las series
# ------------------------------------------------------------------------------

series <- split(datos_hobo, datos_hobo$Etiqueta)

resultados <- map(
  series,
  function(s) {
    map(
      configuraciones,
      function(cfg) {
        ajustar_serie(
          s,
          respuesta = cfg$Respuesta,
          exterior = cfg$Exterior
        )
      }
    )
  }
) %>%
  flatten() %>%
  compact()

# Unificar resultados
resumen_rezagos <- map_dfr(resultados, "resumen") %>%
  arrange(Respuesta, desc(Efecto_acumulado))

coefs_rezagos <- map_dfr(resultados, "coefs")

aic_rezagos <- map_dfr(resultados, "aic")

# ------------------------------------------------------------------------------
# 4. Resultados
# ------------------------------------------------------------------------------

cat("\n=== RESUMEN POR SERIE Y TEMPERATURA ===\n")

print(
  resumen_rezagos %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  n = Inf,
  width = Inf
)

cat("\n=== COEFICIENTES DE CADA REZAGO ===\n")

print(
  coefs_rezagos %>%
    mutate(across(c(Coef, EE), ~ round(.x, 3))),
  n = Inf
)

# Guardar tablas
dir.create(
  "resultados/tablas",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  resumen_rezagos,
  "resultados/tablas/rezagos_resumen_minymax.csv"
)

write_csv(
  coefs_rezagos,
  "resultados/tablas/rezagos_coeficientes_minymax.csv"
)

write_csv(
  aic_rezagos,
  "resultados/tablas/rezagos_aic_minymax.csv"
)

# ------------------------------------------------------------------------------
# Gráfico 1: coeficientes de cada rezago
# ------------------------------------------------------------------------------

g_rezagos <- ggplot(
  coefs_rezagos,
  aes(x = Rezago, y = Coef)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_pointrange(
    aes(
      ymin = Coef - 1.96 * EE,
      ymax = Coef + 1.96 * EE
    ),
    color = "firebrick"
  ) +
  geom_line(
    color = "firebrick",
    alpha = 0.5
  ) +
  facet_wrap(
    vars(Respuesta, Etiqueta),
    ncol = 3,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = 0:max_rezago
  ) +
  labs(
    title = "Efecto de la temperatura exterior sobre la temperatura interior",
    subtitle = "Mejor modelo por AIC (ML). Barras: IC 95 %",
    x = "Rezago (días)",
    y = "Coeficiente (°C interior por °C exterior)"
  ) +
  theme_minimal(base_size = 10)

print(g_rezagos)

# ------------------------------------------------------------------------------
# Gráfico 2: efecto acumulado vs. pendiente del promedio de 3 días
# ------------------------------------------------------------------------------

g_comparacion <- ggplot(
  resumen_rezagos,
  aes(
    x = Pendiente_prom3d,
    y = Efecto_acumulado,
    label = Etiqueta
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_errorbar(
    aes(
      ymin = Efecto_acumulado - 1.96 * EE_acumulado,
      ymax = Efecto_acumulado + 1.96 * EE_acumulado
    ),
    width = 0,
    color = "gray60"
  ) +
  geom_point(
    color = "steelblue",
    size = 2.5
  ) +
  geom_text(
    size = 2.8,
    vjust = -0.8,
    check_overlap = TRUE
  ) +
  facet_wrap(~ Respuesta, scales = "free") +
  labs(
    title = "Efecto acumulado de los rezagos vs. promedio de 3 días",
    x = "Pendiente del promedio de 3 días",
    y = "Suma de los coeficientes de los rezagos"
  ) +
  theme_minimal()

print(g_comparacion)

