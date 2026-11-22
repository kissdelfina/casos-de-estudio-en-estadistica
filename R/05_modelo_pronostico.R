# ==============================================================================
# PRONÓSTICO DE LA MÁXIMA INTERIOR: CON Y SIN LA TMAX DEL MISMO DÍA
# ==============================================================================
# Pregunta: ¿cuánto se pierde al pronosticar la máxima interior de hoy usando
# solo información exterior del pasado (sin la TMAX de hoy, que a la mañana
# todavía no se conoce)?
#
# Para cada serie (sala + subsector) se comparan tres modelos GLS con AR(1):
#
#   con_hoy      temp_max ~ TMAX_l0 + ... + TMAX_lK + sin_year + cos_year
#                Usa la TMAX de hoy: es el techo, no se puede usar para
#                pronosticar en la práctica.
#   sin_hoy      temp_max ~ TMAX_l1 + ... + TMAX_lK + sin_year + cos_year
#                Solo usa la TMAX de ayer hacia atrás: se puede usar a la mañana.
#   sin_exterior temp_max ~ sin_year + cos_year
#                Referencia: solo la época del año, sin datos del exterior.
#
# Evaluación fuera de la muestra, con origen móvil por bloques:
#   - Se entrena con los primeros días de la serie (la mitad, mínimo 40).
#   - Se predicen los 30 días siguientes, se agregan al entrenamiento,
#     se reajusta y se repite hasta el final.
#   - Las predicciones no usan la temperatura interior del período que se
#     predice (predict() de gls devuelve solo X·beta), así que predecir un
#     bloque entero equivale a pronosticar día por día.
#   - La ventana K de cada modelo se elige por AIC (ML) usando solo el primer
#     período de entrenamiento, para que la elección no mire el futuro.
#
# Métricas: RMSE y MAE de las predicciones fuera de la muestra, en °C.
#
# CORRECCIÓN CON EL ERROR DE LOS DÍAS ANTERIORES (modelos "_corr")
#   predict() de gls devuelve solo X·beta y no usa la correlación de los
#   errores. Si se conoce la máxima interior de los días anteriores (por
#   ejemplo, con un datalogger que se lee todos los días), la mejor
#   predicción lineal del modelo es
#       predicho_t = X_t·beta + r' R^-1 e
#   donde e son los errores observados (temp_max - X·beta) de los últimos
#   30 días con datos, R su matriz de correlación y r la correlación de cada
#   uno con el día t, según la estructura de error estimada en el
#   entrenamiento. Con AR(1) y días consecutivos se reduce a
#       predicho_t = X_t·beta + phi · e_(t-1).
#   Se usa el AR(1) de los mismos ajustes de arriba. (Se probó también con
#   ARMA(1,1): mejoraba entre 0 y 9 % y no cambiaba ninguna conclusión, así
#   que se sacó para que el script tarde la mitad.) Se pronostica día por
#   día: para el día t solo se usan errores de días anteriores a t.
#
# Igual que modelo_rezagos.R: se descartan el primer y el último día de cada
# archivo, los rezagos se calculan sobre la serie completa del SMN, y si algo
# falla se avisa. MACN-Or se separa en Antigua y Azara, como en
# 01_modelo_temperatura.R.
# ==============================================================================

library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)

max_rezago   <- 14   # ventana máxima que se prueba (días hacia atrás)
bloque       <- 30   # días que se predicen antes de reajustar
memoria      <- 30   # días hacia atrás que usa la corrección con errores pasados
min_entrena  <- 40   # mínimo de días del primer entrenamiento

# Rutas: funcionan con la estructura nueva (datos/...) y con la vieja
ruta <- function(nueva, vieja) if (file.exists(nueva)) nueva else vieja
ruta_xlsx  <- ruta("datos/smn/smn_datos_meteorologicos_1991_2020.xlsx",
                   "Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx")
ruta_lst   <- ruta("datos/smn/smn_datos_meteorologicos_desde_2021.lst",
                   "Datos TP FINAL/smn_datos_meteorologicos_desde_2021.lst")
ruta_hobo  <- ruta("datos/hobo", "Datos TP FINAL/No HOBO")
dir_tablas <- if (dir.exists("resultados")) "resultados/tablas" else "graficos"
dir_figs   <- if (dir.exists("resultados")) "resultados/figuras" else "graficos"
dir.create(dir_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. SMN: TMAX diaria de la estación 87585 y sus rezagos
# ------------------------------------------------------------------------------
# El xlsx (72 MB) consume mucha memoria al leerse: se lee una sola vez y se
# guarda la TMAX de la estación 87585 en un .rds chico.
ruta_cache <- if (dir.exists("datos")) "datos/procesados/smn_87585_1991_2020.rds" else
  file.path(tempdir(), "smn_87585_1991_2020.rds")
if (file.exists(ruta_cache)) {
  smn_1991 <- readRDS(ruta_cache)
} else {
  smn_1991 <- read_xlsx(ruta_xlsx) %>%
    filter(NRO_OMM == 87585) %>%
    transmute(Date = as.Date(FECHA),
              TMAX = suppressWarnings(as.numeric(gsub(",", ".", TMAX))))
  dir.create(dirname(ruta_cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(smn_1991, ruta_cache)
  gc()
}

smn_2021 <- read.delim(ruta_lst) %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
            TMAX = suppressWarnings(as.numeric(trimws(TMAX))))

smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date)

for (k in 0:max_rezago) {
  smn[[paste0("TMAX_l", k)]] <- dplyr::lag(smn$TMAX, k)
}

# ------------------------------------------------------------------------------
# 2. HOBO: máxima diaria por archivo
# ------------------------------------------------------------------------------
# Misma lógica que 01_modelo_temperatura.R (incluye la separación de MACN-Or)
extraer_etiqueta <- function(filepath) {
  fname <- basename(filepath)
  limpio <- gsub("[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}", "", fname)
  limpio <- gsub("\\.txt$", "", limpio)

  subsector <- "General"
  if (grepl("Anexo", limpio, ignore.case = TRUE)) subsector <- "Anexo"
  if (grepl("Ameghino", limpio, ignore.case = TRUE)) subsector <- "Ameghino"
  if (grepl("coleccion", limpio, ignore.case = TRUE)) subsector <- "Colección"
  if (grepl("laboratorio", limpio, ignore.case = TRUE)) subsector <- "Laboratorio"
  if (grepl("abajo", limpio, ignore.case = TRUE)) subsector <- "Abajo"
  if (grepl("arriba", limpio, ignore.case = TRUE)) subsector <- "Arriba"
  if (grepl("azara", limpio, ignore.case = TRUE) || grepl("04-01_2021", fname)) subsector <- "Azara"
  if (grepl("MACN-Or-03_02_2020", fname)) subsector <- "Antigua"

  sala <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba|azara|antigua)", "",
               limpio, ignore.case = TRUE)
  sala <- gsub("_", "-", sala)
  paste(sala, subsector)
}

leer_hobo <- function(path) {
  d <- read.csv2(path, fileEncoding = "UTF-8-BOM")
  tibble(Date = dmy(d[[1]]), temp = as.numeric(d[[3]])) %>%
    filter(!is.na(Date), !is.na(temp)) %>%
    group_by(Date) %>%
    summarise(temp_max = max(temp), .groups = "drop") %>%
    filter(Date > min(Date), Date < max(Date)) %>%   # primer y último día incompletos
    mutate(Archivo = basename(path), Serie = extraer_etiqueta(path))
}

archivos <- list.files(ruta_hobo, pattern = "\\.txt$", full.names = TRUE)
datos_hobo <- map_dfr(archivos, leer_hobo)

# ------------------------------------------------------------------------------
# 3. Funciones de ajuste y evaluación
# ------------------------------------------------------------------------------
formula_rezagos <- function(desde, hasta) {
  reformulate(c(paste0("TMAX_l", desde:hasta), "sin_year", "cos_year"),
              response = "temp_max")
}

ajustar <- function(f, datos) {
  tryCatch(gls(f, data = datos, correlation = corAR1(form = ~ t), method = "ML"),
           error = function(e) NULL)
}

# Correlación de los errores a distancia h = 0, 1, ..., memoria según el modelo
# (Si la serie tiene huecos, nlme guarda el corAR1 como corARMA(p = 1, q = 0):
# por eso se toman phi y theta por su nombre y no por la clase.)
correlaciones <- function(m) {
  p   <- coef(m$modelStruct$corStruct, unconstrained = FALSE)
  phi <- p[grepl("^Phi", names(p))]
  tht <- p[grepl("^Theta", names(p))]
  as.numeric(ARMAacf(ar = phi, ma = tht, lag.max = memoria))
}

# Predicción corregida para las filas idx_prueba:
#   X·beta + r' R^-1 e, con los errores de días anteriores a cada día
#   (como mucho 'memoria' días atrás). Si no hay días anteriores cercanos,
#   la corrección es 0.
predecir_corregido <- function(m, df, idx_prueba) {
  xb  <- as.numeric(predict(m, newdata = df))
  err <- df$temp_max - xb
  rho <- correlaciones(m)
  corr <- map_dbl(idx_prueba, function(i) {
    previos <- which(df$t < df$t[i] & df$t >= df$t[i] - memoria)
    if (length(previos) == 0) return(0)
    R <- outer(df$t[previos], df$t[previos], function(a, b) rho[abs(a - b) + 1])
    r <- rho[df$t[i] - df$t[previos] + 1]
    sum(r * solve(R, err[previos]))
  })
  xb[idx_prueba] + corr
}

# Elige la ventana K por AIC en los datos de entrenamiento
elegir_k <- function(desde, datos) {
  ks <- max(desde, 1):max_rezago
  if (desde == 0) ks <- 0:max_rezago
  aics <- map_dbl(ks, function(k) {
    m <- ajustar(formula_rezagos(desde, k), datos)
    if (is.null(m)) Inf else AIC(m)
  })
  ks[which.min(aics)]
}

evaluar_serie <- function(df_serie) {
  serie <- df_serie$Serie[1]

  if (anyDuplicated(df_serie$Date) > 0) {
    message("AVISO: '", serie, "' tiene fechas repetidas. Se saltea.")
    return(NULL)
  }

  df <- df_serie %>%
    left_join(smn, by = "Date") %>%
    mutate(sin_year = sin(2 * pi * yday(Date) / 365.25),
           cos_year = cos(2 * pi * yday(Date) / 365.25)) %>%
    drop_na(temp_max, all_of(paste0("TMAX_l", 0:max_rezago))) %>%
    arrange(Date) %>%
    mutate(t = as.integer(Date - min(Date)))

  n <- nrow(df)
  n_ini <- max(min_entrena, floor(n / 2))
  if (n - n_ini < 10) {
    message("AVISO: '", serie, "' tiene solo ", n, " días. Se saltea.")
    return(NULL)
  }

  # Ventana elegida con el primer período de entrenamiento
  entrena_ini <- df[1:n_ini, ]
  k_con <- elegir_k(0, entrena_ini)
  k_sin <- elegir_k(1, entrena_ini)

  formulas <- list(
    con_hoy      = formula_rezagos(0, k_con),
    sin_hoy      = formula_rezagos(1, k_sin),
    sin_exterior = temp_max ~ sin_year + cos_year
  )

  # Origen móvil: entrenar hasta el corte, predecir el bloque siguiente
  cortes <- seq(n_ini, n - 1, by = bloque)
  predicciones <- map_dfr(cortes, function(corte) {
    entrena <- df[1:corte, ]
    prueba  <- df[(corte + 1):min(corte + bloque, n), ]
    # Filas hasta el final del bloque: la corrección usa los errores de los
    # días anteriores, que dentro del bloque ya se conocen al pronosticar cada día
    hasta    <- df[1:min(corte + bloque, n), ]
    idx      <- (corte + 1):min(corte + bloque, n)
    map_dfr(names(formulas), function(nombre) {
      m <- ajustar(formulas[[nombre]], entrena)
      if (is.null(m)) {
        message("AVISO: '", serie, "', modelo ", nombre, ", corte ", corte, ": falló el ajuste.")
        return(NULL)
      }
      bind_rows(
        tibble(Modelo = nombre, predicho = as.numeric(predict(m, newdata = prueba))),
        tibble(Modelo = paste0(nombre, "_corr"), predicho = predecir_corregido(m, hasta, idx))
      ) %>%
        mutate(Serie = serie, Date = rep(prueba$Date, length.out = n()),
               observado = rep(prueba$temp_max, length.out = n()), .before = 1)
    })
  })

  list(
    predicciones = predicciones,
    info = tibble(Serie = serie, N_dias = n, Dias_entrena_inicial = n_ini,
                  Dias_evaluados = n - n_ini, K_con_hoy = k_con, K_sin_hoy = k_sin)
  )
}

# ------------------------------------------------------------------------------
# 4. Correr para todas las series
# ------------------------------------------------------------------------------
series     <- split(datos_hobo, datos_hobo$Serie)
resultados <- map(series, evaluar_serie) %>% compact()

predicciones <- map_dfr(resultados, "predicciones")
info         <- map_dfr(resultados, "info")

errores <- predicciones %>%
  mutate(error = predicho - observado) %>%
  group_by(Serie, Modelo) %>%
  summarise(RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)), .groups = "drop")

resumen <- errores %>%
  select(Serie, Modelo, RMSE) %>%
  pivot_wider(names_from = Modelo, values_from = RMSE, names_prefix = "RMSE_") %>%
  left_join(info, by = "Serie") %>%
  mutate(
    # Cuánto empeora el pronóstico al no tener la TMAX de hoy
    Perdida_sin_hoy_pct = 100 * (RMSE_sin_hoy / RMSE_con_hoy - 1),
    # Cuánto mejora el exterior del pasado frente a solo la época del año
    Mejora_vs_sin_ext_pct = 100 * (1 - RMSE_sin_hoy / RMSE_sin_exterior)
  ) %>%
  select(Serie, N_dias, Dias_evaluados, K_con_hoy, K_sin_hoy,
         RMSE_con_hoy, RMSE_sin_hoy, RMSE_sin_exterior,
         Perdida_sin_hoy_pct, Mejora_vs_sin_ext_pct) %>%
  arrange(desc(Perdida_sin_hoy_pct))

cat("\n=== PRONÓSTICO FUERA DE LA MUESTRA (RMSE en °C) ===\n")
print(resumen %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf, width = Inf)

write_csv(resumen, file.path(dir_tablas, "pronostico_con_y_sin_hoy.csv"))
write_csv(errores, file.path(dir_tablas, "pronostico_errores_por_modelo.csv"))

# Corrección con el error de los días anteriores
resumen_corr <- errores %>%
  select(Serie, Modelo, RMSE) %>%
  pivot_wider(names_from = Modelo, values_from = RMSE, names_prefix = "RMSE_") %>%
  left_join(info %>% select(Serie, Dias_evaluados), by = "Serie") %>%
  mutate(
    # Cuánto baja el error al corregir con los días anteriores (AR(1))
    Mejora_corr_sin_hoy_pct = 100 * (1 - RMSE_sin_hoy_corr / RMSE_sin_hoy),
    # Con corrección: cuánto se pierde sin la TMAX de hoy
    Perdida_sin_hoy_corr_pct = 100 * (RMSE_sin_hoy_corr / RMSE_con_hoy_corr - 1),
    # Con corrección: cuánto aporta el exterior además de la inercia de la sala
    Aporte_exterior_corr_pct = 100 * (1 - RMSE_sin_hoy_corr / RMSE_sin_exterior_corr)
  ) %>%
  select(Serie, Dias_evaluados, RMSE_sin_exterior, RMSE_sin_hoy, RMSE_con_hoy,
         RMSE_sin_exterior_corr, RMSE_sin_hoy_corr, RMSE_con_hoy_corr,
         Mejora_corr_sin_hoy_pct, Perdida_sin_hoy_corr_pct,
         Aporte_exterior_corr_pct) %>%
  arrange(desc(Mejora_corr_sin_hoy_pct))

cat("\n=== CORRECCIÓN CON EL ERROR DE LOS DÍAS ANTERIORES (RMSE en °C) ===\n")
print(resumen_corr %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf, width = Inf)
write_csv(resumen_corr, file.path(dir_tablas, "pronostico_correccion.csv"))

# ------------------------------------------------------------------------------
# 5. Gráfico para la presentación
# ------------------------------------------------------------------------------
# Solo las series con al menos 60 días evaluados: en MACN-Ar (38) y MACN-In (12)
# el entrenamiento tiene menos de 40 días, el seno y el coseno extrapolan mal, y
# el RMSE no es comparable. Quedan en la tabla, pero no en el gráfico.
series_graf <- resumen %>% filter(Dias_evaluados >= 60) %>%
  arrange(Perdida_sin_hoy_pct) %>% pull(Serie)

datos_graf <- errores %>%
  filter(Serie %in% series_graf, Modelo %in% c("sin_exterior", "sin_hoy", "con_hoy")) %>%
  mutate(Serie = factor(Serie, levels = series_graf),
         Modelo = factor(Modelo, levels = c("sin_exterior", "sin_hoy", "con_hoy"),
                         labels = c("Solo época del año", "Sin la TMAX de hoy",
                                    "Con la TMAX de hoy")))

g_rmse <- ggplot(datos_graf, aes(x = RMSE, y = Serie)) +
  geom_line(aes(group = Serie), color = "gray80", linewidth = 1) +
  geom_point(aes(color = Modelo), size = 3.5) +
  scale_color_manual(values = c("Solo época del año" = "gray70",
                                "Sin la TMAX de hoy" = "#2E6FD8",
                                "Con la TMAX de hoy" = "gray25")) +
  labs(
    title = "Sin la TMAX de hoy, el error sube mucho en las 2 salas rápidas y poco en el resto",
    subtitle = "Error de pronóstico fuera de la muestra (origen móvil, bloques de 30 días).\nOrdenado por cuánto empeora el pronóstico sin la TMAX de hoy.",
    x = "RMSE de la máxima interior (°C)", y = NULL, color = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title.position = "plot")

print(g_rmse)
ggsave(file.path(dir_figs, "fig_pronostico_rmse.png"), g_rmse, width = 13, height = 7.5, dpi = 300)

# ------------------------------------------------------------------------------
# 6. Gráfico de la corrección con el error de los días anteriores
# ------------------------------------------------------------------------------
modelos_corr <- c("sin_hoy" = "Sin la TMAX de hoy",
                  "sin_exterior_corr" = "Solo época del año + días anteriores de la sala",
                  "sin_hoy_corr" = "Sin la TMAX de hoy + días anteriores de la sala",
                  "con_hoy_corr" = "Con la TMAX de hoy + días anteriores de la sala")
orden_corr <- resumen_corr %>% filter(Serie %in% series_graf) %>%
  arrange(RMSE_sin_hoy_corr) %>% pull(Serie)

datos_corr <- errores %>%
  filter(Serie %in% series_graf, Modelo %in% names(modelos_corr)) %>%
  mutate(Serie = factor(Serie, levels = rev(orden_corr)),
         Modelo = factor(Modelo, levels = names(modelos_corr), labels = modelos_corr))

g_corr <- ggplot(datos_corr, aes(x = RMSE, y = Serie)) +
  geom_line(aes(group = Serie), color = "gray80", linewidth = 1) +
  geom_point(aes(color = Modelo), size = 3.5) +
  scale_color_manual(values = setNames(c("gray70", "#D2692E", "#2E6FD8", "gray25"), modelos_corr)) +
  guides(color = guide_legend(ncol = 2)) +
  labs(
    title = "Conocer la temperatura de la sala de los días anteriores baja mucho el error",
    subtitle = "Error de pronóstico fuera de la muestra (origen móvil, bloques de 30 días), corrigiendo con\nlos errores de los días anteriores según el AR(1) estimado.",
    x = "RMSE de la máxima interior (°C)", y = NULL, color = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title.position = "plot")

print(g_corr)
ggsave(file.path(dir_figs, "fig_pronostico_correccion.png"), g_corr, width = 13, height = 8, dpi = 300)
