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

max_rezago <- 7   # cantidad máxima de días hacia atrás que se prueban

# ------------------------------------------------------------------------------
# 1. SMN: TMAX diaria de la estación 87585 (Buenos Aires Observatorio)
# ------------------------------------------------------------------------------
smn_1991 <- read_xlsx("Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx") %>%
  filter(NRO_OMM == 87585) %>%
  transmute(Date = as.Date(FECHA),
            TMAX = suppressWarnings(as.numeric(gsub(",", ".", TMAX))))  # "S/D" -> NA

smn_2021 <- read.delim("Datos TP FINAL/smn_datos_meteorologicos_desde_2021.lst") %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
            TMAX = suppressWarnings(as.numeric(trimws(TMAX))))

# Grilla diaria completa: así lag(TMAX, k) es siempre "hace k días", aunque
# falte algún día en el SMN.
smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date)

for (k in 0:max_rezago) {
  smn[[paste0("TMAX_l", k)]] <- dplyr::lag(smn$TMAX, k)
}
smn$TMAX_prom3d <- (smn$TMAX_l0 + smn$TMAX_l1 + smn$TMAX_l2) / 3

# ------------------------------------------------------------------------------
# 2. HOBO: temperatura máxima diaria por archivo
# ------------------------------------------------------------------------------
# Misma lógica de etiquetas que modelo_temperatura.R, para que las series sean
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

  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba)", "", fname, ignore.case = TRUE)
  sala_base <- gsub("-", "_", sala_base)

  paste(sala_base, subsector, sep = " - ")
}

leer_hobo <- function(path) {
  d <- read.csv2(path, fileEncoding = "UTF-8-BOM")
  # Por posición: columna 1 = fecha (dd/mm/aa), columna 3 = temperatura
  tibble(Date = dmy(d[[1]]), temp = as.numeric(d[[3]])) %>%
    filter(!is.na(Date), !is.na(temp)) %>%
    group_by(Date) %>%
    summarise(temp_max = max(temp), .groups = "drop") %>%
    filter(Date > min(Date), Date < max(Date)) %>%   # primer y último día incompletos
    mutate(Archivo = basename(path), Etiqueta = extraer_info_sala(path))
}

archivos <- list.files("Datos TP FINAL/No HOBO", pattern = "\\.txt$", full.names = TRUE)
datos_hobo <- map_dfr(archivos, leer_hobo)

# ------------------------------------------------------------------------------
# 3. Ajuste por serie
# ------------------------------------------------------------------------------
nombres_rezagos <- function(k) paste0("TMAX_l", 0:k)

ajustar_serie <- function(df_serie) {
  etiqueta <- df_serie$Etiqueta[1]

  # Con fechas repetidas (dos dataloggers en la misma serie) el AR(1) no tiene sentido
  if (anyDuplicated(df_serie$Date) > 0) {
    message("AVISO: '", etiqueta, "' tiene fechas repetidas (archivos: ",
            paste(unique(df_serie$Archivo), collapse = ", "), "). Se saltea.")
    return(NULL)
  }

  df <- df_serie %>%
    left_join(smn, by = "Date") %>%
    mutate(yday = yday(Date),
           sin_year = sin(2 * pi * yday / 365.25),
           cos_year = cos(2 * pi * yday / 365.25)) %>%
    drop_na(temp_max, all_of(nombres_rezagos(max_rezago)), TMAX_prom3d) %>%
    arrange(Date) %>%
    mutate(t = as.integer(Date - min(Date)))

  if (nrow(df) < 30) {
    message("AVISO: '", etiqueta, "' tiene solo ", nrow(df), " días completos. Se saltea.")
    return(NULL)
  }

  # Todos los modelos sobre las mismas filas
  formulas <- list(
    sin_exterior = temp_max ~ sin_year + cos_year,
    prom3d       = temp_max ~ TMAX_prom3d + sin_year + cos_year
  )
  for (k in 0:max_rezago) {
    formulas[[paste0("rezagos_0a", k)]] <-
      reformulate(c(nombres_rezagos(k), "sin_year", "cos_year"), response = "temp_max")
  }

  modelos <- imap(formulas, function(f, nombre) {
    tryCatch(
      gls(f, data = df, correlation = corAR1(form = ~ t), method = "ML"),
      error = function(e) {
        message("AVISO: '", etiqueta, "', modelo ", nombre, ": falló el ajuste (",
                conditionMessage(e), ")")
        NULL
      }
    )
  })
  if (any(map_lgl(modelos, is.null))) {
    message("AVISO: '", etiqueta, "' se saltea porque falló algún modelo.")
    return(NULL)
  }

  # --- Tabla de AIC de todos los modelos ---
  aic <- tibble(
    Etiqueta = etiqueta,
    Modelo   = names(modelos),
    AIC      = map_dbl(modelos, AIC),
    phi      = map_dbl(modelos, ~ coef(.x$modelStruct$corStruct, unconstrained = FALSE))
  )

  # --- Mejor cantidad de rezagos ---
  aic_rez  <- aic %>% filter(startsWith(Modelo, "rezagos"))
  mejor_nm <- aic_rez$Modelo[which.min(aic_rez$AIC)]
  k_mejor  <- as.integer(sub("rezagos_0a", "", mejor_nm))
  mejor    <- modelos[[mejor_nm]]

  # Coeficientes de cada rezago en el mejor modelo
  b  <- coef(mejor)
  V  <- vcov(mejor)
  nm <- nombres_rezagos(k_mejor)
  coefs <- tibble(
    Etiqueta = etiqueta,
    Rezago   = 0:k_mejor,
    Coef     = unname(b[nm]),
    EE       = unname(sqrt(diag(V)[nm]))
  )

  # Efecto acumulado = suma de los coeficientes; Var(suma) = 1' V 1
  efecto_acum <- sum(b[nm])
  ee_acum     <- sqrt(sum(V[nm, nm]))

  # ¿Pesos iguales? prom3d está anidado en rezagos_0a2 (b0 = b1 = b2)
  p_pesos_iguales <- anova(modelos$prom3d, modelos$rezagos_0a2)$`p-value`[2]

  resumen <- tibble(
    Etiqueta          = etiqueta,
    N_dias            = nrow(df),
    Rezagos_mejor     = paste0("0 a ", k_mejor),
    Efecto_acumulado  = efecto_acum,
    EE_acumulado      = ee_acum,
    Pendiente_prom3d  = unname(coef(modelos$prom3d)["TMAX_prom3d"]),
    # Negativo = los rezagos separados ajustan mejor que el promedio de 3 días
    dAIC_vs_prom3d    = AIC(mejor) - AIC(modelos$prom3d),
    # Negativo = el exterior aporta (comparado con solo el ciclo anual)
    dAIC_vs_sin_ext   = AIC(mejor) - AIC(modelos$sin_exterior),
    p_pesos_iguales   = p_pesos_iguales,
    # Colinealidad entre la TMAX de hoy y la de ayer
    cor_l0_l1         = cor(df$TMAX_l0, df$TMAX_l1)
  )

  list(resumen = resumen, coefs = coefs, aic = aic)
}

series     <- split(datos_hobo, datos_hobo$Etiqueta)
resultados <- map(series, ajustar_serie) %>% compact()

resumen_rezagos <- map_dfr(resultados, "resumen") %>% arrange(desc(Efecto_acumulado))
coefs_rezagos   <- map_dfr(resultados, "coefs")
aic_rezagos     <- map_dfr(resultados, "aic")

# ------------------------------------------------------------------------------
# 4. Resultados
# ------------------------------------------------------------------------------
cat("\n=== RESUMEN POR SERIE ===\n")
print(resumen_rezagos %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
      n = Inf, width = Inf)

cat("\n=== COEFICIENTES DE CADA REZAGO (mejor modelo) ===\n")
print(coefs_rezagos %>% mutate(across(c(Coef, EE), ~ round(.x, 3))), n = Inf)

# Coeficiente de cada rezago con IC del 95 %
g_rezagos <- ggplot(coefs_rezagos, aes(x = Rezago, y = Coef)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = Coef - 1.96 * EE, ymax = Coef + 1.96 * EE), color = "firebrick") +
  geom_line(color = "firebrick", alpha = 0.5) +
  facet_wrap(~ Etiqueta, ncol = 3) +
  scale_x_continuous(breaks = 0:max_rezago) +
  labs(
    title = "Efecto de la TMAX exterior de cada día anterior sobre la máxima interior",
    subtitle = "Mejor modelo por AIC (ML) en cada serie. Barras: IC 95 %",
    x = "Rezago (días)", y = "Coeficiente (°C interior por °C exterior)"
  ) +
  theme_minimal(base_size = 10)

print(g_rezagos)

# Efecto acumulado de los rezagos contra la pendiente del promedio de 3 días
g_comparacion <- ggplot(resumen_rezagos,
                        aes(x = Pendiente_prom3d, y = Efecto_acumulado, label = Etiqueta)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = Efecto_acumulado - 1.96 * EE_acumulado,
                    ymax = Efecto_acumulado + 1.96 * EE_acumulado), width = 0, color = "gray60") +
  geom_point(color = "steelblue", size = 2.5) +
  geom_text(size = 2.8, vjust = -0.8, check_overlap = TRUE) +
  labs(
    title = "Efecto acumulado de los rezagos vs. pendiente del promedio de 3 días",
    x = "Pendiente de TMAX_prom3d", y = "Suma de los coeficientes de los rezagos"
  ) +
  theme_minimal()

print(g_comparacion)
