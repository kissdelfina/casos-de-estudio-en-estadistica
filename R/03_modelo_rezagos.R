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
#     como prueba de robustez, también se informa el K con menor BIC, que
#     castiga cada parámetro con log(n) en vez de 2 y elige ventanas iguales
#     o más cortas;
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
# 1. SMN: TMAX diaria de la estación 87585 (Buenos Aires Observatorio)
# ------------------------------------------------------------------------------
# Leer el xlsx (72 MB) consume mucha memoria: la primera vez se guarda la
# estación filtrada en un .rds y las corridas siguientes leen ese archivo.
ruta_cache <- "datos/procesados/smn_87585_1991_2020.rds"
if (file.exists(ruta_cache)) {
  smn_1991 <- readRDS(ruta_cache)
} else {
  smn_1991 <- read_xlsx("datos/smn/smn_datos_meteorologicos_1991_2020.xlsx") %>%
    filter(NRO_OMM == 87585) %>%
    transmute(Date = as.Date(FECHA),
              TMAX = suppressWarnings(as.numeric(gsub(",", ".", TMAX))))  # "S/D" -> NA
  dir.create(dirname(ruta_cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(smn_1991, ruta_cache)
  gc()
}

smn_2021 <- read.delim("datos/smn/smn_datos_meteorologicos_desde_2021.lst") %>%
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
  # MACN-Or: son dos dataloggers con fechas superpuestas (colección antigua y
  # Azara), igual que en 01_modelo_temperatura.R y 05_modelo_pronostico.R
  if (grepl("azara", fname, ignore.case = TRUE) || grepl("04-01_2021", basename(filepath))) subsector <- "Azara"
  if (grepl("MACN-Or-03_02_2020", basename(filepath))) subsector <- "Antigua"

  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba|azara|antigua)", "", fname, ignore.case = TRUE)
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

archivos <- list.files("datos/hobo", pattern = "\\.txt$", full.names = TRUE)
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
    BIC      = map_dbl(modelos, BIC),
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

  # Mejor cantidad de rezagos según BIC, y su efecto acumulado
  k_bic   <- as.integer(sub("rezagos_0a", "", aic_rez$Modelo[which.min(aic_rez$BIC)]))
  m_bic   <- modelos[[paste0("rezagos_0a", k_bic)]]
  nm_bic  <- nombres_rezagos(k_bic)

  # ¿Hace falta el AR(1)? El mismo modelo sin correlación entre los errores
  sin_ar1 <- gls(formula(mejor), data = df, method = "ML")

  # ¿Pesos iguales? prom3d está anidado en rezagos_0a2 (b0 = b1 = b2)
  p_pesos_iguales <- anova(modelos$prom3d, modelos$rezagos_0a2)$`p-value`[2]

  resumen <- tibble(
    Etiqueta          = etiqueta,
    N_dias            = nrow(df),
    Rezagos_mejor     = paste0("0 a ", k_mejor),
    Efecto_acumulado  = efecto_acum,
    EE_acumulado      = ee_acum,
    Rezagos_BIC       = paste0("0 a ", k_bic),
    Efecto_acum_BIC   = sum(coef(m_bic)[nm_bic]),
    EE_acum_BIC       = sqrt(sum(vcov(m_bic)[nm_bic, nm_bic])),
    Pendiente_prom3d  = unname(coef(modelos$prom3d)["TMAX_prom3d"]),
    # Negativo = los rezagos separados ajustan mejor que el promedio de 3 días
    dAIC_vs_prom3d    = AIC(mejor) - AIC(modelos$prom3d),
    # Negativo = el exterior aporta (comparado con solo el ciclo anual)
    dAIC_vs_sin_ext   = AIC(mejor) - AIC(modelos$sin_exterior),
    p_pesos_iguales   = p_pesos_iguales,
    # Colinealidad entre la TMAX de hoy y la de ayer
    cor_l0_l1         = cor(df$TMAX_l0, df$TMAX_l1),
    # AR(1): phi estimado, cuánto empeora el AIC sin AR(1) (positivo = el AR(1)
    # ajusta mejor) y EE del efecto acumulado si se ignora la autocorrelación
    phi               = coef(mejor$modelStruct$corStruct, unconstrained = FALSE),
    dAIC_sin_AR1      = AIC(sin_ar1) - AIC(mejor),
    EE_acum_sin_AR1   = sqrt(sum(vcov(sin_ar1)[nm, nm]))
  )

  # Valores ajustados del mejor modelo (parte fija: rezagos + ciclo anual, sin
  # la corrección AR(1)), para graficar contra lo observado
  predicciones <- df %>%
    transmute(Etiqueta = etiqueta, Date, temp_max, TMAX = TMAX_l0,
              predicho = as.numeric(fitted(mejor)))

  list(resumen = resumen, coefs = coefs, aic = aic, predicciones = predicciones)
}

series     <- split(datos_hobo, datos_hobo$Etiqueta)
resultados <- map(series, ajustar_serie) %>% compact()

resumen_rezagos <- map_dfr(resultados, "resumen") %>% arrange(desc(Efecto_acumulado))
coefs_rezagos   <- map_dfr(resultados, "coefs")
aic_rezagos     <- map_dfr(resultados, "aic")
predicciones_rezagos <- map_dfr(resultados, "predicciones")

# ------------------------------------------------------------------------------
# 4. Resultados
# ------------------------------------------------------------------------------
cat("\n=== RESUMEN POR SERIE ===\n")
print(resumen_rezagos %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
      n = Inf, width = Inf)

cat("\n=== COEFICIENTES DE CADA REZAGO (mejor modelo) ===\n")
print(coefs_rezagos %>% mutate(across(c(Coef, EE), ~ round(.x, 3))), n = Inf)

dir.create("resultados/tablas", recursive = TRUE, showWarnings = FALSE)
write_csv(resumen_rezagos, "resultados/tablas/rezagos_resumen.csv")
write_csv(coefs_rezagos,   "resultados/tablas/rezagos_coeficientes.csv")
write_csv(aic_rezagos,     "resultados/tablas/rezagos_aic.csv")

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

# ------------------------------------------------------------------------------
# 5. Figura para el informe: MACN-Pv, observado contra ajustado
# ------------------------------------------------------------------------------
# Los tres dataloggers de Paleovertebrados se midieron en 2019, así que muestran
# el contraste entre sectores sin confundirlo con el año. Se agrega la TMAX
# exterior para ver cuánto de sus picos llega a cada sector.
etiquetas_pv <- c("MACN_Pv - Anexo"            = "Anexo (expuesto)",
                  "MACN_Pv - Ameghino"         = "Ameghino (aislado)",
                  "MACN_Pv - General / Normal" = "General")

datos_pv <- predicciones_rezagos %>%
  filter(Etiqueta %in% names(etiquetas_pv)) %>%
  # Completar los días sin datos con NA para que las líneas no crucen los huecos
  group_by(Etiqueta) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  ungroup() %>%
  left_join(resumen_rezagos %>% select(Etiqueta, Rezagos_mejor, Efecto_acumulado),
            by = "Etiqueta") %>%
  mutate(Panel = paste0(etiquetas_pv[Etiqueta], ": rezagos ", Rezagos_mejor,
                        " días, efecto acumulado ",
                        format(round(Efecto_acumulado, 2), decimal.mark = ",")),
         Panel = factor(Panel, levels = unique(Panel[order(match(Etiqueta, names(etiquetas_pv)))])))

g_pv <- ggplot(datos_pv, aes(x = Date)) +
  geom_line(aes(y = TMAX, color = "Máxima exterior (SMN)"), linewidth = 0.3, na.rm = TRUE) +
  geom_line(aes(y = temp_max, color = "Máxima interior observada"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = predicho, color = "Ajuste del modelo de rezagos"),
            linewidth = 0.5, linetype = "dashed", na.rm = TRUE) +
  scale_color_manual(values = c("Máxima exterior (SMN)"        = "gray75",
                                "Máxima interior observada"    = "black",
                                "Ajuste del modelo de rezagos" = "firebrick"),
                     breaks = c("Máxima exterior (SMN)", "Máxima interior observada",
                                "Ajuste del modelo de rezagos")) +
  facet_wrap(~ Panel, ncol = 1) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(title = "El Anexo sigue los picos del exterior; Ameghino los amortigua",
       subtitle = "Paleovertebrados (MACN-Pv), 2019: temperatura máxima diaria",
       x = NULL, y = "Temperatura máxima (°C)", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.grid.minor = element_blank())

print(g_pv)
dir.create("resultados/figuras", recursive = TRUE, showWarnings = FALSE)
ggsave("resultados/figuras/fig_modelo_pv.png", g_pv, width = 8, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 6. Figura para el informe: cómo se acumula el efecto a lo largo de los rezagos
# ------------------------------------------------------------------------------
# Suma acumulada de los coeficientes hasta cada rezago. Cada coeficiente por
# separado es inestable (rezagos muy correlacionados), pero la suma no, y su
# forma muestra la velocidad: las salas expuestas llegan a su efecto total en
# 1 o 2 días; las aisladas siguen sumando durante una a dos semanas.
nombre_corto <- function(e) {
  e <- sub(" - General / Normal", "", e)
  gsub("_", "-", sub(" - ", " ", e))
}
destacadas <- c("MACN_Pv - Anexo", "BA_v - General / Normal",
                "MACN_Pv - Ameghino", "MACN_Or - Antigua", "BA_Pb - General / Normal")

acumulado <- coefs_rezagos %>%
  arrange(Etiqueta, Rezago) %>%
  group_by(Etiqueta) %>%
  mutate(Acumulado = cumsum(Coef)) %>%
  ungroup() %>%
  mutate(Serie = nombre_corto(Etiqueta),
         Destacada = Etiqueta %in% destacadas)

etiquetas_fin <- acumulado %>%
  filter(Destacada) %>%
  group_by(Serie) %>%
  slice_max(Rezago, n = 1) %>%
  ungroup() %>%
  # El Anexo termina en el rezago 2, entre otras líneas: su etiqueta va arriba
  mutate(arriba = Etiqueta == "MACN_Pv - Anexo",
         x_txt  = if_else(arriba, Rezago - 0.3, Rezago + 0.25),
         y_txt  = if_else(arriba, Acumulado + 0.04, Acumulado))

g_acumulado <- ggplot(acumulado, aes(x = Rezago, y = Acumulado, group = Serie)) +
  geom_line(data = filter(acumulado, !Destacada), color = "gray80", linewidth = 0.5) +
  geom_line(data = filter(acumulado, Destacada), aes(color = Serie), linewidth = 1) +
  geom_point(data = filter(acumulado, Destacada), aes(color = Serie), size = 1.6) +
  geom_text(data = etiquetas_fin, aes(x = x_txt, y = y_txt, label = Serie, color = Serie),
            hjust = 0, size = 3.3) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  scale_x_continuous(breaks = 0:max_rezago, limits = c(0, max_rezago + 3.5)) +
  labs(title = "En las salas expuestas el efecto llega en pocos días;\nen las aisladas se sigue acumulando durante una a dos semanas",
       subtitle = "Suma de los coeficientes hasta cada rezago, en la ventana elegida por AIC. En gris, las demás series.",
       x = "Días hacia atrás (rezago)",
       y = "Efecto acumulado (°C interior por °C exterior)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12))

print(g_acumulado)
ggsave("resultados/figuras/fig_rezagos_acumulado.png", g_acumulado, width = 8, height = 5, dpi = 300)
