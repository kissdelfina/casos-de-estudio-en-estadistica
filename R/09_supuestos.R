# ==============================================================================
# SUPUESTOS DEL MODELO: NORMALIDAD Y VARIANZA CONSTANTE DE LOS RESIDUOS
# ==============================================================================
# Para cada serie ajusta el modelo principal (rezagos 0..K de la TMAX, seno y
# coseno anuales, errores AR(1)) con la ventana K que eligió 03_modelo_rezagos.R
# y revisa dos supuestos sobre los residuos normalizados:
#
#   Normalidad          QQ-plot. Medida: JB/n = asimetría²/6 + curtosis_exceso²/24
#                       (el estadístico de Jarque-Bera dividido por n, para que
#                       no dependa del tamaño de la serie). 0 = normal.
#   Varianza constante  Residuos contra valores ajustados. Medida: cociente entre
#                       el desvío de los residuos en el tercio de ajustados más
#                       alto y en el más bajo. 1 = varianza constante.
#
# Los residuos normalizados ya descuentan la correlación AR(1): si el modelo
# está bien, son aproximadamente N(0, 1) independientes.
#
# Grafica las 3 series que mejor cumplen y las 3 que peor cumplen cada
# supuesto, y además una versión con todas las series (respaldo).
#
# Necesita resultados/tablas/rezagos_resumen.csv (lo genera el 03) y el caché
# datos/procesados/smn_87585_1991_2020.rds.
#
# Con la variable de entorno VARIABLE=min revisa el modelo de la mínima
# interior contra la TMIN, con la ventana K de rezagos_resumen_minymax.csv y
# el caché smn_87585_tmin_1991_2020.rds (los genera modelo_rezagos_minymax.R).
# Las tablas y figuras llevan el sufijo _min:
#   VARIABLE=min Rscript R/09_supuestos.R
# ==============================================================================

library(lubridate)
library(tidyverse)
library(nlme)
library(patchwork)

max_rezago <- 14

variable <- Sys.getenv("VARIABLE", "max")
if (!variable %in% c("max", "min")) stop("VARIABLE tiene que ser max o min")
exterior <- if (variable == "max") "TMAX" else "TMIN"
sufijo   <- if (variable == "max") "" else "_min"

# ------------------------------------------------------------------------------
# 1. Datos (igual que 03_modelo_rezagos.R, usando el caché del SMN)
# ------------------------------------------------------------------------------
ruta_cache <- "datos/procesados/smn_87585_1991_2020.rds"
if (!file.exists(ruta_cache)) stop("Falta ", ruta_cache, ": correr antes R/03_modelo_rezagos.R")
smn_1991 <- readRDS(ruta_cache)
if (variable == "min") {
  ruta_tmin <- "datos/procesados/smn_87585_tmin_1991_2020.rds"
  if (!file.exists(ruta_tmin)) stop("Falta ", ruta_tmin, ": correr antes R/modelo_rezagos_minymax.R")
  smn_1991 <- left_join(smn_1991, readRDS(ruta_tmin), by = "Date")
}

smn_2021 <- read.delim("datos/smn/smn_datos_meteorologicos_desde_2021.lst") %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
            TMAX = suppressWarnings(as.numeric(trimws(TMAX))),
            TMIN = suppressWarnings(as.numeric(trimws(TMIN))))

smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date)
# EXT es la TMAX o la TMIN, según la variable que se revisa
smn$EXT <- smn[[exterior]]
for (k in 0:max_rezago) smn[[paste0("EXT_l", k)]] <- dplyr::lag(smn$EXT, k)

extraer_info_sala <- function(filepath) {
  fname <- basename(filepath)
  original <- fname
  fname <- gsub("[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}", "", fname)
  fname <- gsub("\\.txt$", "", fname)
  subsector <- "General / Normal"
  if (grepl("Anexo", fname, ignore.case = TRUE)) subsector <- "Anexo"
  if (grepl("Ameghino", fname, ignore.case = TRUE)) subsector <- "Ameghino"
  if (grepl("coleccion", fname, ignore.case = TRUE)) subsector <- "Colección"
  if (grepl("laboratorio", fname, ignore.case = TRUE)) subsector <- "Laboratorio"
  if (grepl("abajo", fname, ignore.case = TRUE)) subsector <- "Abajo"
  if (grepl("arriba", fname, ignore.case = TRUE)) subsector <- "Arriba"
  # MACN-Or: dos dataloggers (igual que en el 03)
  if (grepl("azara", fname, ignore.case = TRUE) || grepl("04-01_2021", original)) subsector <- "Azara"
  if (grepl("MACN-Or-03_02_2020", original)) subsector <- "Antigua"
  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba|azara|antigua)", "", fname, ignore.case = TRUE)
  sala_base <- gsub("-", "_", sala_base)
  paste(sala_base, subsector, sep = " - ")
}

leer_hobo <- function(path) {
  d <- read.csv2(path, fileEncoding = "UTF-8-BOM")
  tibble(Date = dmy(d[[1]]), temp = as.numeric(d[[3]])) %>%
    filter(!is.na(Date), !is.na(temp)) %>%
    group_by(Date) %>%
    summarise(temp = if (variable == "max") max(temp) else min(temp), .groups = "drop") %>%
    filter(Date > min(Date), Date < max(Date)) %>%
    mutate(Etiqueta = extraer_info_sala(path))
}

archivos <- list.files("datos/hobo", pattern = "\\.txt$", full.names = TRUE)
datos_hobo <- map_dfr(archivos, leer_hobo)

ventanas <- if (variable == "max") {
  read_csv("resultados/tablas/rezagos_resumen.csv", show_col_types = FALSE)
} else {
  read_csv("resultados/tablas/rezagos_resumen_minymax.csv", show_col_types = FALSE) %>%
    filter(Respuesta == "temp_min")
}
ventanas <- ventanas %>%
  transmute(Etiqueta, K = as.integer(sub("0 a ", "", Rezagos_mejor)))

nombre_corto <- function(x) x %>%
  str_replace("MACN_Pv - General / Normal", "MACN_Pv - General") %>%
  str_replace(" - General / Normal", "") %>%
  str_replace(" - ", " ") %>% str_replace("_", "-")

# ------------------------------------------------------------------------------
# 2. Ajuste y residuos por serie
# ------------------------------------------------------------------------------
ajustar_serie <- function(etiqueta, k) {
  df <- datos_hobo %>%
    filter(Etiqueta == etiqueta) %>%
    left_join(smn, by = "Date") %>%
    mutate(sin_year = sin(2 * pi * yday(Date) / 365.25),
           cos_year = cos(2 * pi * yday(Date) / 365.25)) %>%
    drop_na(temp, all_of(paste0("EXT_l", 0:max_rezago))) %>%
    arrange(Date) %>%
    mutate(t = as.integer(Date - min(Date)))
  f <- reformulate(c(paste0("EXT_l", 0:k), "sin_year", "cos_year"), response = "temp")
  m <- tryCatch(gls(f, data = df, method = "ML", correlation = corAR1(form = ~ t)),
                error = function(e) { message("AVISO: '", etiqueta, "': ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  tibble(Etiqueta = etiqueta, Date = df$Date,
         ajustado = as.numeric(fitted(m)),
         res = as.numeric(resid(m, type = "normalized")))
}

residuos <- map2_dfr(ventanas$Etiqueta, ventanas$K, ajustar_serie) %>%
  mutate(Serie = nombre_corto(Etiqueta))

# ------------------------------------------------------------------------------
# 3. Medidas de cada supuesto
# ------------------------------------------------------------------------------
medidas <- residuos %>%
  group_by(Serie) %>%
  mutate(tercio = ntile(ajustado, 3)) %>%
  summarise(
    n          = n(),
    asimetria  = mean((res - mean(res))^3) / sd(res)^3,
    curtosis   = mean((res - mean(res))^4) / sd(res)^4 - 3,
    JB_n       = asimetria^2 / 6 + curtosis^2 / 24,
    p_JB       = pchisq(n * JB_n, df = 2, lower.tail = FALSE),
    pct_fuera3 = 100 * mean(abs(res) > 3),     # en una normal: 0,27 %
    DE_bajo    = sd(res[tercio == 1]),
    DE_alto    = sd(res[tercio == 3]),
    cociente_DE = DE_alto / DE_bajo,
    .groups = "drop"
  ) %>%
  arrange(JB_n)

cat("\n=== SUPUESTOS DE LOS RESIDUOS NORMALIZADOS ===\n")
print(medidas %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
write_csv(medidas, paste0("resultados/tablas/supuestos_residuos", sufijo, ".csv"))

# Días con residuos extremos (|residuo| > 4): candidatos a errores de medición
# o a eventos puntuales (aire acondicionado, calefacción, puerta abierta, etc.)
extremos <- residuos %>% filter(abs(res) > 4) %>% arrange(Serie, Date) %>%
  select(Serie, Date, ajustado, res)
cat("\n=== RESIDUOS EXTREMOS (|residuo| > 4) ===\n")
print(extremos, n = Inf)
write_csv(extremos, paste0("resultados/tablas/supuestos_extremos", sufijo, ".csv"))

# ------------------------------------------------------------------------------
# 4. Gráficos
# ------------------------------------------------------------------------------
tinta  <- "#2B3640"
gris   <- "#8A949C"
acento <- "#B8501F"
coma   <- function(x, d = 2) format(round(x, d), decimal.mark = ",", nsmall = d, trim = TRUE)

tema <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold", color = tinta),
        text = element_text(color = tinta),
        axis.text = element_text(color = gris),
        plot.subtitle = element_text(color = gris))

# --- Normalidad: QQ-plot ---
etiq_qq <- medidas %>%
  transmute(Serie, Panel = sprintf("%s (n = %d)\nasimetría %s · curtosis exceso %s",
                                   Serie, n, coma(asimetria), coma(curtosis)))

qq_plot <- function(series, subtitulo) {
  d <- residuos %>% filter(Serie %in% series) %>%
    left_join(etiq_qq, by = "Serie") %>%
    mutate(Panel = factor(Panel, levels = etiq_qq$Panel[match(series, etiq_qq$Serie)]))
  ggplot(d, aes(sample = res)) +
    geom_abline(slope = 1, intercept = 0, color = acento, linewidth = 0.8) +
    stat_qq(color = tinta, size = 0.9, alpha = 0.6) +
    facet_wrap(~ Panel, nrow = 1) +
    coord_cartesian(xlim = c(-4, 4), ylim = range(c(-5, 5, residuos$res))) +
    labs(subtitle = subtitulo, x = "Cuantil teórico N(0, 1)", y = "Residuo normalizado") +
    tema
}

mejores_qq <- head(medidas$Serie, 3)
peores_qq  <- rev(tail(medidas$Serie, 3))
g_qq <- qq_plot(mejores_qq, "Las 3 series más cercanas a la normal") /
  qq_plot(peores_qq, "Las 3 más alejadas") +
  plot_annotation(
    title = paste0("Normalidad de los residuos", if (variable == "min") " (mínima)" else "",
                   ": colas más pesadas que la normal en las peores series"),
    subtitle = "Residuos normalizados del modelo con rezagos y AR(1). Si fueran normales, los puntos seguirían la recta.",
    theme = theme(plot.title = element_text(color = tinta, size = 17),
                  plot.subtitle = element_text(color = gris, size = 13)))
ggsave(paste0("resultados/figuras/fig_supuestos_qq", sufijo, ".png"), g_qq, width = 13, height = 9, dpi = 300)

# --- Varianza constante: residuos contra ajustados ---
orden_var <- medidas %>% arrange(abs(log(cociente_DE))) %>% pull(Serie)
etiq_var <- medidas %>%
  transmute(Serie, Panel = sprintf("%s\nDE tercio alto / tercio bajo = %s", Serie, coma(cociente_DE)))

var_plot <- function(series, subtitulo) {
  d <- residuos %>% filter(Serie %in% series) %>%
    left_join(etiq_var, by = "Serie") %>%
    mutate(Panel = factor(Panel, levels = etiq_var$Panel[match(series, etiq_var$Serie)]))
  ggplot(d, aes(ajustado, res)) +
    geom_hline(yintercept = 0, color = gris) +
    geom_hline(yintercept = c(-2, 2), color = gris, linetype = "dashed", linewidth = 0.4) +
    geom_point(color = tinta, size = 0.9, alpha = 0.5) +
    # Suavizado de |residuo|: si sube o baja, la varianza no es constante
    geom_smooth(aes(y = abs(res)), method = "loess", formula = y ~ x, se = FALSE,
                color = acento, linewidth = 1.1) +
    facet_wrap(~ Panel, nrow = 1, scales = "free_x") +
    labs(subtitle = subtitulo, x = "Valor ajustado (°C)", y = "Residuo normalizado") +
    tema
}

mejores_var <- head(orden_var, 3)
peores_var  <- rev(tail(orden_var, 3))
g_var <- var_plot(mejores_var, "Las 3 series con la varianza más pareja") /
  var_plot(peores_var, "Las 3 con la varianza menos pareja") +
  plot_annotation(
    title = paste0("Varianza de los residuos según el valor ajustado", if (variable == "min") " (mínima)" else ""),
    subtitle = "Línea naranja: promedio de |residuo|. Si la varianza es constante, queda plana (cerca de 0,8).",
    theme = theme(plot.title = element_text(color = tinta, size = 17),
                  plot.subtitle = element_text(color = gris, size = 13)))
ggsave(paste0("resultados/figuras/fig_supuestos_varianza", sufijo, ".png"), g_var, width = 13, height = 9, dpi = 300)

# --- Respaldo: todas las series ---
g_qq_todas <- residuos %>%
  left_join(etiq_qq, by = "Serie") %>%
  mutate(Panel = factor(Panel, levels = etiq_qq$Panel)) %>%
  ggplot(aes(sample = res)) +
  geom_abline(slope = 1, intercept = 0, color = acento, linewidth = 0.8) +
  stat_qq(color = tinta, size = 0.7, alpha = 0.6) +
  facet_wrap(~ Panel, ncol = 4) +
  labs(title = "QQ-plot de los residuos normalizados, de la más normal a la menos",
       x = "Cuantil teórico N(0, 1)", y = "Residuo normalizado") +
  tema
ggsave(paste0("resultados/figuras/fig_supuestos_qq_todas", sufijo, ".png"), g_qq_todas, width = 15, height = 11, dpi = 300)

g_var_todas <- residuos %>%
  left_join(etiq_var, by = "Serie") %>%
  mutate(Panel = factor(Panel, levels = etiq_var$Panel[match(orden_var, etiq_var$Serie)])) %>%
  ggplot(aes(ajustado, res)) +
  geom_hline(yintercept = 0, color = gris) +
  geom_point(color = tinta, size = 0.7, alpha = 0.5) +
  geom_smooth(aes(y = abs(res)), method = "loess", formula = y ~ x, se = FALSE,
              color = acento, linewidth = 1) +
  facet_wrap(~ Panel, ncol = 4, scales = "free_x") +
  labs(title = "Residuos contra ajustados, de la varianza más pareja a la menos",
       x = "Valor ajustado (°C)", y = "Residuo normalizado") +
  tema
ggsave(paste0("resultados/figuras/fig_supuestos_varianza_todas", sufijo, ".png"), g_var_todas, width = 15, height = 11, dpi = 300)
