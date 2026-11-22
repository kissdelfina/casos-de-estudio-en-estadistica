# ==============================================================================
# RESIDUOS EN EL TIEMPO: SIN ESTRUCTURA, AR(1) Y ARMA(1,1)
# ==============================================================================
# Para una serie, ajusta el modelo de rezagos (con la ventana K que eligió
# 03_modelo_rezagos.R) con tres estructuras de error:
#
#   Sin estructura   gls sin correlación (equivale a mínimos cuadrados)
#   AR(1)            corAR1(form = ~ t)
#   ARMA(1,1)        corARMA(form = ~ t, p = 1, q = 1)
#
# y grafica, para cada una, los residuos normalizados contra la fecha y su ACF.
#
# Por qué residuos normalizados: los residuos crudos (observado - predicho) son
# casi iguales en los tres modelos, porque la estructura de error casi no cambia
# los coeficientes. Los normalizados descuentan la correlación que el modelo ya
# explica; si la estructura alcanza, deberían parecer ruido: sin rachas largas
# de un mismo signo y con la ACF dentro de las bandas. En el modelo sin
# estructura, el residuo normalizado es el crudo dividido por sigma.
#
# La serie se elige con la variable de entorno SERIE (por defecto, BA-Pb, que
# es continua y de un solo archivo). Ejemplo:
#   SERIE="MACN_Pv - Ameghino" Rscript R/08_residuos_tiempo.R
#
# Con SIN_HOY=1 se usa el modelo "sin hoy" del pronóstico (05): rezagos de 1 a
# K, sin la TMAX del mismo día. La ventana K se elige por AIC (ML, AR(1)) entre
# 1 y 14 con todos los datos de la serie. Ejemplo:
#   SERIE="MACN_Pv - Anexo" SIN_HOY=1 Rscript R/08_residuos_tiempo.R
#
# Necesita resultados/tablas/rezagos_resumen.csv, que genera el 03.
# ==============================================================================

library(lubridate)
library(tidyverse)
library(nlme)
library(patchwork)

serie_elegida <- Sys.getenv("SERIE", "BA_Pb - General / Normal")
sin_hoy <- Sys.getenv("SIN_HOY", "0") == "1"
desde <- if (sin_hoy) 1 else 0   # primer rezago del modelo
max_rezago <- 14

# ------------------------------------------------------------------------------
# 1. Datos (igual que 03_modelo_rezagos.R, usando el caché del SMN)
# ------------------------------------------------------------------------------
ruta_cache <- "datos/procesados/smn_87585_1991_2020.rds"
if (!file.exists(ruta_cache)) stop("Falta ", ruta_cache, ": correr antes R/03_modelo_rezagos.R")
smn_1991 <- readRDS(ruta_cache)

smn_2021 <- read.delim("datos/smn/smn_datos_meteorologicos_desde_2021.lst") %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
            TMAX = suppressWarnings(as.numeric(trimws(TMAX))))

smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date)
for (k in 0:max_rezago) smn[[paste0("TMAX_l", k)]] <- dplyr::lag(smn$TMAX, k)

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
  tibble(Date = dmy(d[[1]]), temp = as.numeric(d[[3]])) %>%
    filter(!is.na(Date), !is.na(temp)) %>%
    group_by(Date) %>%
    summarise(temp_max = max(temp), .groups = "drop") %>%
    filter(Date > min(Date), Date < max(Date)) %>%
    mutate(Etiqueta = extraer_info_sala(path))
}

archivos <- list.files("datos/hobo", pattern = "\\.txt$", full.names = TRUE)
archivos <- archivos[map_chr(archivos, extraer_info_sala) == serie_elegida]
if (length(archivos) == 0) stop("No hay archivos para la serie '", serie_elegida, "'")

k <- read_csv("resultados/tablas/rezagos_resumen.csv", show_col_types = FALSE) %>%
  filter(Etiqueta == serie_elegida) %>%
  pull(Rezagos_mejor) %>% sub("0 a ", "", .) %>% as.integer()

df <- map_dfr(archivos, leer_hobo) %>%
  left_join(smn, by = "Date") %>%
  mutate(sin_year = sin(2 * pi * yday(Date) / 365.25),
         cos_year = cos(2 * pi * yday(Date) / 365.25)) %>%
  drop_na(temp_max, all_of(paste0("TMAX_l", 0:max_rezago))) %>%
  arrange(Date) %>%
  mutate(t = as.integer(Date - min(Date)))

# ------------------------------------------------------------------------------
# 2. Los tres modelos
# ------------------------------------------------------------------------------
formula_k <- function(k) reformulate(c(paste0("TMAX_l", desde:k), "sin_year", "cos_year"),
                                     response = "temp_max")

# Sin hoy: la K del 03 es la del modelo con hoy, así que se vuelve a elegir
if (sin_hoy) {
  aics <- map_dbl(1:max_rezago, ~ AIC(gls(formula_k(.x), data = df, method = "ML",
                                          correlation = corAR1(form = ~ t))))
  k <- (1:max_rezago)[which.min(aics)]
}
cat("Serie:", serie_elegida, "| modelo:", if (sin_hoy) "sin hoy" else "con hoy",
    "| rezagos de", desde, "a", k, "\n")
f <- formula_k(k)
modelos <- list(
  "Sin estructura" = gls(f, data = df, method = "ML"),
  "AR(1)"          = gls(f, data = df, method = "ML", correlation = corAR1(form = ~ t)),
  "ARMA(1,1)"      = gls(f, data = df, method = "ML", correlation = corARMA(form = ~ t, p = 1, q = 1))
)

# Ljung-Box a 14 días sobre una grilla diaria (respeta los huecos)
grilla_res <- function(res, t) { g <- rep(NA_real_, max(t) + 1); g[t + 1] <- res; g }
ljung_box <- function(res, t, lag = 14, fitdf = 0) {
  r <- acf(grilla_res(res, t), lag.max = lag, na.action = na.pass, plot = FALSE)$acf[-1]
  n <- length(res)
  q <- n * (n + 2) * sum(r^2 / (n - seq_len(lag)))
  pchisq(q, df = lag - fitdf, lower.tail = FALSE)
}
params_cor <- c("Sin estructura" = 0, "AR(1)" = 1, "ARMA(1,1)" = 2)

residuos <- imap_dfr(modelos, ~ tibble(Modelo = .y, Date = df$Date,
                                        res = resid(.x, type = "normalized")))
acfs <- imap_dfr(modelos, function(m, nombre) {
  r <- acf(grilla_res(resid(m, type = "normalized"), df$t), lag.max = 21,
           na.action = na.pass, plot = FALSE)$acf[-1]
  tibble(Modelo = nombre, Rezago = 1:21, ACF = r)
})
resumen <- imap_dfr(modelos, function(m, nombre) {
  cor_txt <- if (is.null(m$modelStruct$corStruct)) "" else
    paste0(" · parámetros: ", paste(format(round(coef(m$modelStruct$corStruct, unconstrained = FALSE), 2),
                                           decimal.mark = ","), collapse = " / "))
  tibble(Modelo = nombre, AIC = AIC(m),
         p_LB = ljung_box(resid(m, type = "normalized"), df$t, 14, params_cor[[nombre]]),
         cor_txt = cor_txt)
}) %>%
  mutate(dAIC = AIC - min(AIC),
         Titulo = sprintf("%s — ΔAIC = %s · Ljung-Box (14 días) p %s%s", Modelo,
                          format(round(dAIC, 1), decimal.mark = ",", trim = TRUE),
                          ifelse(p_LB < 0.001, "< 0,001", paste("=", format(round(p_LB, 3), decimal.mark = ","))),
                          cor_txt))
# Desvío de los residuos crudos (°C): cuánto se equivoca el modelo dentro de la muestra
resumen$DE_residuo <- map_dbl(modelos, ~ sd(resid(.x, type = "response")))
print(resumen %>% select(Modelo, dAIC, p_LB, DE_residuo))

orden <- names(modelos)
etiquetas <- setNames(resumen$Titulo, resumen$Modelo)
residuos <- residuos %>% mutate(Modelo = factor(Modelo, levels = orden, labels = etiquetas[orden]))
acfs     <- acfs     %>% mutate(Modelo = factor(Modelo, levels = orden, labels = etiquetas[orden]))

# ------------------------------------------------------------------------------
# 3. Gráfico: residuos contra la fecha (izquierda) y su ACF (derecha)
# ------------------------------------------------------------------------------
tinta  <- "#2B3640"
gris   <- "#8A949C"
acento <- "#B8501F"
banda  <- 1.96 / sqrt(nrow(df))
nombre_serie <- serie_elegida %>% str_replace(" - General / Normal", "") %>%
  str_replace(" - ", " ") %>% str_replace("_", "-")

tema <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold", color = tinta),
        text = element_text(color = tinta),
        axis.text = element_text(color = gris))

g_tiempo <- ggplot(residuos, aes(Date, res)) +
  geom_hline(yintercept = 0, color = gris) +
  geom_hline(yintercept = c(-2, 2), color = gris, linetype = "dashed", linewidth = 0.4) +
  geom_line(color = tinta, linewidth = 0.4) +
  geom_point(color = tinta, size = 0.8) +
  # Suavizado: si queda una onda lenta, la estructura no alcanza
  geom_smooth(method = "loess", span = 0.15, se = FALSE, color = acento, linewidth = 1.1,
              formula = y ~ x) +
  facet_wrap(~ Modelo, ncol = 1) +
  # Meses en español sin depender del idioma del sistema
  scale_x_date(labels = function(x) paste(c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago",
                                            "sep", "oct", "nov", "dic")[month(x)], year(x))) +
  labs(x = NULL, y = "Residuo normalizado") +
  tema

g_acf <- ggplot(acfs, aes(Rezago, ACF)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -banda, ymax = banda,
           fill = gris, alpha = 0.18) +
  geom_hline(yintercept = 0, color = gris) +
  geom_segment(aes(xend = Rezago, yend = 0), color = tinta, linewidth = 1) +
  geom_point(aes(color = abs(ACF) > banda), size = 2) +
  scale_color_manual(values = c(`TRUE` = acento, `FALSE` = tinta), guide = "none") +
  facet_wrap(~ Modelo, ncol = 1, labeller = labeller(Modelo = function(x) rep(" ", length(x)))) +
  scale_x_continuous(breaks = c(1, 7, 14, 21)) +
  labs(x = "Rezago (días)", y = "ACF de los residuos") +
  tema

g <- (g_tiempo | g_acf) +
  plot_layout(widths = c(2.6, 1)) +
  plot_annotation(
    title = paste0(nombre_serie, ": residuos del modelo ",
                   if (sin_hoy) paste0("sin la TMAX de hoy (rezagos de 1 a ", k, ")")
                   else paste0("con la TMAX de hoy (rezagos de 0 a ", k, ")"),
                   ", con tres estructuras de error"),
    subtitle = paste0("Izquierda: residuos normalizados en el tiempo; línea naranja: suavizado (si ondula, queda correlación sin modelar). ",
                      "Líneas punteadas: ±2.\nDerecha: autocorrelación de esos residuos; franja gris: ±1,96/√n, ",
                      "lo esperable si fueran ruido; en naranja, los rezagos que se salen."),
    theme = theme(plot.title = element_text(face = "bold", size = 16, color = tinta),
                  plot.subtitle = element_text(size = 11, color = gris))
  )

nombre_archivo <- paste0("resultados/figuras/fig_residuos_tiempo_",
                         gsub("[^A-Za-z]+", "_", nombre_serie),
                         if (sin_hoy) "_sin_hoy" else "", ".png")
dir.create("resultados/figuras", recursive = TRUE, showWarnings = FALSE)
ggsave(nombre_archivo, g, width = 15, height = 9, dpi = 200, bg = "white")
cat("Figura guardada en", nombre_archivo, "\n")
