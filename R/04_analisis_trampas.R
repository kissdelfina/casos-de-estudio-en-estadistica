# ==============================================================================
# CAPTURAS DE INSECTOS EN TRAMPAS: MACN-Ba, MACN-En y MACN-Ma (2019-2020)
# ==============================================================================
# Cada fila es una revisión de trampa: cuántos insectos (coleópteros o
# lepidópteros) se juntaron desde la revisión anterior. Los intervalos entre
# revisiones son desiguales, así que se modela la cantidad con un offset
# log(días de exposición): el modelo describe la tasa diaria de capturas.
#
# Modelos (en este orden):
#   1. Poisson:  cantidad ~ temp_media_prev + taxon + sala + offset(log(dias_exp))
#      Se mide la dispersión con el estadístico de Pearson (chi^2 / gl). En un
#      Poisson debería dar cerca de 1.
#   2. Binomial negativa (MASS::glm.nb), con los mismos términos: agrega un
#      parámetro theta, con Var(Y) = mu + mu^2 / theta, para la sobredispersión.
#      Se compara con la Poisson por AIC y por razón de verosimilitud (Poisson
#      es el caso límite theta -> infinito, en el borde del espacio de
#      parámetros, así que el p-valor se divide por 2).
#   3. Cuasi-Poisson, solo como control: es lo que usaba la versión anterior
#      de este script.
#
# temp_media_prev: temperatura media exterior (SMN, estación 87585) promedio
# de los días de exposición. En la primera revisión de cada trampa y cuando el
# intervalo supera los 60 días (la revisión después del invierno) se usan los
# 14 días previos. Las revisiones con más de 100 días de exposición o sin
# revisión anterior no entran en el modelo.
# ==============================================================================

library(tidyverse)
library(lubridate)
library(readxl)
library(patchwork)

dir_tablas  <- "resultados/tablas"
dir_figuras <- "resultados/figuras"
dir.create(dir_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figuras, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Capturas
# ------------------------------------------------------------------------------
# datos/trampas/<taxon>/<sala>_<taxon>.csv, con columnas fecha y cantidad
archivos_capturas <- list.files("datos/trampas", recursive = TRUE, full.names = TRUE,
                                pattern = "_(coleopteros|lepidopteros)[.]csv$")

df_capturas <- map_dfr(archivos_capturas, function(f) {
  read.csv(f) %>%
    mutate(fecha = as.Date(fecha),
           sala  = gsub("_.*", "", basename(f)),
           taxon = str_to_title(basename(dirname(f))))
})

resumen_observaciones <- df_capturas %>%
  group_by(sala, taxon) %>%
  summarise(n_fechas = n(), total_insectos = sum(cantidad),
            primera_fecha = min(fecha), ultima_fecha = max(fecha), .groups = "drop")
print(as.data.frame(resumen_observaciones))

df_analisis <- df_capturas %>%
  arrange(sala, taxon, fecha) %>%
  group_by(sala, taxon) %>%
  mutate(fecha_ant   = lag(fecha),
         dias_exp    = as.numeric(fecha - fecha_ant),
         tasa_diaria = round(cantidad / dias_exp, 2)) %>%
  ungroup()

cat("\n=== INTERVALOS DE MUESTREO ===\n")
cat("Mínimo:", min(df_analisis$dias_exp, na.rm = TRUE), "días\n")
cat("Mediana:", median(df_analisis$dias_exp, na.rm = TRUE), "días\n")
cat("Máximo:", max(df_analisis$dias_exp, na.rm = TRUE), "días (revisión tras el invierno)\n")

# ------------------------------------------------------------------------------
# 2. SMN: temperatura media y humedad diarias
# ------------------------------------------------------------------------------
# Caché propia (la de 01/03/05 solo guarda TMAX). La primera vez lee el xlsx de
# 72 MB, que tarda unos minutos.
ruta_cache <- "datos/procesados/smn_87585_tmedia_hum_1991_2020.rds"
if (file.exists(ruta_cache)) {
  smn <- readRDS(ruta_cache)
} else {
  smn <- read_xlsx("datos/smn/smn_datos_meteorologicos_1991_2020.xlsx") %>%
    filter(NRO_OMM == 87585) %>%
    transmute(Date    = as.Date(FECHA),
              TMEDIA  = suppressWarnings(as.numeric(gsub(",", ".", TMEDIA))),   # "S/D" -> NA
              HUM_REL = suppressWarnings(as.numeric(gsub(",", ".", HUM_REL))))
  dir.create(dirname(ruta_cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(smn, ruta_cache)
  gc()
}

promedio_previo <- function(variable, fecha, dias) {
  mean(smn[[variable]][smn$Date <= fecha & smn$Date >= fecha - dias], na.rm = TRUE)
}

df_analisis <- df_analisis %>%
  mutate(dias_ventana    = if_else(is.na(dias_exp) | dias_exp > 60, 14, dias_exp),
         temp_media_prev = pmap_dbl(list(fecha, dias_ventana), ~ promedio_previo("TMEDIA", ..1, ..2)),
         hum_media_prev  = pmap_dbl(list(fecha, dias_ventana), ~ promedio_previo("HUM_REL", ..1, ..2)))

cat("\n=== CORRELACIONES (CANTIDAD Y TASA DIARIA vs TEMPERATURA PREVIA) ===\n")
df_analisis %>%
  group_by(sala, taxon) %>%
  summarise(n = n(),
            r_cantidad_vs_temp = round(cor(cantidad, temp_media_prev), 3),
            r_tasa_vs_temp     = round(cor(tasa_diaria, temp_media_prev, use = "complete.obs"), 3),
            .groups = "drop") %>%
  as.data.frame() %>% print()

# ------------------------------------------------------------------------------
# 3. Modelos de conteo con offset
# ------------------------------------------------------------------------------
df_modelo <- df_analisis %>% filter(!is.na(dias_exp), dias_exp < 100)
cat("\nRevisiones en el modelo:", nrow(df_modelo), "\n")

formula_trampas <- cantidad ~ temp_media_prev + taxon + sala + offset(log(dias_exp))

# 3.1 Poisson y su dispersión
mod_pois <- glm(formula_trampas, family = poisson, data = df_modelo)
chi2_pearson <- sum(residuals(mod_pois, type = "pearson")^2)
gl_pois      <- df.residual(mod_pois)
dispersion   <- chi2_pearson / gl_pois
p_dispersion <- pchisq(chi2_pearson, gl_pois, lower.tail = FALSE)

cat("\n=== POISSON ===\n")
print(summary(mod_pois))
cat(sprintf("Dispersión (chi^2 de Pearson / gl) = %.1f / %d = %.2f, p = %.2g\n",
            chi2_pearson, gl_pois, dispersion, p_dispersion))

# 3.2 Binomial negativa
mod_nb <- MASS::glm.nb(formula_trampas, data = df_modelo)
cat("\n=== BINOMIAL NEGATIVA ===\n")
print(summary(mod_nb))
dispersion_nb <- sum(residuals(mod_nb, type = "pearson")^2) / df.residual(mod_nb)

# Razón de verosimilitud Poisson vs binomial negativa (theta en el borde: p / 2)
lrt_nb <- as.numeric(2 * (logLik(mod_nb) - logLik(mod_pois)))
p_lrt  <- pchisq(lrt_nb, df = 1, lower.tail = FALSE) / 2

# 3.3 Cuasi-Poisson (control)
mod_qp <- glm(formula_trampas, family = quasipoisson, data = df_modelo)

# 3.4 Sensibilidad: ¿agrega algo la humedad exterior en la binomial negativa?
mod_nb_hum <- MASS::glm.nb(update(formula_trampas, . ~ . + hum_media_prev), data = df_modelo)

# ------------------------------------------------------------------------------
# 4. Tablas
# ------------------------------------------------------------------------------
efecto_temp <- function(modelo, nombre) {
  b  <- coef(modelo)["temp_media_prev"]
  ee <- sqrt(vcov(modelo)["temp_media_prev", "temp_media_prev"])
  tibble(Modelo = nombre, beta = b, EE = ee,
         p = 2 * pnorm(-abs(b / ee)),
         Cambio_pct    = 100 * (exp(b) - 1),
         IC95_inf_pct  = 100 * (exp(b - 1.96 * ee) - 1),
         IC95_sup_pct  = 100 * (exp(b + 1.96 * ee) - 1))
}

comparacion <- tibble(
  Modelo     = c("Poisson", "Binomial negativa"),
  AIC        = c(AIC(mod_pois), AIC(mod_nb)),
  LogLik     = c(as.numeric(logLik(mod_pois)), as.numeric(logLik(mod_nb))),
  Dispersion = c(dispersion, dispersion_nb),
  theta      = c(NA, mod_nb$theta),
  EE_theta   = c(NA, mod_nb$SE.theta)
)
efectos <- bind_rows(efecto_temp(mod_pois, "Poisson"),
                     efecto_temp(mod_nb, "Binomial negativa"),
                     efecto_temp(mod_qp, "Cuasi-Poisson")) %>%
  # en la cuasi-Poisson el p-valor es con t; se recalcula con summary()
  mutate(p = if_else(Modelo == "Cuasi-Poisson",
                     summary(mod_qp)$coefficients["temp_media_prev", 4], p))

coefs_nb <- tibble(
  Termino   = names(coef(mod_nb)),
  beta      = unname(coef(mod_nb)),
  EE        = unname(sqrt(diag(vcov(mod_nb)))),
  exp_beta  = exp(beta),
  p         = 2 * pnorm(-abs(beta / EE))
)

cat("\n=== COMPARACIÓN POISSON vs BINOMIAL NEGATIVA ===\n")
print(comparacion)
cat(sprintf("Razón de verosimilitud = %.1f, p = %.2g (theta en el borde, p/2)\n", lrt_nb, p_lrt))
cat(sprintf("Dispersión de Poisson: %.2f (p = %.2g)\n", dispersion, p_dispersion))
cat("\n=== EFECTO DE LA TEMPERATURA (% de cambio en la tasa por °C) ===\n")
print(efectos %>% mutate(across(where(is.double), ~ signif(.x, 3))))
cat("\n=== COEFICIENTES DE LA BINOMIAL NEGATIVA ===\n")
print(coefs_nb %>% mutate(across(where(is.double), ~ signif(.x, 3))))
cat(sprintf("\nHumedad en la binomial negativa: Delta AIC = %.1f, p = %.2g\n",
            AIC(mod_nb_hum) - AIC(mod_nb),
            summary(mod_nb_hum)$coefficients["hum_media_prev", 4]))

write_csv(comparacion, file.path(dir_tablas, "trampas_comparacion_modelos.csv"))
write_csv(efectos,     file.path(dir_tablas, "trampas_efecto_temperatura.csv"))
write_csv(coefs_nb,    file.path(dir_tablas, "trampas_coeficientes_nb.csv"))

# ------------------------------------------------------------------------------
# 5. Gráficos
# ------------------------------------------------------------------------------
colores_taxon <- c("Coleopteros" = "#D95F02", "Lepidopteros" = "#1B9E77")

# Evolución temporal de capturas
g1 <- ggplot(df_analisis, aes(x = fecha, y = cantidad, color = taxon)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  facet_wrap(~ sala, ncol = 1, scales = "free_y") +
  scale_color_manual(values = colores_taxon) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(title = "Evolución Temporal de Capturas en las 3 Salas",
       subtitle = "Notar los saltos temporales irregulares entre revisiones de trampas",
       x = "Fecha de Revisión", y = "Cantidad de Insectos Capturados", color = "Taxón") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
print(g1)

# Temperatura previa vs. captura diaria
g2 <- ggplot(df_analisis %>% filter(!is.na(tasa_diaria), dias_exp < 100),
             aes(x = temp_media_prev, y = tasa_diaria, color = taxon)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ sala, scales = "free_y") +
  scale_color_manual(values = colores_taxon) +
  labs(title = "Temperatura Previa vs. Tasa Diaria de Capturas (Capturas / Días)",
       subtitle = "Corregir por días de exposición permite comparar períodos de distinta duración",
       x = "Temperatura Media del Intervalo de Exposición (°C)",
       y = "Tasa Diaria de Capturas (insectos/día)", color = "Taxón") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
print(g2)

ggsave(file.path(dir_figuras, "fig_trampas_temporal.png"), g1, width = 8, height = 5, dpi = 300)
ggsave(file.path(dir_figuras, "fig_trampas_temperatura.png"), g2, width = 8, height = 4, dpi = 300)
