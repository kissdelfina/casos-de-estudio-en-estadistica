# Análisis exploratorio y diagnóstico predictivo de capturas de insectos en: MACN-Ba, MACN-En y MACN-Ma
library(tidyverse)
library(lubridate)
library(readxl)
library(patchwork)

archivos_capturas <- list.files("capturas", recursive = TRUE, full.names = TRUE, pattern = "[.]csv$")

df_capturas <- map_dfr(archivos_capturas, function(f) {
  partes <- strsplit(f, "/")[[1]]
  taxon  <- partes[2]
  sala   <- gsub("_.*", "", partes[3])
  
  read.csv(f) %>%
    mutate(
      fecha = as.Date(fecha),
      sala = sala,
      taxon = taxon
    )
})

resumen_observaciones <- df_capturas %>%
  group_by(sala, taxon) %>%
  summarise(
    n_fechas = n(),
    total_insectos = sum(cantidad),
    primera_fecha = min(fecha),
    ultima_fecha = max(fecha),
    .groups = "drop"
  )
print(as.data.frame(resumen_observaciones))


# cada cuanto midio
df_analisis <- df_capturas %>%
  arrange(sala, taxon, fecha) %>%
  group_by(sala, taxon) %>%
  mutate(
    fecha_ant = lag(fecha),
    dias_exp = as.numeric(fecha - fecha_ant),
    tasa_diaria = round(cantidad / dias_exp, 2)
  ) %>%
  ungroup()

cat("INTERVALOS DE MUESTREO")
cat("Mínimo:", min(df_analisis$dias_exp, na.rm=TRUE), "días\n")
cat("Mediana:", median(df_analisis$dias_exp, na.rm=TRUE), "días\n")
cat("Máximo:", max(df_analisis$dias_exp, na.rm=TRUE), "días (revisión tras el invierno en BA)\n")


# temperatura media 
smn_1991 <- read_xlsx("Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx") %>%
  filter(trimws(as.character(NRO_OMM)) == "87585") %>%
  mutate(
    Date = as.Date(FECHA),
    TMEDIA = as.numeric(gsub(",", ".", na_if(TMEDIA, "S/D"))),
    HUM_REL = as.numeric(gsub(",", ".", na_if(HUM_REL, "S/D")))
  ) %>%
  select(Date, TMEDIA, HUM_REL)

df_analisis <- df_analisis %>%
  rowwise() %>%
  mutate(
    # Si es la primera fecha o intervalo > 60 días, usamos los 14 días previos
    dias_ventana = if (is.na(dias_exp) || dias_exp > 60) 14 else dias_exp,
    temp_media_prev = mean(smn_1991$TMEDIA[smn_1991$Date <= fecha & smn_1991$Date >= (fecha - dias_ventana)], na.rm = TRUE),
    hum_media_prev  = mean(smn_1991$HUM_REL[smn_1991$Date <= fecha & smn_1991$Date >= (fecha - dias_ventana)], na.rm = TRUE)
  ) %>%
  ungroup()


# correlacion
cat("\n=== CORRELACIONES (CANTIDAD Y TASA DIARIA vs TEMPERATURA PREVIA) ===\n")
tabla_correlacion <- df_analisis %>%
  group_by(sala, taxon) %>%
  summarise(
    n = n(),
    r_cantidad_vs_temp = round(cor(cantidad, temp_media_prev), 3),
    r_tasa_vs_temp     = round(cor(tasa_diaria, temp_media_prev, use = "complete.obs"), 3),
    .groups = "drop"
  )
print(as.data.frame(tabla_correlacion))


# modelo de conteos con offset temporal
df_modelo <- df_analisis %>%
  filter(!is.na(dias_exp) & dias_exp < 100)

mod_pred <- glm(cantidad ~ temp_media_prev + taxon + sala + offset(log(dias_exp)), 
                family = quasipoisson, data = df_modelo)

cat("\n=== MODELO DE TASA DE CAPTURA (QUASIPOISSON CON OFFSET TEMPORAL) ===\n")
print(summary(mod_pred))

beta_temp <- coef(mod_pred)["temp_media_prev"]
cat(sprintf("\nPor cada 1°C de aumento en la temperatura media, la tasa de capturas aumenta un %.1f%% (exp(beta) = %.3f)\n", 
            (exp(beta_temp) - 1) * 100, exp(beta_temp)))


# graficos

# evolución temporal de capturas 
g1 <- ggplot(df_analisis, aes(x = fecha, y = cantidad, color = taxon)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  facet_wrap(~ sala, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Coleopteros" = "#D95F02", "Lepidopteros" = "#1B9E77")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(
    title = "Evolución Temporal de Capturas en las 3 Salas",
    subtitle = "Notar los saltos temporales irregulares entre revisiones de trampas",
    x = "Fecha de Revisión",
    y = "Cantidad de Insectos Capturados",
    color = "Taxón"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

print(g1)

# temperatura previa vs. captura diaria
g2 <- ggplot(df_analisis %>% filter(!is.na(tasa_diaria) & dias_exp < 100), 
             aes(x = temp_media_prev, y = tasa_diaria, color = taxon)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ sala, scales = "free_y") +
  scale_color_manual(values = c("Coleopteros" = "#D95F02", "Lepidopteros" = "#1B9E77")) +
  labs(
    title = "Temperatura Previa vs. Tasa Diaria de Capturas (Capturas / Días)",
    subtitle = "Corregir por días de exposición permite comparar períodos de distinta duración",
    x = "Temperatura Media del Intervalo de Exposición (°C)",
    y = "Tasa Diaria de Capturas (insectos/día)",
    color = "Taxón"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

print(g2)

# Guardar figuras para el informe
ggsave("fig_trampas_temporal.png", g1, width = 8, height = 5, dpi = 300)
ggsave("fig_trampas_temperatura.png", g2, width = 8, height = 4, dpi = 300)
