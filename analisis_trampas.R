# ==============================================================================
# SCRIPT CORREGIDO: GLS CON INERCIA TÉRMICA Y CICLO ESTACIONAL
# ==============================================================================

library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. PREPARACIÓN DE DATOS Y RETARDOS
# ------------------------------------------------------------------------------
path_macn_ma <- "~/Documents/UBA/Casos de Estudio en Estadística/Casos de Estudio en Estadística/Datos TP FINAL/MACN-Ma/MACN-Ma-03_02_2020-22_12_2020.txt"
targets_raw <- read.delim(path_macn_ma, sep = ";")

targets <- targets_raw %>%
  mutate(
    Date = as.Date(parse_date_time(Date, orders = c("dmy", "ymd", "mdy"))),
    temp_max = as.numeric(gsub(",", ".", Temp...C))
  ) %>%
  group_by(Date) %>%
  summarise(temp_max = max(temp_max, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Date))

smn_1991 <- read_xlsx('Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx')
smn_clean <- smn_1991 %>%
  filter(NRO_OMM == '87585') %>%
  mutate(
    Date = as.Date(FECHA),
    TMAX = as.numeric(gsub(",", ".", TMAX)),
    HUM_REL = as.numeric(gsub(",", ".", HUM_REL))
  ) %>%
  select(Date, TMAX, HUM_REL)

grilla_fechas <- data.frame(Date = seq(min(targets$Date), max(targets$Date), by = "day"))

df_modelo <- grilla_fechas %>%
  left_join(targets, by = "Date") %>%
  left_join(smn_clean, by = "Date") %>%
  arrange(Date) %>%
  mutate(
    TMAX_prom3d = (TMAX + dplyr::lag(TMAX, 1) + dplyr::lag(TMAX, 2)) / 3,
    
    # --------------------------------------------------------------------------
    # VARIABLE CLAVE: Componentes armónicos estacionales (Ciclo Anual)
    # --------------------------------------------------------------------------
    yday = yday(Date),
    sin_year = sin(2 * pi * yday / 365.25),
    cos_year = cos(2 * pi * yday / 365.25)
  ) %>%
  select(Date, temp_max, TMAX_prom3d, sin_year, cos_year) %>%
  na.omit() %>%
  mutate(t = as.integer(Date - min(Date)))

# ------------------------------------------------------------------------------
# 2. AJUSTE DEL MODELO GLS ESTACIONAL
# ------------------------------------------------------------------------------
m_estacional <- gls(
  temp_max ~ TMAX_prom3d + sin_year + cos_year,
  data = df_modelo,
  correlation = corAR1(form = ~ t)
)

summary(m_estacional)

# ------------------------------------------------------------------------------
# 3. PREDICCIONES Y EVALUACIÓN DE ESCALA
# ------------------------------------------------------------------------------
df_modelo$predicho <- predict(m_estacional, newdata = df_modelo)

r_val  <- cor(df_modelo$temp_max, df_modelo$predicho)
r2_val <- r_val^2

cat(sprintf("\nCorrelación Observado vs Predicho: %.4f\n", r_val))
cat(sprintf("R² del modelo: %.4f\n", r2_val))

# ------------------------------------------------------------------------------
# 4. GRÁFICOS DE CONTROL
# ------------------------------------------------------------------------------
# Serie temporal (Verificar que la curva roja acompañe la caída invernal)
ggplot(df_modelo, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Sensor (Observado)"), linewidth = 0.8) +
  geom_line(aes(y = predicho, color = "Modelo GLS Corregido"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Sensor (Observado)" = "black", "Modelo GLS Corregido" = "firebrick")) +
  labs(
    title = "MACN-Ma (2020): Ajuste Corregido (Estacionalidad + Inercia)",
    subtitle = sprintf("R² = %.3f | Corr = %.3f", r2_val, r_val),
    x = "Fecha", y = "Temperatura (°C)", color = "Serie"
  ) +
  theme_minimal()

# Scatter Plot (Verificar que los puntos sigan la diagonal 1:1)
ggplot(df_modelo, aes(x = predicho, y = temp_max)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Observado vs. Predicho (Escala Corregida)",
    x = "Predicho por Modelo GLS (°C)",
    y = "Observado por Sensor (°C)"
  ) +
  theme_minimal()