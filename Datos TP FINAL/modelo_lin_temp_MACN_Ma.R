library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)
library(ggplot2)

# 1. CARGA Y PREPROCESAMIENTO
path_macn_ma <- "~/Documents/UBA/Casos de Estudio en Estadística/Casos de Estudio en Estadística/Datos TP FINAL/MACN-Ma/MACN-Ma-03_02_2020-22_12_2020.txt"
targets_raw <- read.delim(path_macn_ma, sep = ";")

targets <- targets_raw %>%
  mutate(
    Date = as.Date(parse_date_time(Date, orders = c("dmy", "ymd", "mdy"))),
    temp_max = as.numeric(gsub(",", ".", Temp...C))
  ) %>%
  group_by(Date) %>%
  summarise(temp_max = max(temp_max, na.rm = TRUE), .groups = "drop")

smn_1991 <- read_xlsx('Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx')
smn_1991 <- smn_1991[smn_1991$NRO_OMM == '87585', ]
smn_1991$Date <- as.Date(smn_1991$FECHA)

# 2. CREACIÓN DE GRILLA TEMPORAL COMPLETA Y REZAGOS (LAGS)
# Se crea un calendario diario continuo para asegurar que lag(1) sea exactamente el día anterior
df_grid <- data.frame(Date = seq(min(targets$Date), max(targets$Date), by = "day"))

df_completo <- df_grid %>%
  left_join(targets, by = "Date") %>%
  left_join(smn_1991, by = "Date") %>%
  arrange(Date) %>%
  mutate(
    TMAX = as.numeric(gsub(",", ".", TMAX)),
    HUM_REL = as.numeric(gsub(",", ".", HUM_REL)),
    
    # Rezagos de TMAX (1 día atrás y 2 días atrás)
    TMAX_lag1 = dplyr::lag(TMAX, 1),
    TMAX_lag2 = dplyr::lag(TMAX, 2),
    
    # Promedio acumulado de TMAX de los últimos 2 días (inercia térmica)
    TMAX_mean2d = (TMAX + dplyr::lag(TMAX, 1)) / 2,
    
    # Índice temporal para estructura corAR1
    t = as.integer(Date - min(Date))
  )

# Filtrar filas sin datos (los primeros días pierden valor por el lag)
df_model <- na.omit(df_completo[, c("Date", "t", "temp_max", "TMAX", "TMAX_lag1", "TMAX_lag2", "TMAX_mean2d")])

# 3. COMPARACIÓN DE MODELOS GLS CON DIFERENTES ESTRUCTURAS DE DESFASE
m_mismo_dia <- gls(temp_max ~ TMAX, 
                   data = df_model, correlation = corAR1(form = ~ t))

m_lag1      <- gls(temp_max ~ TMAX_lag1, 
                   data = df_model, correlation = corAR1(form = ~ t))

m_combinado <- gls(temp_max ~ TMAX + TMAX_lag1, 
                   data = df_model, correlation = corAR1(form = ~ t))

m_prom2d    <- gls(temp_max ~ TMAX_mean2d, 
                   data = df_model, correlation = corAR1(form = ~ t))

# Comparar AIC (Menor AIC indica mejor ajuste)
model_comp <- AIC(m_mismo_dia, m_lag1, m_combinado, m_prom2d)
print("Comparación de Criterio de Akaike (AIC):")
print(model_comp)

# Select el modelo con menor AIC para predecir (ej. m_prom2d o m_lag1)
modelo_optimo <- m_prom2d  # Cambiar según el resultado del AIC
df_model$predichos <- predict(modelo_optimo, newdata = df_model)

# 4. GRÁFICO DE AJUSTE CON DESFASE
ggplot(df_model, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Observado (Sensor)"), linewidth = 0.8) +
  geom_line(aes(y = predichos, color = "Modelo GLS Desfasado"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Observado (Sensor)" = "black", "Modelo GLS Desfasado" = "firebrick")) +
  labs(
    title = "MACN-Ma (2020): Ajuste con Temperatura Exterior Desfasada / Acumulada",
    subtitle = paste("Correlación con predicho:", round(cor(df_model$temp_max, df_model$predichos), 3)),
    x = "Fecha",
    y = "Temperatura (°C)",
    color = "Serie"
  ) +
  theme_minimal()