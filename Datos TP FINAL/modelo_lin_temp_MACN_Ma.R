library(readxl)
library(lubridate)
require(tidyverse)
require(patchwork)
require(emmeans)
require(ggplot2)
library(performance)
library(lmerTest)
library(nlme)
library(mgcv)

options(emmeans = list(emmeans = list(infer = c(TRUE, FALSE)),
                       contrast = list(infer = c(TRUE, TRUE))))

# Cargar datos meteorológicos SMN
datos_smn      <- read.delim('Datos TP FINAL/smn_datos_meteorologicos_desde_2021.lst')
datos_smn_1991 <- read_xlsx('Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx')

# Filtrar por la estación 87585 (Aeroparque / CABA)
datos_smn <- datos_smn[datos_smn$NRO_OMM == '  87585', ]
datos_smn$FECHA <- as.Date(datos_smn$FECHA, format = '%d/%m/%Y')

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$NRO_OMM == '87585', ]
datos_smn_1991$FECHA <- as.Date(datos_smn_1991$FECHA)

# Filtrar rango de años relevante (el archivo objetivo es del 2020)
datos_smn_1991 <- datos_smn_1991[datos_smn_1991$FECHA >= as.Date('2014-01-01'), ]
names(datos_smn_1991)[14] <- "NUB_TOT" 

smn <- rbind(datos_smn, datos_smn_1991)
smn$FECHA <- ymd(smn$FECHA)
smn <- smn %>% rename(Date = FECHA)

path_macn_ma <- "Datos TP FINAL/MACN-Ma/MACN-Ma-03_02_2020-22_12_2020.txt"

# Si el archivo está en la carpeta de trabajo actual, puedes usar ruta relativa:
# path_macn_ma <- "Datos TP FINAL/MACN-Ma/MACN-Ma-03_02_2020-22_12_2020.txt"

targets <- read.delim(path_macn_ma, sep = ";")

# Extraer temperatura máxima por día
targets <- targets %>%
  group_by(Date) %>%
  summarise(temp_max = max(Temp...C, na.rm = TRUE), .groups = "drop")

targets$Date <- dmy(targets$Date)

# Merge con datos del SMN
df <- targets %>% 
  left_join(smn, by = "Date")

# Ordenar por fecha y crear índice temporal 't' (contempla huecos en la serie)
df <- df[order(df$Date), ]
df$t <- as.integer(df$Date - min(df$Date))

# Asegurar formatos numéricos
df$TMAX     <- as.numeric(df$TMAX)
df$HUM_REL  <- as.numeric(df$HUM_REL)
df$temp_max <- as.numeric(gsub(",", ".", df$temp_max))


m0 <- gls(temp_max ~ TMAX * HUM_REL, 
          data = df, 
          correlation = corAR1(form = ~ t),
          na.action = na.omit)

summary(m0)

# Predicciones marginales con emmeans
efectos_humedad <- emmeans(m0, ~ TMAX * HUM_REL, 
                           data = df,
                           at = list(TMAX = seq(min(df$TMAX, na.rm = TRUE), 
                                                max(df$TMAX, na.rm = TRUE), 
                                                length.out = 50),
                                     HUM_REL = quantile(df$HUM_REL, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)))

df_pred <- as.data.frame(efectos_humedad)
df_pred$HUM_REL_f <- factor(round(df_pred$HUM_REL, 1), 
                            labels = c("Baja (Q25)", "Media (Q50)", "Alta (Q75)"))

# Gráfico 1: Efecto de TMAX según Humedad Relativa
ggplot(df, aes(x = TMAX, y = temp_max)) +
  geom_point(alpha = 0.5, color = "gray30") +
  geom_line(data = df_pred, aes(y = emmean, color = HUM_REL_f), linewidth = 1.2) +
  geom_ribbon(data = df_pred, aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = HUM_REL_f), 
              alpha = 0.15, color = NA) +
  labs(
    title = "MACN-Ma (2020): Efecto de TMAX sobre temp_max según la Humedad Relativa",
    subtitle = "Ajustado por modelo GLS con autocorrelación corAR1",
    x = "Temperatura Máxima SMN (TMAX °C)",
    y = "Temperatura Máxima Registrada (temp_max °C)",
    color = "Humedad Relativa",
    fill = "Humedad Relativa"
  ) +
  theme_minimal()

# Gráfico 2: Ajuste del Modelo GLS m0 en el tiempo
df_complete <- na.omit(df[, c("Date", "t", "temp_max", "TMAX", "HUM_REL")])
df_complete$predichos <- predict(m0, newdata = df_complete)

ggplot(df_complete, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Observado"), linewidth = 0.8) +
  geom_line(aes(y = predichos, color = "Modelo GLS"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Observado" = "black", "Modelo GLS" = "red")) +
  labs(
    title = "MACN-Ma (2020): Ajuste del Modelo GLS en el Tiempo",
    x = "Fecha",
    y = "Temperatura Máxima (°C)",
    color = "Leyenda"
  ) +
  theme_minimal()

df$yday <- yday(df$Date)

m_gam <- gamm(temp_max ~ s(yday, bs = "cc", k = 12) + s(TMAX) + HUM_REL,
              data = df,
              correlation = corAR1(form = ~ t))

summary(m_gam$gam)
plot(m_gam$gam, pages = 1)


m_simple <- gls(temp_max ~ TMAX, 
                data = df, 
                correlation = corAR1(form = ~ t),
                na.action = na.omit)

summary(m_simple)

# Comparación por criterio de información de Akaike
AIC(m0, m_simple)

# Predicciones del modelo simple
df_clean <- na.omit(df[, c("Date", "t", "temp_max", "TMAX")])
df_clean$pred <- predict(m_simple, newdata = df_clean)

# Gráfico 3: Serie temporal del modelo simple vs observado
ggplot(df_clean, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Observado (Sensor)"), linewidth = 0.8) +
  geom_line(aes(y = pred, color = "Modelo GLS (TMAX)"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Observado (Sensor)" = "black", "Modelo GLS (TMAX)" = "firebrick")) +
  labs(
    x = "Fecha", 
    y = "Temperatura (°C)", 
    color = "Serie",
    title = "MACN-Ma (2020): Ajuste del Modelo GLS Simple en el Tiempo"
  ) +
  theme_minimal()