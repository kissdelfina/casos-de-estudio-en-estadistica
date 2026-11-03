#y ~ normal

#x correlacionados por tiempo (+huecos)

#Usar gls para ver matriz de autocorrelacion empezamos
#meter como predictora el mes con spline cíclico -<si usamos GAM
#paquete de r lubridate

install.packages("readxl")

library(readxl)
library(lubridate)
require(tidyverse)
require(patchwork)
require(lubridate)
require(emmeans)
require(ggplot2)
options(emmeans = list(emmeans = list(infer = c(TRUE, FALSE)),
                       contrast = list(infer = c(TRUE, TRUE))))
library(performance)
library(lmerTest)
library(nlme)  

#estaciones <- read_excel('/home/Estudiante/Descargas/Datos-diarios-2021-2026-27102026/smn_estaciones_meteorologicas21-26.xlsx')
datos_smn <- read.delim('Datos TP FINAL/smn_datos_meteorologicos_desde_2021.lst')
datos_smn_1991 <- read_xlsx('Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx')

datos_smn <- datos_smn[datos_smn$NRO_OMM=='  87585',]
datos_smn$FECHA <- as.Date(datos_smn$FECHA, format='%d/%m/%Y')

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$NRO_OMM =='87585',]
datos_smn_1991$FECHA <- as.Date(datos_smn_1991$FECHA)

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$FECHA >= as.Date('2014-01-01'),]
names(datos_smn_1991)[14] <-paste("NUB_TOT") 

smn <- rbind(datos_smn,datos_smn_1991)


#MACN_Ar
targets <- read.delim('Datos TP FINAL/No HOBO/MACN_Ar-05_04_2016-23_06_2016.txt', sep = ";") 

targets <- targets %>%
  group_by(Date) %>%
  summarise(temp_max = max(Temp...C, na.rm = TRUE), .groups = "drop")

targets$Date <- dmy(targets$Date)
smn$FECHA <- ymd(smn$FECHA)

smn <- smn %>% rename(Date = FECHA)

df <- targets %>% 
  left_join(smn, by = "Date")

df <- df[order(df$Date), ]
df$t <- as.integer(df$Date - min(df$Date))   # índice temporal entero (admite días faltantes)
df$TMAX <- as.numeric(df$TMAX)
df$HUM_REL <- as.numeric(df$HUM_REL)
df$temp_max <- as.numeric(gsub(",", ".", df$temp_max))


m0 <- gls(temp_max ~ TMAX*HUM_REL, 
          data = df, 
          correlation = corAR1(form = ~ t))

summary(m0)

library(emmeans)
library(ggplot2)

# predicciones marginales para un rango de TMAX y 3 niveles de HUM_REL
efectos_humedad <- emmeans(m0, ~ TMAX * HUM_REL, 
                           data = df, # <--- Agregar esto
                           at = list(TMAX = seq(min(df$TMAX, na.rm=TRUE), 
                                                max(df$TMAX, na.rm=TRUE), 
                                                length.out = 50),
                                     HUM_REL = quantile(df$HUM_REL, probs = c(0.25, 0.50, 0.75), na.rm=TRUE)))

df_pred <- as.data.frame(efectos_humedad)

df_pred$HUM_REL_f <- factor(round(df_pred$HUM_REL, 1), 
                            labels = c("Baja (Q25)", "Media (Q50)", "Alta (Q75)"))

# 2. Graficar datos observados + rectas del modelo con intervalos de confianza
ggplot(df, aes(x = TMAX, y = temp_max)) +
  geom_point(alpha = 0.5, color = "gray30") +
  geom_line(data = df_pred, aes(y = emmean, color = HUM_REL_f), size = 1.2) +
  geom_ribbon(data = df_pred, aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = HUM_REL_f), 
              alpha = 0.15, color = NA) +
  labs(
    title = "Efecto de TMAX sobre temp_max según la Humedad Relativa",
    subtitle = "Ajustado por modelo GLS con autocorrelación corAR1",
    x = "Temperatura Máxima SMN (TMAX °C)",
    y = "Temperatura Máxima Registrada (temp_max °C)",
    color = "Humedad Relativa",
    fill = "Humedad Relativa"
  ) +
  theme_minimal()



# Calcular predicciones para los datos originales (excluyendo filas con NAs)
df_complete <- na.omit(df[, c("Date", "t", "temp_max", "TMAX", "HUM_REL")])
df_complete$predichos <- predict(m0, newdata = df_complete)

ggplot(df_complete, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Observado"), size = 0.8) +
  geom_line(aes(y = predichos, color = "Modelo GLS"), size = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Observado" = "black", "Modelo GLS" = "red")) +
  labs(
    title = "Ajuste del Modelo GLS en el Tiempo",
    x = "Fecha",
    y = "Temperatura Máxima (°C)",
    color = "Leyenda"
  ) +
  theme_minimal()


library(mgcv)

# Convertir la fecha a mes numérico o día del año (1 a 365)
df$yday <- yday(df$Date)

# Modelo GAM con spline cíclico para la época del año y autocorrelación AR1
m_gam <- gamm(temp_max ~ s(yday, bs = "cc", k = 12) + s(TMAX) + HUM_REL,
              data = df,
              correlation = corAR1(form = ~ t))

summary(m_gam$gam)
plot(m_gam$gam, pages = 1)

# No hay curvatura ni no-linealidad en el efecto de TMAX sobre temp_max.
# La época del año (yday) no aporta información adicional porque la variación estacional ya está capturada por TMAX.
# HUM_REL Su inclusión como término lineal simple no está sumando poder explicativo al modelo.
# El modelo explica el 45% de la variancia ($R^2_{\text{adj}} = 0.45$) con $n = 80$ observaciones.
# Dado que s(yday) no aporta al ajuste y s(TMAX) es totalmente lineal (edf = 1), no necesitas un GAM complejo ni un spline cíclico.

# Modelo solo con TMAX y estructura temporal AR(1)
m_simple <- gls(temp_max ~ TMAX, 
                data = df, 
                correlation = corAR1(form = ~ t),
                na.action = na.omit)

summary(m_simple)

AIC(m0, m_simple)

# Predicciones del modelo simple
df_clean <- na.omit(df[, c("Date", "t", "temp_max", "TMAX")])
df_clean$pred <- predict(m_simple, newdata = df_clean)

# Gráfico de observados vs predichos en el tiempo
ggplot(df_clean, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Observado (Sensor)"), size = 0.8) +
  geom_line(aes(y = pred, color = "Modelo GLS (TMAX)"), size = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Observado (Sensor)" = "black", "Modelo GLS (TMAX)" = "firebrick")) +
  labs(x = "Fecha", y = "Temperatura (°C)", color = "Serie",
       title = "Ajuste del Modelo GLS Simple en el Tiempo") +
  theme_minimal()
