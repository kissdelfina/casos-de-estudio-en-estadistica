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
datos_smn <- read.delim('D:/Casos de Estadistica/TP/smn_datos_meteorologicos_desde_2021.lst')
datos_smn_1991 <- read_xlsx('D:/Casos de Estadistica/TP/smn_datos_meteorologicos_1991_2020.xlsx')

datos_smn <- datos_smn[datos_smn$NRO_OMM=='  87585',]
datos_smn$FECHA <- as.Date(datos_smn$FECHA, format='%d/%m/%Y')

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$NRO_OMM =='87585',]
datos_smn_1991$FECHA <- as.Date(datos_smn_1991$FECHA)

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$FECHA >= as.Date('2014-01-01'),]
names(datos_smn_1991)[14] <-paste("NUB_TOT") 

smn <- rbind(datos_smn,datos_smn_1991)


#MACN_Ar
targets <- read.delim('D:/Casos de Estadistica/TP/Datos TP FINAL/No HOBO/MACN_Ar-05_04_2016-23_06_2016.txt', sep = ";") 

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
