#y ~ normal

#x correlacionados por tiempo (+huecos)

#Usar gls para ver matriz de autocorrelacion empezamos
#meter como predictora el mes con spline cíclico -<si usamos GAM
#paquete de r lubridate

install.packages("readxl")

library(readxl)
library(lubridate)

#estaciones <- read_excel('/home/Estudiante/Descargas/Datos-diarios-2021-2026-27102026/smn_estaciones_meteorologicas21-26.xlsx')
datos_smn <- read.delim('/home/Estudiante/Descargas/Datos-diarios-2021-2026-27102026/smn_datos_meteorologicos_desde_2021.lst')
datos_smn_1991 <- read_xlsx('/home/Estudiante/Descargas/smn_datos_meteorologicos_1991_2020/smn_datos_meteorologicos_1991_2020.xlsx')

datos_smn <- datos_smn[datos_smn$NRO_OMM=='  87585',]
datos_smn$FECHA <- as.Date(datos_smn$FECHA, format='%d/%m/%Y')

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$NRO_OMM =='87585',]
datos_smn_1991$FECHA <- as.Date(datos_smn_1991$FECHA)

datos_smn_1991 <- datos_smn_1991[datos_smn_1991$FECHA >= as.Date('2014-01-01'),]
names(datos_smn_1991)[14] <-paste("NUB_TOT") 

smn <- rbind(datos_smn,datos_smn_1991)

targets <- 




