library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tidyr)

BA_v_2018_9 <- read.csv2("Datos TP FINAL/BA_v-21_05_2018-27_12_2019.txt")
BA_v_2020 <- read.csv2("Datos TP FINAL/BA_v-02_03_2020-22_12_2020.txt")

# Función para procesar y generar el gráfico apilado
plot_hobo <- function(data, titulo, breaks_x = "3 months") {
  df <- data %>%
    rename(Fecha = Date, Hora = Time..GMT.03.00, Temp = Temp...C, RH = RH..., DewPt = DewPt...C) %>%
    mutate(Datetime = dmy_hms(paste(Fecha, Hora)))
  
  p1 <- ggplot(df, aes(x = Datetime)) +
    geom_line(aes(y = Temp, color = "Temp (°C)"), linewidth = 0.6) +
    geom_line(aes(y = DewPt, color = "Punto Rocío (°C)"), linewidth = 0.6, linetype = "dashed") +
    scale_color_manual(values = c("Temp (°C)" = "#D95F02", "Punto Rocío (°C)" = "#7570B3")) +
    scale_x_datetime(date_labels = "%b %Y", expand = c(0.01, 0.01)) +
    labs(title = titulo, y = "Temp (°C)", x = NULL, color = NULL) +
    theme_minimal(11) +
    theme(legend.position = "top", legend.justification = "right", legend.margin = margin(b = -10),
          plot.title = element_text(face = "bold"), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), panel.grid.minor = element_blank())
  
  p2 <- ggplot(df, aes(x = Datetime, y = RH)) +
    geom_line(color = "#1B9E77", linewidth = 0.6) +
    scale_x_datetime(date_labels = "%b %Y", date_breaks = breaks_x, expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(y = "RH (%)", x = "Fecha / Hora") +
    theme_minimal(11) +
    theme(panel.grid.minor = element_blank())
  
  return(p1 / p2)
}

# Generar ambos gráficos
g_2018_9 <- plot_hobo(BA_v_2018_9, "Microclima HOBO (2018 - 2019)", "3 months")
g_2020   <- plot_hobo(BA_v_2020,   "Microclima HOBO (2020)",        "1 month")

# Visualizar
g_2018_9
g_2020

# --------------------------------------------------------------------------------------------------
# coleopteros
# 1. Dataset de capturas (BA-v)
BA_v_capturas <- data.frame(
  Fecha = dmy(c("15/04/19", "29/04/19", "16/05/19", "03/06/19", "11/11/19", 
                "25/11/19", "10/12/19", "27/12/19", "18/02/20", "01/04/20", "10/04/20")),
  Total_Capturas = c(20, 13, 2, 0, 7, 89, 58, 16, 115, 19, 2)
)

# 2. Promedios diarios del HOBO (2018-2019 y 2020)
hobo_diario <- bind_rows(BA_v_2018_9, BA_v_2020) %>%
  rename(Fecha = Date, Temp = Temp...C, RH = RH..., DewPt = DewPt...C) %>%
  mutate(Fecha = dmy(Fecha)) %>%
  group_by(Fecha) %>%
  summarise(
    Temp = mean(Temp, na.rm = TRUE),
    RH = mean(RH, na.rm = TRUE),
    DewPt = mean(DewPt, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Cruzar datos por fecha
df_corr <- inner_join(BA_v_capturas, hobo_diario, by = "Fecha")

# Graficos de dispersion

# Gráfico 1: Humedad Relativa vs. Capturas
g_rh <- ggplot(df_corr, aes(x = RH, y = Total_Capturas)) +
  geom_point(size = 3.5, color = "#1B9E77") +
  geom_smooth(method = "lm", color = "#0D5C43", fill = "#B2E2E2", alpha = 0.3) +
  labs(
    title = "Humedad Relativa vs. Capturas",
    x = "Humedad Relativa Promedio (%)",
    y = "Total Coleópteros"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Gráfico 2: Punto de Rocío vs. Capturas
g_dewpt <- ggplot(df_corr, aes(x = DewPt, y = Total_Capturas)) +
  geom_point(size = 3.5, color = "#7570B3") +
  geom_smooth(method = "lm", color = "#3F3A78", fill = "#CBC9E2", alpha = 0.3) +
  labs(
    title = "Punto de Rocío vs. Capturas",
    x = "Punto de Rocío Promedio (°C)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Gráfico 3: Temperatura vs. Capturas
g_temp <- ggplot(df_corr, aes(x = Temp, y = Total_Capturas)) +
  geom_point(size = 3.5, color = "#D95F02") +
  geom_smooth(method = "lm", color = "#963B00", fill = "#FDBB84", alpha = 0.3) +
  labs(
    title = "Temperatura vs. Capturas",
    x = "Temperatura Promedio (°C)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# 4. Combinar los tres gráficos horizontalmente
(g_rh | g_dewpt | g_temp)


# Interaccion
g_interaccion_2d <- ggplot(df_corr, aes(x = Temp, y = RH)) +
  geom_point(aes(size = Total_Capturas, fill = Total_Capturas), shape = 21, color = "black", alpha = 0.85) +
  scale_fill_gradient(low = "#FFEDA0", high = "#E31A1C", name = "Total Capturas") +
  scale_size_continuous(range = c(3, 10), name = "Total Capturas") +
  geom_text(aes(label = Total_Capturas), vjust = -1.2, size = 3, fontface = "bold") +
  labs(
    title = "Interacción Climática: Temperatura vs. Humedad Relativa",
    subtitle = "El tamaño y color representan la cantidad de coleópteros capturados",
    x = "Temperatura Promedio (°C)",
    y = "Humedad Relativa Promedio (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(g_interaccion_2d)
