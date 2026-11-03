library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)
library(ggplot2)
library(scales)

# smn
smn_1991 <- read_xlsx('Datos TP FINAL/smn_datos_meteorologicos_1991_2020.xlsx') %>%
  filter(NRO_OMM == '87585') %>%
  mutate(Date = as.Date(FECHA), TMAX = as.numeric(gsub(",", ".", TMAX))) %>%
  select(Date, TMAX)

smn_2021 <- read.delim('Datos TP FINAL/smn_datos_meteorologicos_desde_2021.lst') %>%
  filter(trimws(NRO_OMM) == '87585') %>%
  mutate(Date = as.Date(FECHA, format = '%d/%m/%Y'), TMAX = as.numeric(gsub(",", ".", TMAX))) %>%
  select(Date, TMAX)

smn_clean <- bind_rows(smn_1991, smn_2021) %>%
  distinct(Date, .keep_all = TRUE) %>%
  arrange(Date)

# salas y sub-sectores
extraer_info_sala <- function(filepath) {
  fname <- basename(filepath)
  fname_clean <- gsub("[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}", "", fname)
  fname_clean <- gsub("\\.txt$", "", fname_clean)
  
  # Extraer subsector si existe (Anexo, Ameghino, Coleccion, Laboratorio, Abajo, Arriba)
  subsector <- "General / Normal"
  if (grepl("Anexo", fname_clean, ignore.case = TRUE)) subsector <- "Anexo"
  if (grepl("Ameghino", fname_clean, ignore.case = TRUE)) subsector <- "Ameghino"
  if (grepl("coleccion", fname_clean, ignore.case = TRUE)) subsector <- "Colección"
  if (grepl("laboratorio", fname_clean, ignore.case = TRUE)) subsector <- "Laboratorio"
  if (grepl("abajo", fname_clean, ignore.case = TRUE)) subsector <- "Abajo"
  if (grepl("arriba", fname_clean, ignore.case = TRUE)) subsector <- "Arriba"
  
  # Extraer nombre base de la sala (ej: MACN_Pv, BA_v, MACN_Ict)
  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba)", "", fname_clean, ignore.case = TRUE)
  sala_base <- gsub("-", "_", sala_base)
  
  return(list(
    Sala_Base = sala_base,
    Subsector = subsector,
    Etiqueta_Unica = paste(sala_base, subsector, sep = " - ")
  ))
}

path_carpeta <- "~/Documents/UBA/Casos de Estudio en Estadística/Casos de Estudio en Estadística/Datos TP FINAL/No HOBO"
archivos_txt <- list.files(path = path_carpeta, pattern = "\\.txt$", full.names = TRUE)

datos_sensores <- map_dfr(archivos_txt, function(path) {
  df <- read.delim(path, sep = ";")
  if (!"Date" %in% colnames(df) || !"Temp...C" %in% colnames(df)) return(NULL)
  
  info <- extraer_info_sala(path)
  
  df %>%
    mutate(
      Date = as.Date(parse_date_time(Date, orders = c("dmy", "ymd", "mdy"))),
      temp_max = as.numeric(gsub(",", ".", Temp...C)),
      Sala_Base = info$Sala_Base,
      Subsector = info$Subsector,
      Etiqueta_Unica = info$Etiqueta_Unica
    ) %>%
    group_by(Sala_Base, Subsector, Etiqueta_Unica, Date) %>%
    summarise(temp_max = max(temp_max, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(Date) & !is.infinite(temp_max))
})

# gls
procesar_serie <- function(df_serie, smn_datos) {
  etiqueta <- unique(df_serie$Etiqueta_Unica)
  sala <- unique(df_serie$Sala_Base)
  subsector <- unique(df_serie$Subsector)
  
  grilla <- data.frame(
    Date = seq(min(df_serie$Date), max(df_serie$Date), by = "day"),
    Sala_Base = sala,
    Subsector = subsector,
    Etiqueta_Unica = etiqueta
  )
  
  df_comp <- grilla %>%
    left_join(smn_datos, by = "Date") %>%
    left_join(df_serie, by = c("Sala_Base", "Subsector", "Etiqueta_Unica", "Date")) %>%
    arrange(Date) %>%
    mutate(
      TMAX_prom3d = (TMAX + dplyr::lag(TMAX, 1) + dplyr::lag(TMAX, 2)) / 3,
      yday = yday(Date),
      sin_year = sin(2 * pi * yday / 365.25),
      cos_year = cos(2 * pi * yday / 365.25)
    )
  
  df_train <- df_comp %>%
    filter(!is.na(temp_max) & !is.na(TMAX_prom3d)) %>%
    mutate(t = as.integer(Date - min(Date)))
  
  if (nrow(df_train) < 30) return(NULL)
  
  modelo <- tryCatch({
    gls(
      temp_max ~ TMAX_prom3d + sin_year + cos_year,
      data = df_train,
      correlation = corAR1(form = ~ t)
    )
  }, error = function(e) NULL)
  
  if (is.null(modelo)) return(NULL)
  
  df_pred <- df_comp %>% filter(!is.na(TMAX_prom3d))
  df_pred$predicho <- predict(modelo, newdata = df_pred)
  
  r_val <- cor(df_train$temp_max, predict(modelo, newdata = df_train))
  
  list(
    predicciones = df_pred,
    metricas = data.frame(
      Etiqueta = etiqueta,
      Sala_Base = sala,
      Subsector = subsector,
      N_obs = nrow(df_train),
      R2 = round(r_val^2, 4),
      Correlacion = round(r_val, 4)
    )
  )
}

series_list <- split(datos_sensores, datos_sensores$Etiqueta_Unica)
resultados <- map(series_list, ~ procesar_serie(.x, smn_clean)) %>% compact()

metricas_todas <- map_dfr(resultados, "metricas")
predicciones_todas <- map_dfr(resultados, "predicciones")

# salas y sub-sectores
g_general <- ggplot(predicciones_todas, aes(x = Date)) +
  geom_line(aes(y = temp_max, color = "Observado"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = predicho, color = "Modelo GLS"), linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("Observado" = "black", "Modelo GLS" = "firebrick")) +
  facet_wrap(~ Etiqueta_Unica, scales = "free", ncol = 3) +
  scale_x_date(date_labels = "%b %Y") +
  labs(
    title = "Ajuste del Modelo GLS en Todos los Sectores ('No HOBO')",
    subtitle = "Escala de tiempo liberada independientemente por sector",
    x = "Fecha",
    y = "Temperatura (°C)",
    color = "Serie"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(g_general)

# MACN_Pv
pred_pv <- predicciones_todas %>% filter(Sala_Base == "MACN_Pv")

g_pv <- ggplot(pred_pv, aes(x = Date, group = Subsector, color = Subsector)) +
  geom_line(aes(y = temp_max), linewidth = 0.8, na.rm = TRUE) +
  geom_line(aes(y = predicho), linewidth = 0.6, linetype = "dashed", alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(
    title = "Comparación Térmica Interna en MACN_Pv",
    subtitle = "Línea continua = Mediciones Reales | Línea punteada = Modelo GLS",
    x = "Fecha",
    y = "Temperatura (°C)",
    color = "Subsector"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(g_pv)