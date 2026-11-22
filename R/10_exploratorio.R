# ==============================================================================
# ANÁLISIS EXPLORATORIO
# ==============================================================================
# Caracteriza las 14 series (una por datalogger) que se usan en los modelos:
#   - período, cantidad de archivos, intervalo de muestreo y días con datos;
#   - media, DE y rango de la máxima diaria interior;
#   - DE del cambio de la máxima de un día al siguiente, en el interior y en el
#     exterior (SMN) para los mismos días. A diferencia de la DE de la serie,
#     no depende del ciclo anual, así que mide directamente cuánto amortigua
#     cada sala las variaciones del exterior;
#   - línea de tiempo de los períodos medidos;
#   - resumen de las capturas por sala y taxón.
#
# Las series se arman igual que en 03_modelo_rezagos.R: se descartan el primer
# y el último día de cada archivo, y MACN-Or se separa en Antigua y Azara.
# ==============================================================================

library(lubridate)
library(tidyverse)

dir_tablas  <- "resultados/tablas"
dir_figuras <- "resultados/figuras"
dir.create(dir_tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figuras, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. SMN: TMAX diaria de la estación 87585 (caché que genera el 01 o el 03)
# ------------------------------------------------------------------------------
smn_1991 <- readRDS("datos/procesados/smn_87585_1991_2020.rds")
smn_2021 <- read.delim("datos/smn/smn_datos_meteorologicos_desde_2021.lst") %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
            TMAX = suppressWarnings(as.numeric(trimws(TMAX))))
smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 2. HOBO: máxima diaria por archivo y serie
# ------------------------------------------------------------------------------
nombre_serie <- function(filepath) {
  fname <- basename(filepath)
  limpio <- gsub("[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}", "", fname)
  limpio <- gsub("\\.txt$", "", limpio)

  sector <- ""
  if (grepl("Anexo", limpio, ignore.case = TRUE)) sector <- "Anexo"
  if (grepl("Ameghino", limpio, ignore.case = TRUE)) sector <- "Ameghino"
  if (grepl("coleccion", limpio, ignore.case = TRUE)) sector <- "Colección"
  if (grepl("laboratorio", limpio, ignore.case = TRUE)) sector <- "Laboratorio"
  if (grepl("abajo", limpio, ignore.case = TRUE)) sector <- "Abajo"
  if (grepl("arriba", limpio, ignore.case = TRUE)) sector <- "Arriba"
  if (grepl("azara", limpio, ignore.case = TRUE) || grepl("04-01_2021", fname)) sector <- "Azara"
  if (grepl("MACN-Or-03_02_2020", fname)) sector <- "Antigua"
  if (grepl("^MACN.Pv", limpio) && sector == "") sector <- "General"

  sala <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba|azara|antigua)", "",
               limpio, ignore.case = TRUE)
  sala <- gsub("_", "-", sala)
  trimws(paste(sala, sector))
}

leer_hobo <- function(path) {
  d <- read.csv2(path, fileEncoding = "UTF-8-BOM")
  lecturas <- tibble(Date = dmy(d[[1]]),
                     hora = d[[2]],
                     temp = as.numeric(d[[3]])) %>%
    filter(!is.na(Date), !is.na(temp))

  # Intervalo de muestreo más frecuente del archivo, en minutos
  momentos <- dmy_hms(paste(d[[1]], d[[2]]), quiet = TRUE)
  if (all(is.na(momentos))) momentos <- dmy_hm(paste(d[[1]], d[[2]]), quiet = TRUE)
  difs <- round(as.numeric(diff(momentos), units = "mins"))
  intervalo <- as.numeric(names(which.max(table(difs))))

  lecturas %>%
    group_by(Date) %>%
    summarise(temp_max = max(temp), .groups = "drop") %>%
    filter(Date > min(Date), Date < max(Date)) %>%   # primer y último día incompletos
    mutate(Archivo = basename(path), Serie = nombre_serie(path), Intervalo_min = intervalo)
}

archivos  <- list.files("datos/hobo", pattern = "\\.txt$", full.names = TRUE)
datos     <- map_dfr(archivos, leer_hobo) %>% left_join(smn, by = "Date")

# ------------------------------------------------------------------------------
# 3. Tabla descriptiva por serie
# ------------------------------------------------------------------------------
formato_intervalo <- function(x) {
  x <- sort(unique(x))
  paste(ifelse(x >= 60, paste0(x / 60, " h"), paste0(x, " min")), collapse = " y ")
}

descriptiva <- datos %>%
  arrange(Serie, Date) %>%
  group_by(Serie) %>%
  mutate(
    # Cambio de un día al siguiente, solo entre días consecutivos
    consecutivo = as.numeric(Date - lag(Date)) == 1,
    dif_int = if_else(consecutivo, temp_max - lag(temp_max), NA_real_),
    dif_ext = if_else(consecutivo, TMAX - lag(TMAX), NA_real_)
  ) %>%
  summarise(
    Archivos   = n_distinct(Archivo),
    Desde      = min(Date),
    Hasta      = max(Date),
    Intervalo  = formato_intervalo(Intervalo_min),
    N_dias     = n(),
    Media_int  = mean(temp_max),
    DE_int     = sd(temp_max),
    Min_int    = min(temp_max),
    Max_int    = max(temp_max),
    DE_cambio_int = sd(dif_int, na.rm = TRUE),
    DE_cambio_ext = sd(dif_ext, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Cociente_cambio = DE_cambio_int / DE_cambio_ext) %>%
  arrange(desc(Cociente_cambio))

cat("\n=== DESCRIPTIVA POR SERIE ===\n")
print(descriptiva %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf, width = Inf)
write_csv(descriptiva, file.path(dir_tablas, "exploratorio_series.csv"))

# ------------------------------------------------------------------------------
# 4. Línea de tiempo de los períodos medidos
# ------------------------------------------------------------------------------
# Un segmento por tramo continuo de registros (los huecos de más de 7 días cortan)
tramos <- datos %>%
  arrange(Serie, Date) %>%
  group_by(Serie) %>%
  mutate(tramo = cumsum(c(TRUE, diff(Date) > 7))) %>%
  group_by(Serie, tramo) %>%
  summarise(inicio = min(Date), fin = max(Date), .groups = "drop") %>%
  left_join(descriptiva %>% select(Serie, Cociente_cambio), by = "Serie") %>%
  mutate(Serie = fct_reorder(Serie, inicio, .fun = min, .desc = TRUE))

g_tiempo <- ggplot(tramos, aes(y = Serie)) +
  geom_segment(aes(x = inicio, xend = fin, yend = Serie), linewidth = 3.5,
               color = "steelblue", lineend = "butt") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Las salas se midieron en años distintos, con pocas superposiciones",
       subtitle = "Períodos con registros de cada datalogger, de 2014 a 2025",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

print(g_tiempo)
ggsave(file.path(dir_figuras, "fig_linea_tiempo.png"), g_tiempo, width = 8, height = 4.5, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Capturas por sala y taxón
# ------------------------------------------------------------------------------
archivos_capturas <- list.files("datos/trampas", recursive = TRUE, full.names = TRUE,
                                pattern = "_(coleopteros|lepidopteros)[.]csv$")
capturas <- map_dfr(archivos_capturas, function(f) {
  read.csv(f) %>%
    mutate(fecha = as.Date(fecha),
           sala  = gsub("_.*", "", basename(f)),
           taxon = str_to_title(basename(dirname(f))))
})

resumen_capturas <- capturas %>%
  group_by(sala, taxon) %>%
  summarise(Revisiones = n(), Total = sum(cantidad),
            Mediana = median(cantidad), Maximo = max(cantidad),
            Desde = min(fecha), Hasta = max(fecha), .groups = "drop")

cat("\n=== CAPTURAS POR SALA Y TAXÓN ===\n")
print(resumen_capturas)
write_csv(resumen_capturas, file.path(dir_tablas, "exploratorio_capturas.csv"))
