# ==============================================================================
# EFECTO DEL DÍA DE LA SEMANA
# ==============================================================================
# En el diagnóstico de residuos (06_diagnostico_residuos.R) aparece un pico de
# autocorrelación a los 7 días en algunas salas (MACN-Ar, Ict Colección, BA-v).
# El clima no tiene ciclo semanal, pero el uso del edificio sí: en días hábiles
# hay gente, luces, puertas abiertas y quizás climatización.
#
# Para cada serie, con la ventana K que eligió 03_modelo_rezagos.R y errores
# AR(1), se comparan tres modelos:
#
#   sin_dias  temp_max ~ TMAX_l0 + ... + TMAX_lK + sin_year + cos_year
#   habil     ... + habil        (1 de lunes a viernes, 0 sábado y domingo)
#   dia       ... + dia          (factor de 7 niveles, referencia: domingo)
#
# Se informa el ΔAIC contra el modelo sin días, cuántos °C más tiene
# la sala en un día hábil, el coeficiente de cada día respecto del domingo, si
# cambia el efecto acumulado de la TMAX exterior y si baja la autocorrelación
# de los residuos a los 7 días.
#
# Los feriados todavía no se tienen en cuenta: cuentan como días hábiles, lo
# que, si acaso, achica el efecto estimado.
#
# Necesita resultados/tablas/rezagos_resumen.csv, que genera el 03.
# ==============================================================================

library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)

max_rezago <- 14  # el mismo que en 03_modelo_rezagos.R

# ------------------------------------------------------------------------------
# 1. Datos (igual que 03_modelo_rezagos.R)
# ------------------------------------------------------------------------------
ruta_cache <- "datos/procesados/smn_87585_1991_2020.rds"
if (file.exists(ruta_cache)) {
  smn_1991 <- readRDS(ruta_cache)
} else {
  smn_1991 <- read_xlsx("datos/smn/smn_datos_meteorologicos_1991_2020.xlsx") %>%
    filter(NRO_OMM == 87585) %>%
    transmute(Date = as.Date(FECHA),
              TMAX = suppressWarnings(as.numeric(gsub(",", ".", TMAX))))
  dir.create(dirname(ruta_cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(smn_1991, ruta_cache)
  gc()
}

smn_2021 <- read.delim("datos/smn/smn_datos_meteorologicos_desde_2021.lst") %>%
  filter(trimws(NRO_OMM) == "87585") %>%
  transmute(Date = as.Date(trimws(FECHA), format = "%d/%m/%Y"),
            TMAX = suppressWarnings(as.numeric(trimws(TMAX))))

smn <- bind_rows(smn_1991, smn_2021) %>%
  filter(!is.na(Date)) %>%
  distinct(Date, .keep_all = TRUE) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  arrange(Date)

for (k in 0:max_rezago) {
  smn[[paste0("TMAX_l", k)]] <- dplyr::lag(smn$TMAX, k)
}
smn$TMAX_prom3d <- (smn$TMAX_l0 + smn$TMAX_l1 + smn$TMAX_l2) / 3

extraer_info_sala <- function(filepath) {
  fname <- basename(filepath)
  fname <- gsub("[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}[-_]\\d{1,2}[-_]\\d{1,2}[-_]\\d{2,4}", "", fname)
  fname <- gsub("\\.txt$", "", fname)

  subsector <- "General / Normal"
  if (grepl("Anexo", fname, ignore.case = TRUE)) subsector <- "Anexo"
  if (grepl("Ameghino", fname, ignore.case = TRUE)) subsector <- "Ameghino"
  if (grepl("coleccion", fname, ignore.case = TRUE)) subsector <- "Colección"
  if (grepl("laboratorio", fname, ignore.case = TRUE)) subsector <- "Laboratorio"
  if (grepl("abajo", fname, ignore.case = TRUE)) subsector <- "Abajo"
  if (grepl("arriba", fname, ignore.case = TRUE)) subsector <- "Arriba"

  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba)", "", fname, ignore.case = TRUE)
  sala_base <- gsub("-", "_", sala_base)

  paste(sala_base, subsector, sep = " - ")
}

leer_hobo <- function(path) {
  d <- read.csv2(path, fileEncoding = "UTF-8-BOM")
  tibble(Date = dmy(d[[1]]), temp = as.numeric(d[[3]])) %>%
    filter(!is.na(Date), !is.na(temp)) %>%
    group_by(Date) %>%
    summarise(temp_max = max(temp), .groups = "drop") %>%
    filter(Date > min(Date), Date < max(Date)) %>%
    mutate(Archivo = basename(path), Etiqueta = extraer_info_sala(path))
}

archivos <- list.files("datos/hobo", pattern = "\\.txt$", full.names = TRUE)
datos_hobo <- map_dfr(archivos, leer_hobo)

# Ventana elegida por el 03 para cada serie
ventanas <- read_csv("resultados/tablas/rezagos_resumen.csv", show_col_types = FALSE) %>%
  transmute(Etiqueta, K = as.integer(sub("0 a ", "", Rezagos_mejor)))

# ------------------------------------------------------------------------------
# 2. Funciones
# ------------------------------------------------------------------------------
nombres_rezagos <- function(k) paste0("TMAX_l", 0:k)
dias_es <- c("Lunes", "Martes", "Miércoles", "Jueves", "Viernes", "Sábado", "Domingo")

formula_k <- function(k, extra = NULL) {
  reformulate(c(nombres_rezagos(k), "sin_year", "cos_year", extra), response = "temp_max")
}

ajustar <- function(f, df) {
  tryCatch(gls(f, data = df, correlation = corAR1(form = ~ t), method = "ML"),
           error = function(e) {
             message("  AVISO: falló el ajuste (", conditionMessage(e), ")")
             NULL
           })
}

# Autocorrelación de los residuos normalizados a los 7 días, sobre una grilla
# diaria para respetar los huecos
acf_7 <- function(m, t) {
  grilla <- rep(NA_real_, max(t) + 1)
  grilla[t + 1] <- resid(m, type = "normalized")
  acf(grilla, lag.max = 7, na.action = na.pass, plot = FALSE)$acf[8]
}

analizar_serie <- function(df_serie) {
  etiqueta <- df_serie$Etiqueta[1]
  k <- ventanas$K[ventanas$Etiqueta == etiqueta]
  if (length(k) == 0 || anyDuplicated(df_serie$Date) > 0) return(NULL)

  df <- df_serie %>%
    left_join(smn, by = "Date") %>%
    mutate(sin_year = sin(2 * pi * yday(Date) / 365.25),
           cos_year = cos(2 * pi * yday(Date) / 365.25),
           # wday con week_start = 1: 1 = lunes, ..., 7 = domingo (no depende del idioma)
           num_dia  = wday(Date, week_start = 1),
           habil    = as.integer(num_dia <= 5),
           dia      = relevel(factor(dias_es[num_dia], levels = dias_es), ref = "Domingo")) %>%
    drop_na(temp_max, all_of(nombres_rezagos(max_rezago)), TMAX_prom3d) %>%
    arrange(Date) %>%
    mutate(t = as.integer(Date - min(Date)))

  message("Serie: ", etiqueta, " (", nrow(df), " días, K = ", k, ")")

  modelos <- list(
    sin_dias = ajustar(formula_k(k), df),
    habil    = ajustar(formula_k(k, "habil"), df),
    dia      = ajustar(formula_k(k, "dia"), df)
  )
  if (any(map_lgl(modelos, is.null))) return(NULL)

  nm <- nombres_rezagos(k)
  comparacion <- imap_dfr(modelos, function(m, nombre) {
    tibble(Etiqueta = etiqueta, Modelo = nombre, AIC = AIC(m),
           Efecto_acum = sum(coef(m)[nm]),
           EE_acum = sqrt(sum(vcov(m)[nm, nm])),
           phi = coef(m$modelStruct$corStruct, unconstrained = FALSE),
           ACF_7 = acf_7(m, df$t))
  }) %>%
    mutate(dAIC = AIC - AIC[Modelo == "sin_dias"])

  m_h <- modelos$habil
  efecto_habil <- tibble(
    Etiqueta = etiqueta, N_dias = nrow(df), K = k,
    Habil = unname(coef(m_h)["habil"]),
    EE_habil = unname(sqrt(diag(vcov(m_h))["habil"])),
    # Test de Wald conjunto para los 6 coeficientes del día (vía razón de verosimilitud)
    p_dia = anova(modelos$sin_dias, modelos$dia)$`p-value`[2]
  )

  m_d <- modelos$dia
  nms <- paste0("dia", setdiff(dias_es, "Domingo"))
  coefs_dia <- tibble(
    Etiqueta = etiqueta,
    Dia = c(setdiff(dias_es, "Domingo"), "Domingo"),
    Coef = c(unname(coef(m_d)[nms]), 0),
    EE = c(unname(sqrt(diag(vcov(m_d))[nms])), 0)
  )

  list(comparacion = comparacion, efecto_habil = efecto_habil, coefs_dia = coefs_dia)
}

# ------------------------------------------------------------------------------
# 3. Correr y resumir
# ------------------------------------------------------------------------------
series     <- split(datos_hobo, datos_hobo$Etiqueta)
resultados <- map(series, analizar_serie) %>% compact()

comparacion  <- map_dfr(resultados, "comparacion")
efecto_habil <- map_dfr(resultados, "efecto_habil")
coefs_dia    <- map_dfr(resultados, "coefs_dia")

resumen <- efecto_habil %>%
  left_join(comparacion %>% select(Etiqueta, Modelo, dAIC, Efecto_acum, ACF_7) %>%
              pivot_wider(names_from = Modelo, values_from = c(dAIC, Efecto_acum, ACF_7)),
            by = "Etiqueta") %>%
  select(Etiqueta, N_dias, K, Habil, EE_habil, p_dia,
         dAIC_habil, dAIC_dia,
         Efecto_acum_sin_dias, Efecto_acum_dia, ACF_7_sin_dias, ACF_7_dia) %>%
  arrange(dAIC_dia)

cat("\n=== DÍA DE LA SEMANA (ΔAIC contra el modelo sin días; negativo = mejora) ===\n")
print(resumen %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)

cat("\n=== COEFICIENTE DE CADA DÍA RESPECTO DEL DOMINGO (°C) ===\n")
print(coefs_dia %>% select(-EE) %>% mutate(Coef = round(Coef, 2)) %>%
        pivot_wider(names_from = Dia, values_from = Coef), n = Inf, width = Inf)

dir.create("resultados/tablas", recursive = TRUE, showWarnings = FALSE)
dir.create("resultados/figuras", recursive = TRUE, showWarnings = FALSE)
write_csv(resumen,     "resultados/tablas/dia_semana_resumen.csv")
write_csv(comparacion, "resultados/tablas/dia_semana_modelos.csv")
write_csv(coefs_dia,   "resultados/tablas/dia_semana_coeficientes.csv")

# Coeficiente de cada día respecto del domingo, con IC del 95 %
datos_graf <- coefs_dia %>%
  mutate(Serie = Etiqueta %>%
           str_replace("MACN_Pv - General / Normal", "MACN-Pv General") %>%
           str_replace(" - General / Normal", "") %>%
           str_replace(" - ", " ") %>%
           str_replace("_", "-"),
         Dia = factor(Dia, levels = dias_es))

g_dias <- ggplot(datos_graf, aes(Dia, Coef)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_pointrange(aes(ymin = Coef - 1.96 * EE, ymax = Coef + 1.96 * EE),
                  color = "firebrick", size = 0.3) +
  geom_line(aes(group = 1), color = "firebrick", alpha = 0.5) +
  facet_wrap(~ Serie, ncol = 4) +
  scale_x_discrete(labels = c("L", "Ma", "Mi", "J", "V", "S", "D")) +
  labs(title = "Diferencia de la máxima interior respecto del domingo",
       subtitle = "Modelo de rezagos con AR(1) más el día de la semana. Barras: IC 95 %",
       x = NULL, y = "°C respecto del domingo") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(g_dias)
ggsave("resultados/figuras/fig_dia_semana.png", g_dias, width = 13, height = 9, dpi = 200)
