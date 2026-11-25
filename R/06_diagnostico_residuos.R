# ==============================================================================
# DIAGNÓSTICO DE RESIDUOS DEL MODELO DE REZAGOS: ¿ALCANZA EL AR(1)?
# ==============================================================================
# El modelo de 03_modelo_rezagos.R supone errores AR(1). Si quedara
# autocorrelación sin modelar, los errores estándar del efecto acumulado
# serían demasiado chicos. Para cada serie, con la ventana K que eligió el 03:
#
#   1. Se ajusta el mismo modelo con tres estructuras de error:
#        AR(1)       corAR1(form = ~ t)              (la que usamos)
#        AR(2)       corARMA(form = ~ t, p = 2)
#        ARMA(1,1)   corARMA(form = ~ t, p = 1, q = 1)
#      y se comparan por AIC (ML).
#   2. Se mira la autocorrelación de los residuos normalizados del AR(1)
#      (prueba de Ljung-Box a 7 y 14 días). Los residuos normalizados ya
#      descuentan la correlación modelada: si el AR(1) alcanza, deberían
#      parecerse a ruido blanco.
#   3. Se compara el efecto acumulado y su EE entre el AR(1) y la mejor
#      estructura, y se vuelve a elegir K con la mejor estructura, para ver si
#      cambian las conclusiones.
#
# Lee los datos igual que 03_modelo_rezagos.R (mismas series, mismas filas) y
# necesita resultados/tablas/rezagos_resumen.csv, que genera ese script.
# ==============================================================================

library(readxl)
library(lubridate)
library(tidyverse)
library(nlme)

max_rezago <- 14  # el mismo que en 03_modelo_rezagos.R

# ------------------------------------------------------------------------------
# 1. Datos (copiado de 03_modelo_rezagos.R)
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
  # MACN-Or: son dos dataloggers con fechas superpuestas (colección antigua y
  # Azara), igual que en 03_modelo_rezagos.R
  if (grepl("azara", fname, ignore.case = TRUE) || grepl("04-01_2021", basename(filepath))) subsector <- "Azara"
  if (grepl("MACN-Or-03_02_2020", basename(filepath))) subsector <- "Antigua"

  sala_base <- gsub("[-_](Anexo|Ameghino|coleccion|laboratorio|abajo|arriba|azara|antigua)", "", fname, ignore.case = TRUE)
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
formula_k <- function(k) reformulate(c(nombres_rezagos(k), "sin_year", "cos_year"),
                                     response = "temp_max")

estructuras <- list(
  "AR(1)"     = function() corAR1(form = ~ t),
  "AR(2)"     = function() corARMA(form = ~ t, p = 2, q = 0),
  "ARMA(1,1)" = function() corARMA(form = ~ t, p = 1, q = 1)
)

# Se prueba primero el optimizador por defecto (nlminb) y, si falla, optim.
ajustar <- function(k, estructura, df, avisar = TRUE) {
  for (opt in c("nlminb", "optim")) {
    m <- tryCatch(
      gls(formula_k(k), data = df, correlation = estructuras[[estructura]](),
          method = "ML", control = glsControl(msMaxIter = 200, opt = opt)),
      error = function(e) {
        if (avisar) message("  ", estructura, " con ", opt, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(m)) return(m)
  }
  NULL
}

efecto_acumulado <- function(m, k) {
  nm <- nombres_rezagos(k)
  c(efecto = sum(coef(m)[nm]), ee = sqrt(sum(vcov(m)[nm, nm])))
}

# Ljung-Box sobre los residuos normalizados. La serie tiene huecos (días sin
# datos y archivos de años distintos), así que los residuos se ubican en una
# grilla diaria y las autocorrelaciones se calculan solo con pares de días
# que existen. fitdf descuenta los parámetros de la estructura de error.
ljung_box <- function(res, t, lag, fitdf) {
  grilla <- rep(NA_real_, max(t) + 1)
  grilla[t + 1] <- res
  r <- acf(grilla, lag.max = lag, na.action = na.pass, plot = FALSE)$acf[-1]
  n <- sum(!is.na(grilla))
  q <- n * (n + 2) * sum(r^2 / (n - seq_len(lag)))
  pchisq(q, df = lag - fitdf, lower.tail = FALSE)
}

acf_residuos <- function(res, t, lag = 14) {
  grilla <- rep(NA_real_, max(t) + 1)
  grilla[t + 1] <- res
  acf(grilla, lag.max = lag, na.action = na.pass, plot = FALSE)$acf[-1]
}

diagnosticar_serie <- function(df_serie) {
  etiqueta <- df_serie$Etiqueta[1]
  k <- ventanas$K[ventanas$Etiqueta == etiqueta]
  if (length(k) == 0 || anyDuplicated(df_serie$Date) > 0) return(NULL)

  df <- df_serie %>%
    left_join(smn, by = "Date") %>%
    mutate(sin_year = sin(2 * pi * yday(Date) / 365.25),
           cos_year = cos(2 * pi * yday(Date) / 365.25)) %>%
    drop_na(temp_max, all_of(nombres_rezagos(max_rezago)), TMAX_prom3d) %>%
    arrange(Date) %>%
    mutate(t = as.integer(Date - min(Date)))

  message("Serie: ", etiqueta, " (", nrow(df), " días, K = ", k, ")")

  modelos <- map(names(estructuras), ~ ajustar(k, .x, df)) %>% set_names(names(estructuras))
  fallidos <- names(modelos)[map_lgl(modelos, is.null)]
  if (length(fallidos) > 0) message("  AVISO: no convergió ", paste(fallidos, collapse = ", "))
  modelos <- compact(modelos)
  if (!"AR(1)" %in% names(modelos)) {
    message("  AVISO: no convergió el AR(1). Se saltea la serie.")
    return(NULL)
  }

  params_cor <- c("AR(1)" = 1, "AR(2)" = 2, "ARMA(1,1)" = 2)
  comparacion <- imap_dfr(modelos, function(m, nombre) {
    ef  <- efecto_acumulado(m, k)
    res <- resid(m, type = "normalized")
    tibble(Etiqueta = etiqueta, Estructura = nombre, AIC = AIC(m),
           Efecto = ef[["efecto"]], EE = ef[["ee"]],
           p_LB7  = ljung_box(res, df$t, 7,  params_cor[[nombre]]),
           p_LB14 = ljung_box(res, df$t, 14, params_cor[[nombre]]),
           Parametros_cor = paste(round(coef(m$modelStruct$corStruct, unconstrained = FALSE), 3),
                                  collapse = " / "))
  }) %>%
    mutate(dAIC = AIC - min(AIC))

  mejor <- comparacion$Estructura[which.min(comparacion$AIC)]

  # Si gana otra estructura, ¿cambia la ventana elegida?
  k_mejor_estr <- k
  if (mejor != "AR(1)") {
    aics <- map_dbl(0:max_rezago, function(kk) {
      m <- ajustar(kk, mejor, df, avisar = FALSE)
      if (is.null(m)) Inf else AIC(m)
    })
    k_mejor_estr <- (0:max_rezago)[which.min(aics)]
  }

  acfs <- imap_dfr(modelos, function(m, nombre) {
    tibble(Etiqueta = etiqueta, Estructura = nombre, Rezago = 1:14,
           ACF = acf_residuos(resid(m, type = "normalized"), df$t),
           n = nrow(df))
  })

  list(comparacion = comparacion %>% mutate(K_AR1 = k, Mejor = mejor,
                                            K_con_mejor = k_mejor_estr),
       acfs = acfs)
}

# ------------------------------------------------------------------------------
# 3. Correr y resumir
# ------------------------------------------------------------------------------
series     <- split(datos_hobo, datos_hobo$Etiqueta)
resultados <- map(series, diagnosticar_serie) %>% compact()

comparacion <- map_dfr(resultados, "comparacion")
acfs        <- map_dfr(resultados, "acfs")

cat("\n=== COMPARACIÓN DE ESTRUCTURAS DE ERROR (K fijo, el elegido con AR(1)) ===\n")
print(comparacion %>%
        select(Etiqueta, Estructura, dAIC, Parametros_cor, Efecto, EE, p_LB7, p_LB14) %>%
        mutate(across(where(is.double), ~ round(.x, 3))),
      n = Inf, width = Inf)

resumen <- comparacion %>%
  group_by(Etiqueta) %>%
  summarise(
    K_AR1          = first(K_AR1),
    Mejor          = first(Mejor),
    dAIC_AR1       = dAIC[Estructura == "AR(1)"],
    Efecto_AR1     = Efecto[Estructura == "AR(1)"],
    EE_AR1         = EE[Estructura == "AR(1)"],
    Efecto_mejor   = Efecto[Estructura == first(Mejor)],
    EE_mejor       = EE[Estructura == first(Mejor)],
    K_con_mejor    = first(K_con_mejor),
    p_LB14_AR1     = p_LB14[Estructura == "AR(1)"],
    p_LB14_mejor   = p_LB14[Estructura == first(Mejor)],
    .groups = "drop"
  ) %>%
  mutate(Cambio_EE_pct = 100 * (EE_mejor / EE_AR1 - 1)) %>%
  arrange(desc(dAIC_AR1))

cat("\n=== RESUMEN: ¿CAMBIA ALGO SI NO USAMOS AR(1)? ===\n")
print(resumen %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)

dir.create("resultados/tablas", recursive = TRUE, showWarnings = FALSE)
dir.create("resultados/figuras", recursive = TRUE, showWarnings = FALSE)
write_csv(comparacion, "resultados/tablas/residuos_estructuras.csv")
write_csv(resumen,     "resultados/tablas/residuos_resumen.csv")

# ACF de los residuos normalizados: AR(1) contra la mejor estructura.
# Banda: ±1,96/√n, lo esperable si fueran ruido blanco.
datos_acf <- acfs %>%
  left_join(resumen %>% select(Etiqueta, Mejor), by = "Etiqueta") %>%
  filter(Estructura == "AR(1)" | Estructura == Mejor) %>%
  mutate(Serie = Etiqueta %>%
           str_replace("MACN_Pv - General / Normal", "MACN-Pv General") %>%
           str_replace(" - General / Normal", "") %>%
           str_replace(" - ", " ") %>%
           str_replace("_", "-"))

g_acf <- ggplot(datos_acf, aes(Rezago, ACF, color = Estructura)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_hline(aes(yintercept = 1.96 / sqrt(n)), linetype = "dashed", color = "gray60") +
  geom_hline(aes(yintercept = -1.96 / sqrt(n)), linetype = "dashed", color = "gray60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ Serie, ncol = 4) +
  scale_x_continuous(breaks = c(1, 7, 14)) +
  scale_color_manual(values = c("AR(1)" = "firebrick", "AR(2)" = "steelblue",
                                "ARMA(1,1)" = "darkgreen")) +
  labs(title = "Autocorrelación de los residuos normalizados",
       subtitle = "Si la estructura de error alcanza, los puntos quedan dentro de las líneas punteadas (±1,96/√n)",
       x = "Rezago (días)", y = "ACF", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

print(g_acf)
ggsave("resultados/figuras/fig_residuos_acf.png", g_acf, width = 13, height = 9, dpi = 200)
