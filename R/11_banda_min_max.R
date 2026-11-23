# ==============================================================================
# BANDA DEL DÍA: ENTRE QUÉ TEMPERATURAS VA A ESTAR LA SALA
# ==============================================================================
# Junta los pronósticos de la máxima y de la mínima interior (05, corrido con
# VARIABLE=max y con VARIABLE=min) y arma, para cada día, una banda
#
#     [ inferior de la mínima ; superior de la máxima ]
#
# que debería contener todas las temperaturas del día con un 95 % de
# confianza. Cada extremo es un límite de un solo lado al 97,5 %, así que,
# por Bonferroni, los dos se cumplen juntos al menos el 95 % de los días.
#
# Se evalúa fuera de la muestra (las predicciones del 05 ya lo son): la
# cobertura es el % de días en que la mínima observada quedó por encima del
# límite inferior y la máxima observada, por debajo del superior.
#
# Dos formas de calcular los límites:
#   teorico    predicho ± 1,96·EE, con el EE del modelo (normal).
#   calibrado  predicho + q·EE, donde q es el cuantil 2,5 % (mínima) o 97,5 %
#              (máxima) de los errores estandarizados de los días ANTERIORES
#              ya evaluados de la misma serie. No mira el futuro. Mientras
#              haya menos de 30 días anteriores se usa 1,96.
#
# Modelos (ver 05):
#   sin_hoy_corr  el que se puede usar a la mañana: exterior de ayer hacia
#                 atrás y temperatura de la sala de los días anteriores.
#   con_hoy_corr  agrega la TMAX/TMIN de hoy. En la práctica se reemplazaría
#                 por el pronóstico del SMN para hoy; es el techo.
#
# Necesita resultados/tablas/pronostico_predicciones.csv y
# pronostico_predicciones_min.csv.
# ==============================================================================

library(tidyverse)

dir_tablas <- "resultados/tablas"
dir_figs   <- "resultados/figuras"
min_calib  <- 30   # días anteriores necesarios para calibrar

leer <- function(sufijo) {
  ruta <- file.path(dir_tablas, paste0("pronostico_predicciones", sufijo, ".csv"))
  if (!file.exists(ruta)) stop("Falta ", ruta, ": correr antes R/05_modelo_pronostico.R",
                               if (sufijo == "_min") " con VARIABLE=min")
  read_csv(ruta, show_col_types = FALSE)
}

modelos <- c("sin_hoy_corr", "con_hoy_corr")

pred <- inner_join(
  leer("") %>% filter(Modelo %in% modelos) %>%
    select(Serie, Modelo, Date, obs_max = observado, pred_max = predicho, ee_max = ee),
  leer("_min") %>% filter(Modelo %in% modelos) %>%
    select(Serie, Modelo, Date, obs_min = observado, pred_min = predicho, ee_min = ee),
  by = c("Serie", "Modelo", "Date")
) %>% arrange(Serie, Modelo, Date)

# Cuantil de los errores estandarizados de los días anteriores (sin el actual)
cuantil_previo <- function(z, p, defecto) {
  map_dbl(seq_along(z), function(i) {
    previos <- z[seq_len(i - 1)]
    if (length(previos) < min_calib) defecto else unname(quantile(previos, p))
  })
}

banda <- pred %>%
  group_by(Serie, Modelo) %>%
  mutate(
    z_max = (obs_max - pred_max) / ee_max,
    z_min = (obs_min - pred_min) / ee_min,
    q_max = cuantil_previo(z_max, 0.975, qnorm(0.975)),
    q_min = cuantil_previo(z_min, 0.025, qnorm(0.025)),
    inf_teorico   = pred_min - qnorm(0.975) * ee_min,
    sup_teorico   = pred_max + qnorm(0.975) * ee_max,
    inf_calibrado = pred_min + q_min * ee_min,
    sup_calibrado = pred_max + q_max * ee_max
  ) %>%
  ungroup()

cobertura <- banda %>%
  pivot_longer(c(starts_with("inf_"), starts_with("sup_")),
               names_to = c(".value", "Metodo"), names_sep = "_") %>%
  group_by(Serie, Modelo, Metodo) %>%
  summarise(
    Dias = n(),
    Cobertura_banda = 100 * mean(obs_min >= inf & obs_max <= sup),
    Cobertura_min   = 100 * mean(obs_min >= inf),
    Cobertura_max   = 100 * mean(obs_max <= sup),
    Ancho_medio     = mean(sup - inf),
    Rango_obs_medio = mean(obs_max - obs_min),
    .groups = "drop"
  )

cat("\n=== COBERTURA DE LA BANDA (% de días; objetivo: 95 %) ===\n")
print(cobertura %>% mutate(across(where(is.double), ~ round(.x, 1))), n = Inf, width = Inf)
write_csv(cobertura, file.path(dir_tablas, "banda_cobertura.csv"))
write_csv(banda, file.path(dir_tablas, "banda_por_dia.csv"))

# ------------------------------------------------------------------------------
# Gráfico: la banda de una serie a lo largo del tiempo
# ------------------------------------------------------------------------------
serie_graf <- Sys.getenv("SERIE", "BA-Pb General")
datos_graf <- banda %>% filter(Serie == serie_graf, Modelo == "sin_hoy_corr")

if (nrow(datos_graf) > 0) {
  cob <- cobertura %>% filter(Serie == serie_graf, Modelo == "sin_hoy_corr", Metodo == "calibrado")
  fuera <- datos_graf %>% filter(obs_min < inf_calibrado | obs_max > sup_calibrado)

  g <- ggplot(datos_graf, aes(x = Date)) +
    geom_ribbon(aes(ymin = inf_calibrado, ymax = sup_calibrado), fill = "#2E6FD8", alpha = 0.18) +
    geom_linerange(aes(ymin = obs_min, ymax = obs_max), color = "gray25", linewidth = 0.5) +
    geom_linerange(data = fuera, aes(ymin = obs_min, ymax = obs_max), color = "#D2692E", linewidth = 0.9) +
    labs(
      title = paste0(serie_graf, ": la banda pronosticada a la mañana contiene el día ",
                     round(cob$Cobertura_banda), " % de las veces"),
      subtitle = paste0("Banda: mínima y máxima pronosticadas sin el SMN de hoy, con los días anteriores\n",
                        "de la sala (95 %, calibrada). Líneas: rango observado de cada día; en naranja, los que se salen."),
      x = NULL, y = "Temperatura interior (°C)"
    ) +
    theme_minimal(base_size = 16) +
    theme(panel.grid.minor = element_blank(), plot.title.position = "plot")

  print(g)
  ggsave(file.path(dir_figs, "fig_banda_min_max.png"), g, width = 13, height = 6.5, dpi = 300)
} else {
  message("AVISO: no hay predicciones para '", serie_graf, "'. No se hace el gráfico.")
}
