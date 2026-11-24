# ==============================================================================
# Esquema de la validación con origen móvil (R/05_modelo_pronostico.R)
# ==============================================================================
# Dibuja, para MACN-Or Antigua (693 días), cómo crece el entrenamiento y qué
# bloque de 30 días se pronostica en cada iteración. No lee datos: usa las
# mismas reglas que R/05 (primera mitad, mínimo 40 días, bloques de 30).

library(tidyverse)

n           <- 693   # días de MACN-Or Antigua con datos completos
min_entrena <- 40
bloque      <- 30
n_ini       <- max(min_entrena, floor(n / 2))
cortes      <- seq(n_ini, n - 1, by = bloque)

tramos <- map_dfr(seq_along(cortes), function(i) {
  corte <- cortes[i]
  tibble(
    iteracion = i,
    tipo  = c("Entrenamiento", "Pronóstico (30 días)"),
    desde = c(1, corte + 1),
    hasta = c(corte, min(corte + bloque, n))
  )
})

g <- ggplot(tramos) +
  geom_segment(aes(x = desde, xend = hasta, y = iteracion, yend = iteracion, color = tipo),
               linewidth = 3.2) +
  geom_vline(xintercept = n_ini + 0.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = n_ini - 8, y = 0.1, size = 3.4, hjust = 1,
           label = "K se elige por AIC solo con el primer entrenamiento") +
  scale_color_manual(values = c("Entrenamiento" = "gray75",
                                "Pronóstico (30 días)" = "#2E6FD8")) +
  scale_y_reverse(breaks = seq_along(cortes), expand = expansion(add = c(1.2, 0.6))) +
  labs(x = "Día de la serie", y = "Iteración", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("resultados/figuras/fig_origen_movil.png", g, width = 8, height = 3.6, dpi = 300)
