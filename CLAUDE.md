# CLAUDE.md

Este archivo orienta a Claude Code (claude.ai/code) cuando trabaja con el código de este repositorio.

## Idioma

Todo el trabajo se hace en español: las respuestas, los comentarios del código, los nombres de variables, los mensajes de commit y la documentación.

## Proyecto

Trabajo práctico final de la materia universitaria "Casos de Estudio en Estadística". Es un análisis en R de datos de microclima (temperatura, humedad relativa y punto de rocío) registrados con dataloggers HOBO en salas de colecciones de museo, sobre todo del MACN. Esos datos se relacionan con el clima exterior del Servicio Meteorológico Nacional (SMN) y con las capturas de insectos en trampas. La pregunta central es cómo afectan las condiciones exteriores a las interiores de cada sala.

**Antes de trabajar, leer `docs/`:**
- `descripcion_tp.md`: objetivo del TP.
- `modelo_MACN_Ar.md`: el primer modelo y sus problemas.
- `resultados_y_siguientes_pasos.md`: resultados por sala, problemas encontrados en el código y la lista de próximos pasos. Es el estado actual del análisis.

No hay build, linter ni tests. Los scripts se corren de forma interactiva en RStudio, abriendo `Casos de Estudio en Estadística.Rproj`, con la raíz del repo como directorio de trabajo. Para correr un script completo: `Rscript <archivo>.R` desde la raíz del repo. Paquetes usados: tidyverse, lubridate, readxl, nlme, mgcv, emmeans, ggplot2, patchwork, scales, performance, lmerTest.

## Datos (`Datos TP FINAL/`)

- Hay una subcarpeta por colección o sitio: `BA-v`, `BA-Pb` y `MACN-*` (Ar, He, Ict, In, Ma, Or, Pv). MACN = Museo Argentino de Ciencias Naturales.
- Los `*.hobo` son exportaciones binarias crudas de HOBOware y R no las puede leer. Hay que usar las versiones `.txt`.
- Los `*.txt` son exportaciones HOBO ya limpias. Usan `;` como separador, **coma decimal** y BOM UTF-8. El encabezado es `"Date";"Time, GMT-03:00";"Temp, °C";"RH, %";"DewPt, °C"`, que R convierte en las columnas `Date`, `Time..GMT.03.00`, `Temp...C`, `RH...` y `DewPt...C`. Las fechas son `dd/mm/aa`. Hay que leerlos con `read.csv2`, o convertir a número **antes** de resumir. Con `read.delim(sep=";")` la temperatura queda como texto y `max()` la compara alfabéticamente: `"9,9"` le gana a `"23,1"`. Ese error ya produjo resultados incorrectos.
- Los nombres siguen el patrón `<SITIO>-<inicio dd_mm_aaaa>-<fin dd_mm_aaaa>[-<ubicación>].txt`, pero el rango del nombre no siempre es el real: `MACN-Or-03_02_2020-28_12_2020.txt` llega hasta el 28/12/2021. Una sala puede tener varios dataloggers (Ict colección/laboratorio, In arriba/abajo, Pv Anexo/Ameghino, Or antigua/Azara). La unidad de análisis es el datalogger.
- El primer y el último día de cada archivo están incompletos.
- `No HOBO/` tiene copias idénticas byte a byte de todos los `.txt` limpios, juntadas en una sola carpeta. Si se modifica un `.txt`, hay que actualizar las dos copias. Hay una diferencia de nombre: en `No HOBO/` está `BA_v-10_12_14-9_05_2015.txt` y en `BA-v/` está `BA_v-10_12_2014-9_05_2015.txt`.
- El intervalo de muestreo cambia entre archivos (15 min, 30 min, 2 h o 4 h), y a veces dentro de un mismo archivo. Eso afecta la máxima diaria.
- `general_trampas.csv` tiene las capturas de coleópteros y lepidópteros por trampa (filas) y fecha (columnas), entre 2019 y 2020.
- `smn_datos_meteorologicos_desde_2021.lst` tiene los datos diarios del SMN para todas las estaciones. Está separado por tabulaciones, con campos rellenados con espacios. Se usa la estación 87585 (Buenos Aires Observatorio): hay que filtrar con `trimws(NRO_OMM) == "87585"` y aplicar `as.numeric` a las variables.
- `smn_datos_meteorologicos_1991_2020.xlsx` (72 MB, commiteado en `64969b7`) tiene los datos diarios del SMN de 1991 a 2020. `LEAME_SMN.xlsx` es su descripción y está solo en forma local. Diferencias con el `.lst`: `NRO_OMM` es numérico, la última columna se llama `NUB_TOTAL` en vez de `NUB_TOT`, y las variables vienen como texto con `"S/D"` para "sin dato".

## Scripts

- `modelo_temperatura.R`: **el modelo principal.** Lee todos los `.txt` de `No HOBO/` y arma una serie por sala/subsector con `extraer_info_sala()`. A cada serie le ajusta `gls(temp_max ~ TMAX_prom3d + sin_year + cos_year, correlation = corAR1(form = ~ t))`, donde `TMAX_prom3d` es el promedio de la TMAX exterior de hoy y los 2 días anteriores, y seno/coseno modelan el ciclo anual. Devuelve un R² por serie y gráficos.
- `modelo_rezagos.R`: la alternativa al promedio de 3 días. Pone la TMAX de hoy y de hasta 7 días atrás como predictoras separadas, elige por AIC (ML) cuántos días usar y reporta el efecto acumulado (la suma de los coeficientes) con su EE. Usa rutas relativas y se corre solo, sin depender de otros scripts; tarda varios minutos por la lectura del xlsx. Los resultados están en la sección 3 de `docs/resultados_y_siguientes_pasos.md`.
- `exploracion_datos.R`: los gráficos de BA-v y de capturas, y `evaluar_ventanas_inercia()`, que compara promedios de 1 a 5 días por sala. Depende de objetos creados en `modelo_temperatura.R`.
- `modelo_lin_temp_MACN_Ar.R` y `Datos TP FINAL/modelo_lin_temp_MACN_Ma.R`: los primeros modelos de una sola sala (`temp_max ~ TMAX*HUM_REL`, GLS con AR1, más un GAM). Todavía tienen el error del `max()` sobre texto.
- `analisis_trampas.R`: a pesar del nombre, es el GLS estacional de MACN-Ma. Todavía no hay un análisis de trampas.

## Detalles a tener en cuenta

- Algunos scripts usan rutas absolutas (`~/Documents/UBA/...`, `D:/...`). Hay que usar rutas relativas al repo.
- Algunos scripts llaman a `install.packages()` al principio.
