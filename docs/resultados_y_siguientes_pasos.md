# Resultados y siguientes pasos

Estado al 17/11/2026. Las secciones 1 y 2 vienen de la revisión del commit `64969b7` ("modelo_temperatura tiene gls con tiempo defasado"). Las secciones 3 a 5 se agregaron después: el modelo con rezagos separados (`R/03_modelo_rezagos.R`), el pronóstico (`R/05_modelo_pronostico.R`), problemas nuevos en los datos, estructuras de correlación y el diagnóstico de residuos (`R/06_diagnostico_residuos.R`).

## 1. MACN-Ar: el modelo original

`gls(temp_max ~ TMAX * HUM_REL, corAR1)` con 80 días (05/04 al 23/06/2016).

| | Original | Con el `max()` corregido |
|---|---|---|
| TMAX | −0,05 (p = 0,90) | 0,05 (p = 0,88) |
| HUM_REL | −0,05 (p = 0,52) | −0,04 (p = 0,59) |
| TMAX × HUM_REL | 0,003 (p = 0,46) | 0,002 (p = 0,53) |
| φ (AR1) | 0,87 | 0,89 |

No da nada significativo, por dos razones:

- **`max()` sobre texto:** la temperatura se lee como texto, y `"9,928"` le gana a `"23,1"`. El 16 % de los máximos diarios estaba mal.
- **La interacción sin centrar:** las correlaciones entre los coeficientes llegan a −0,985, así que el modelo no puede separar los términos. Sacarla no empeora el ajuste (AIC 250,7 contra 252,3).

Con el `max()` corregido y sin la interacción:

| Predictora | Coeficiente | AIC |
|---|---|---|
| TMAX del día (+ HUM_REL) | 0,25 (p < 0,001); HUM_REL p = 0,75 | — |
| TMAX del día | 0,27 | 243,6 |
| TMAX del día anterior | 0,19 | 254,0 |
| TMAX de hace 2 días | 0,01 | 262,5 |
| **Promedio de TMAX de 3 días** | **0,76** | **222,9** |

- La humedad exterior no aporta.
- El edificio responde con inercia: con el promedio de 3 días el efecto es mucho mayor y φ baja de 0,89 a 0,56.
- Sin AR(1) el ajuste es mucho peor (AIC 309 contra 270).

## 2. Todas las salas: `modelo_temperatura.R`

El grupo generalizó el modelo a todos los `.txt` de `No HOBO/`:

```r
gls(temp_max ~ TMAX_prom3d + sin_year + cos_year, correlation = corAR1(form = ~ t))
```

La pendiente de `TMAX_prom3d` dice cuántos °C sube la máxima interior por cada °C que sube el promedio de 3 días de la máxima exterior, una vez descontado el ciclo anual.

| Serie | Días | Período | Pendiente (EE) | φ | R² |
|---|---|---|---|---|---|
| MACN-Ar | 78 | abr–jun 2016 | **0,62** (0,10) | 0,59 | 0,83 |
| MACN-Pv Anexo | 366 | 2019 | **0,54** (0,03) | 0,48 | 0,95 |
| BA-v | 1275 | 2014–2020 | **0,40** (0,01) | 0,83 | 0,97 |
| MACN-Ict Laboratorio | 398 | 2022–2023 | 0,22 (0,01) | 0,87 | 0,98 |
| MACN-Ma | 334 | 2016 y 2020 | 0,19 (0,01) | 0,93 | 0,97 |
| MACN-Pv (sin subsector) | 305 | 2019 | 0,17 (0,02) | 0,83 | 0,96 |
| MACN-In Arriba | 52 | mar–may 2025 | 0,15 (0,03) | 0,97 | 0,94 |
| MACN-Pv Ameghino | 269 | 2019 | 0,12 (0,01) | 0,91 | 0,95 |
| MACN-In Abajo | 52 | mar–may 2025 | 0,12 (0,04) | 0,91 | 0,94 |
| MACN-Ict Colección | 712 | 2022–2025 | 0,10 (0,01) | 0,96 | 0,95 |
| MACN-He | 196 | 2014 y 2016 | 0,09 (0,01) | 0,97 | 0,81 |
| BA-Pb | 325 | 2022–2023 | 0,09 (0,01) | 0,96 | 0,96 |
| MACN-Or | — | 2020–2021 | **falla** | — | — |

Todas las pendientes son significativas, pero hay dos grupos:

- **Salas expuestas** (0,4–0,6): MACN-Ar, MACN-Pv Anexo y BA-v. Siguen bastante al clima exterior.
- **Salas aisladas** (≈ 0,1–0,2): el resto. El exterior les llega muy amortiguado. Además, un φ ≈ 0,9 indica que su temperatura depende sobre todo de su propia historia.

### Ventana de inercia por sala (`exploracion_datos.R`)

Para cada sala se compararon promedios de TMAX de 1 a 5 días, eligiendo la ventana con menor AIC:

| Ventana ganadora | Series |
|---|---|
| 2 días | MACN-Ar, MACN-Pv Anexo, BA-v |
| 3 días | MACN-Ict Laboratorio, MACN-Ma |
| 4 días | MACN-Pv Ameghino, MACN-Pv (sin subsector), BA-Pb |
| 5 días (el máximo probado) | MACN-He, MACN-Ict Colección, MACN-In Abajo y Arriba, MACN-Or |

Es coherente con las pendientes: **las salas expuestas responden en ~2 días y las aisladas tardan 4–5 días o más**. Como 5 días es el máximo probado, en las salas aisladas la inercia real puede ser mayor.

Hay que tomar la comparación con cuidado por tres razones:

- `gls` ajusta por REML, y los AIC de REML no son comparables entre modelos con distintas predictoras. Hay que reajustar con `method = "ML"`.
- Cuando el AR(1) falla, el script ajusta sin AR(1) sin avisar. Eso pasa en MACN-Or, por las fechas repetidas, así que su resultado no es comparable con el de las otras salas.
- El R² casi no cambia entre ventanas, porque lo domina el ciclo anual.

### Problemas encontrados

1. **MACN-Or falla sin aviso.** `MACN-Or-03_02_2020-28_12_2020.txt` llega en realidad hasta el 28/12/2021, y se superpone con `MACN_Or-04-01_2021-28_12_2021.txt`. Son dos dataloggers distintos (por los `.hobo`, colección "antigua" y colección "Azara") que quedan en una misma serie con fechas repetidas, y el `corAR1` no lo admite. `compact()` descarta la serie sin decir nada.
2. **Series que mezclan archivos.** `extraer_info_sala()` junta todos los archivos sin subsector de una sala:
   - MACN-Pv "General" une dos archivos de 2019 que pueden ser dataloggers distintos.
   - MACN-Ma une 12 días de 2016 con 2020.
   - MACN-He une 2014 con 2016.
   - BA-v une 4 archivos de 2014 a 2020.

   Hay que confirmar si cada caso es el mismo datalogger en la misma ubicación.
3. **El R² está inflado.** Con series de casi un año, el seno y el coseno solos ya explican casi toda la variación. Un R² de 0,95 no significa que el exterior explique el 95 %. Para medir el aporte del exterior hay que comparar contra un modelo sin `TMAX_prom3d`.
4. **Series cortas** (MACN-Ar con 78 días, MACN-In con 52): el ciclo anual no se puede estimar con menos de un año. En ese caso el seno y el coseno funcionan como una tendencia y le "roban" efecto a TMAX.
5. **Períodos distintos:** las salas se midieron en años diferentes, así que comparar pendientes mezcla el efecto de la sala con el del año.
6. **Días incompletos:** el primer y el último día de cada archivo tienen menos registros, así que la máxima de esos días está subestimada.

### Otros problemas en el commit

- `modelo_lin_temp_MACN_Ar.R` y `Datos TP FINAL/modelo_lin_temp_MACN_Ma.R` siguen calculando `max()` sobre texto. `modelo_temperatura.R` y `analisis_trampas.R` ya lo hacen bien.
- `analisis_trampas.R` es en realidad el GLS estacional de MACN-Ma y no analiza las trampas. Además usa una ruta absoluta propia (`~/Documents/UBA/...`), igual que `modelo_temperatura.R`.
- ~~`Datos TP FINAL/plot_zoom_png?width=1470&height=923` tiene un `?` en el nombre.~~ Ya se borró en `5f2b7ba`.
- El xlsx del SMN 1991–2020 quedó commiteado (72 MB). Funciona, pero cualquier cambio a ese archivo va a hacer crecer el repo.

## 3. Rezagos separados: `R/03_modelo_rezagos.R`

En vez del promedio de 3 días, cada día anterior entra como predictora separada:

```r
gls(temp_max ~ TMAX_l0 + TMAX_l1 + ... + TMAX_lK + sin_year + cos_year,
    correlation = corAR1(form = ~ t), method = "ML")
```

`TMAX_lk` es la máxima exterior de hace k días. Para cada serie se prueba K = 0 a 14 y se elige el de menor AIC. (Hasta el 16/11 se probaba hasta 7, y en 8 de 12 series ganaba el máximo; con 14 ninguna serie elige el máximo.) Todos los modelos se ajustan por ML y sobre las mismas filas, así que los AIC son comparables. Además, el script:

- calcula los rezagos sobre la serie completa del SMN, así que no pierde los primeros días de cada serie;
- descarta el primer y el último día de cada archivo, que están incompletos;
- avisa cuando una serie falla, en vez de descartarla en silencio;
- usa rutas relativas y `read.csv2`.

**Efecto acumulado:** es la suma de los coeficientes, β₀ + … + β_K. Indica cuánto sube la máxima interior si la exterior sube 1 °C y se mantiene así K+1 días seguidos. Equivale también al efecto total, sumado sobre los días siguientes, de un solo día 1 °C más caluroso. Se puede comparar directamente con la pendiente del promedio de 3 días, que es el caso particular β₀ = β₁ = β₂. Los rezagos están muy correlacionados entre sí (r ≈ 0,85 entre días consecutivos), así que cada β por separado es inestable, pero la suma no: su error estándar es √(1ᵀV1).

| Serie | Rezagos elegidos | Ventanas con ΔAIC < 2 | Efecto acumulado (EE) | Pendiente del promedio de 3 días | ΔAIC contra el promedio |
|---|---|---|---|---|---|
| BA-v | 0 a 7 | 7–8 | 0,63 (0,02) | 0,40 | −657 |
| MACN-Ar | 0 a 2 | 1–3 | 0,60 (0,10) | 0,61 | −5 |
| MACN-Pv Anexo | 0 a 2 | 2–3 | 0,55 (0,02) | 0,54 | −152 |
| MACN-Pv Ameghino | 0 a **13** | 13–14 | 0,50 (0,04) | 0,13 | −55 |
| MACN-Ict Laboratorio | 0 a **11** | 10–13 | 0,50 (0,04) | 0,22 | −116 |
| MACN-Ma | 0 a **10** | 7–11 | 0,47 (0,04) | 0,20 | −93 |
| MACN-Ict Colección | 0 a **10** | 10–12 | 0,39 (0,03) | 0,11 | −104 |
| MACN-In Arriba | 0 a 5 | 4–6 | 0,28 (0,05) | 0,14 | −7 |
| MACN-Pv (sin subsector) | 0 a 5 | 5–6 | 0,27 (0,03) | 0,17 | −22 |
| MACN-He | 0 a 9 | 7–10 | 0,27 (0,03) | 0,10 | −24 |
| BA-Pb | 0 a 7 | 7–9 | 0,25 (0,03) | 0,09 | −74 |
| MACN-In Abajo | 0 a 4 | 4–6 | 0,22 (0,04) | 0,13 | −8 |
| MACN-Or | se saltea por las fechas repetidas | — | — | — | — |

- **Los rezagos separados ajustan mejor que el promedio en todas las salas.** El test de razón de verosimilitud rechaza los pesos iguales (p < 0,05) en 8 de 12 series. Las excepciones son MACN-He, Pv Ameghino, Pv sin subsector e In Arriba, con p entre 0,07 y 0,10.
- **El promedio de 3 días subestima el efecto en las salas aisladas:** el efecto acumulado es entre 2 y 4 veces la pendiente del promedio. Con ventanas de hasta 14 días, **la diferencia entre salas expuestas y aisladas está sobre todo en la velocidad, no en la magnitud**: las expuestas reciben 0,55–0,63 °C por °C exterior en 2 a 7 días, y varias aisladas reciben 0,39–0,50 en 10 a 13 días. Las que reciben menos (0,22–0,28: In, Pv sin subsector, He, BA-Pb) siguen siendo las más aisladas.
- **Las salas responden con perfiles distintos.** En las expuestas (MACN-Ar, Pv Anexo) casi todo el efecto (86–93 %) llega el mismo día y el siguiente. BA-v también tiene su mayor peso el mismo día, pero el efecto sigue hasta 7 días después. En las aisladas, solo el 24–42 % llega en los dos primeros días, y el resto se reparte a lo largo de una a dos semanas.
- **Ventanas ambiguas:** en MACN-Ma (7–11) y MACN-He (7–10), el AIC es casi plano en un rango amplio, así que la duración exacta de la inercia es incierta. En Pv Ameghino, en cambio, cualquier ventana de menos de 13 días ajusta claramente peor.
- El efecto acumulado es condicional al ciclo anual: mide la respuesta a días más calurosos o más fríos de lo normal para la época, no todo el vínculo entre interior y exterior.
- Tablas completas en `resultados/tablas/rezagos_resumen.csv`, `rezagos_coeficientes.csv` y `rezagos_aic.csv` (AIC y BIC de cada ventana).

**Robustez con BIC.** El BIC castiga cada día agregado con log(n) (de 4 a 7 según la serie) en vez de 2, así que elige ventanas iguales o más cortas:

| Serie | Ventana AIC → BIC | Efecto acumulado AIC → BIC |
|---|---|---|
| BA-v | 7 → 6 | 0,63 → 0,61 |
| MACN-Ar | 2 → 1 | 0,60 → 0,51 |
| MACN-Pv Anexo | 2 → 2 | 0,55 → 0,55 |
| MACN-Pv Ameghino | 13 → **6** | 0,50 → **0,29** |
| MACN-Ict Laboratorio | 11 → **6** | 0,50 → **0,38** |
| MACN-Ma | 10 → 6 | 0,47 → 0,38 |
| MACN-Ict Colección | 10 → 7 | 0,39 → 0,31 |
| MACN-In Arriba | 5 → 4 | 0,28 → 0,24 |
| MACN-Pv (sin subsector) | 5 → 5 | 0,27 → 0,27 |
| MACN-He | 9 → 5 | 0,27 → 0,18 |
| BA-Pb | 7 → 6 | 0,25 → 0,22 |
| MACN-In Abajo | 4 → 4 | 0,22 → 0,22 |

- **Se mantiene con los dos criterios:** las salas rápidas (Pv Anexo, Ar) usan 1–2 días y las aisladas, 5–7 o más. La diferencia de velocidad no depende del criterio.
- **Depende del criterio:** la cola larga de las salas aisladas (los días 7 a 13). Cada uno de esos días aporta un coeficiente chico que el BIC no considera suficiente, pero juntos suman entre 0,1 y 0,2 °C/°C. Por eso el efecto acumulado de las aisladas queda en un rango: Ameghino entre 0,29 y 0,50, Ict Laboratorio entre 0,38 y 0,50. La afirmación "las aisladas reciben casi tanto como las expuestas, pero más lento" vale con AIC y no con BIC.
- **Evidencia fuera de la muestra a favor de las ventanas largas:** al pasar el pronóstico de un máximo de 7 a 14 días, las series que eligieron ventanas más largas pronosticaron mejor en 8 de 10 casos (por ejemplo, Ameghino sin la TMAX de hoy: RMSE de 0,89 a 0,68 °C; Or Antigua con hoy: de 0,60 a 0,51). La cola larga ayuda a predecir días nuevos, que es lo que el AIC intenta optimizar.
- Un modelo con forma impuesta para los rezagos (por ejemplo, decaimiento geométrico, con la máxima interior de ayer como predictora) estimaría la cola con un solo parámetro y evitaría tener que elegir entre AIC y BIC.

### Pronóstico con y sin la TMAX de hoy: `R/05_modelo_pronostico.R`

Se compara, fuera de la muestra, el error de pronosticar la máxima interior con la TMAX de hoy (el techo, porque a la mañana no se conoce), sin ella (solo la de ayer hacia atrás) y sin datos del exterior (solo la época del año). Origen móvil: se entrena con la primera mitad de la serie (mínimo 40 días), se predicen los 30 días siguientes, se reentrena y se repite. La ventana se elige por AIC, entre 0 y 14 días, solo con el primer período de entrenamiento. Pérdida = 100·(RMSE sin hoy / RMSE con hoy − 1).

| Serie | Días evaluados | K con hoy | K sin hoy | RMSE con hoy | RMSE sin hoy | Pérdida sin hoy | Mejora contra solo la época |
|---|---|---|---|---|---|---|---|
| BA-v | 635 | 12 | 1 | 0,91 | 1,51 | +66 % | 12 % |
| MACN-Pv Anexo | 183 | 2 | 1 | 1,01 | 1,57 | +55 % | 27 % |
| MACN-Ict Laboratorio | 199 | 11 | 7 | 0,57 | 0,71 | +24 % | 34 % |
| MACN-Pv Ameghino | 135 | 11 | 11 | 0,58 | 0,68 | +17 % | 40 % |
| MACN-Ict Colección | 355 | 11 | 9 | 1,01 | 1,18 | +17 % | 19 % |
| MACN-He | 97 | 7 | 5 | 1,17 | 1,31 | +12 % | 17 % |
| BA-Pb | 163 | 8 | 8 | 0,69 | 0,76 | +10 % | 31 % |
| MACN-Or Antigua | 347 | 14 | 10 | 0,51 | 0,56 | +10 % | 38 % |
| MACN-Ma | 166 | 7 | 5 | 0,87 | 0,84 | −3 % | 29 % |
| MACN-Pv (sin subsector) | 152 | 5 | 5 | 0,72 | 0,68 | −6 % | 22 % |
| MACN-Or Azara | 179 | 11 | 11 | 0,54 | 0,48 | −12 % | 35 % |

MACN-Ar (38 días evaluados) y MACN-In (12) quedan en la tabla del CSV pero no son comparables: el primer entrenamiento es tan corto que el seno y el coseno extrapolan mal.

- **Solo dos salas dependen de verdad de la TMAX de hoy:** BA-v (+66 %) y Pv Anexo (+55 %). Son las de respuesta rápida.
- **Cambio respecto de la corrida con máximo 7:** Pv Ameghino bajaba de +49 % a +17 %. Con 7 días, el modelo sin hoy elegía solo 2 rezagos y no captaba la inercia; con 14 elige 11 y pronostica casi igual de bien sin la TMAX de hoy. Ict Laboratorio pasa de +28 % a +24 %, Ict Colección de +9 % a +17 % y Or Antigua de +3 % a +10 %.
- En todas las salas, la TMAX de los días anteriores mejora el pronóstico entre 12 % y 40 % respecto de usar solo la época del año.
- La ventana que elige el pronóstico puede diferir de la del 03 porque se elige con la mitad de los datos (por ejemplo, Ameghino: 11 contra 13).
- `predict()` de `gls` no usa la corrección AR(1): no aprovecha la temperatura interior de ayer. Un modelo que sí la use (por ejemplo, con la máxima interior de ayer como predictora) debería pronosticar mejor, sobre todo en las salas con φ alto.

## 4. Problemas nuevos en los datos

1. **El intervalo de muestreo cambia entre series y dentro de una misma serie:** 15 min (MACN-He 2014), 30 min (BA-v 2014–2017, MACN-Ar, MACN-He 2016, MACN-In), 2 h (MACN-Ma 2016) y 4 h (el resto). Con 6 lecturas por día la máxima diaria queda subestimada, y eso sesga la comparación entre salas. `BA_v-16_09_2016-17_05_2017.txt` pasa de 30 min a 4 h el 10/01/2017, y la serie BA-v completa mezcla los dos intervalos.
2. **MACN-He podría juntar dos ubicaciones.** Los `.hobo` se llaman `Escalera_Emergencia` e `Interior`, y en invierno de 2014 la temperatura va de 16,5 a 22,7 °C, mientras que en 2016 va de 13,1 a 17,0 °C.
3. **BA-v 2016–17 podría cambiar de datalogger el 10/01/2017.** En esa fecha cambia el intervalo, y hay dos `.hobo` del mismo período: `BA 10-01-2017 al 17-05-2017` y `BA interior compactos 10-01-17 al 17-05-17`.
4. **MACN-Ict Colección une dos archivos separados por unos 11 meses** (02/2023 a 01/2024).
5. **Hay `.hobo` que no tienen su `.txt`.** El caso más importante es `BA_11-12-2019 al 02-03-2020.hobo`: los datos de BA-v terminan el 27/12/2019, así que falta el período de la captura más grande (18/02/2020, 154 coleópteros en BA). Tampoco tienen `.txt`: BA interior compactos 2015–16, Ict laboratorio 2024–25, Pv 2013, Ma "tipos" y Or "armario".
6. **Nombres de archivo:**
   - En `No HOBO/` se llama `BA_v-10_12_14-9_05_2015.txt` y en `BA-v/` se llama `BA_v-10_12_2014-9_05_2015.txt`. El contenido es idéntico.
   - `MACN_Ar-abril-junio_2026.hobo` dice 2026 y los datos son de 2016.
   - El `.hobo` de BA-v 2014 se llama "duda fechas".
7. **`general_trampas.csv` está mal formado.** El encabezado tiene 25 fechas (50 columnas), pero las filas tienen 5 columnas más con datos (por ejemplo `17,17,0`) que no tienen fecha. Además hay 3 filas vacías al final. Tampoco está documentado qué significan los códigos `BAR`/`BAV` y `En`.
8. **Los conteos de `exploracion_datos.R` no salen del CSV.** Están escritos a mano (`20, 13, 2, 0, 7, 89, 58, 16, 115, 19, 2`) y no coinciden con ninguna suma del CSV, ni BAR, ni BAV, ni BA. Además, cruzan cada captura con el clima del mismo día de recolección, cuando debería usarse el período entre dos recolecciones.
9. `Datos TP FINAL/modelo_lin_temp_MACN_Ma.R` llama "Aeroparque" a la estación 87585, que es Buenos Aires Observatorio.

## 5. Estructuras de correlación

- **Entre días consecutivos de una misma serie:** se modela con `corAR1(form = ~ t)`. Como `t` cuenta días, se respetan los huecos. El AR(1) no alcanza en la mayoría de las series, pero cambiarlo casi no mueve los resultados (ver abajo).
- **Entre predictoras:** los rezagos de TMAX están muy correlacionados entre sí, y también con el seno y el coseno. Por eso se informa el efecto acumulado y no cada coeficiente por separado.
- **Entre salas medidas al mismo tiempo: no se tiene en cuenta.** Las salas medidas en el mismo período comparten el clima y factores que el modelo no ve, como el viento, el sol o la climatización del edificio. Pasa con Pv en 2019, Ict en 2022–2023 e In en 2025. Mientras cada sala se ajuste por separado no es un problema, pero para comparar salas hace falta un modelo conjunto con un efecto aleatorio de fecha o un modelo multivariado.
- **Trampas:** las capturas van a estar correlacionadas entre trampas de la misma sala, entre recolecciones sucesivas de la misma trampa y entre coleópteros y lepidópteros. Hace falta un modelo mixto para conteos, por ejemplo `glmmTMB(capturas ~ clima + (1 | sala/trampa), family = nbinom2)`.

### Diagnóstico de residuos: `R/06_diagnostico_residuos.R`

Con la ventana elegida en el 03, cada serie se ajusta con tres estructuras de error, AR(1), AR(2) y ARMA(1,1), y se comparan por AIC. Además, se prueba si queda autocorrelación en los residuos normalizados (Ljung-Box a 7 y 14 días, calculado sobre una grilla diaria para respetar los huecos) y, si gana otra estructura, se vuelve a elegir la ventana con ella.

| Serie | Mejor estructura | ΔAIC del AR(1) | Efecto AR(1) | Efecto con la mejor | Cambio del EE | K con AR(1) → con la mejor |
|---|---|---|---|---|---|---|
| MACN-Pv (sin subsector) | ARMA(1,1) | 53,5 | 0,27 | 0,27 | +9 % | 5 → 5 |
| BA-Pb | ARMA(1,1) | 38,3 | 0,25 | 0,27 | +15 % | 7 → 14 |
| MACN-Ict Colección | ARMA(1,1) | 25,0 | 0,39 | 0,39 | +12 % | 10 → 10 |
| MACN-Ma | ARMA(1,1) | 23,5 | 0,47 | 0,45 | +12 % | 10 → 7 |
| MACN-In Abajo | AR(2) | 18,4 | 0,22 | 0,19 | +38 % | 4 → 1 |
| MACN-In Arriba | AR(2) | 14,8 | 0,28 | 0,26 | +36 % | 5 → 2 |
| MACN-Ar | ARMA(1,1) | 12,6 | 0,60 | 0,64 | −1 % | 2 → 2 |
| MACN-Ict Laboratorio | ARMA(1,1) | 4,2 | 0,50 | 0,50 | 0 % | 11 → 11 |
| MACN-Pv Ameghino | ARMA(1,1) | 2,0 | 0,50 | 0,50 | −2 % | 13 → 13 |
| BA-v | AR(1) | 0 | 0,63 | — | — | 7 |
| MACN-He | AR(1) | 0 | 0,27 | — | — | 9 |
| MACN-Pv Anexo | AR(1) | 0 | 0,55 | — | — | 2 |

- **El AR(1) no alcanza en 9 de 12 series:** en casi todas ajusta mejor un ARMA(1,1). En MACN-Ict Colección el AR(2) no converge.
- **Las conclusiones no cambian.** El efecto acumulado se mueve como mucho 0,04 °C/°C y el EE crece entre 0 y 15 % (36–38 % en MACN-In, que tiene solo 52 días). Las ventanas largas de Ameghino, Ict Laboratorio e Ict Colección se mantienen. Cambian las ventanas de las series ya ambiguas: Ma pasa de 10 a 7 (estaba en el rango 7–11), BA-Pb de 7 a 14, y las de MACN-In se acortan mucho, pero con 52 días no son confiables.
- **Queda autocorrelación en BA-v y en MACN-Ict Colección** aun con la mejor estructura (Ljung-Box p < 0,001). En BA-v los valores son chicos (|ACF| ≤ 0,10) y significativos porque la serie tiene 1269 días; además mezcla archivos de años e intervalos de muestreo distintos.
- **Aparece un ciclo semanal:** la ACF de los residuos tiene un pico en el rezago 7 en MACN-Ar (0,42 con AR(1), 0,32 con ARMA(1,1)), MACN-Ict Colección (0,18) y, más débil, en BA-v (0,10). Sugiere un efecto del día de la semana (luces, gente, climatización en días hábiles), que el modelo no incluye.
- Figura: `resultados/figuras/fig_residuos_acf.png`. Tablas: `resultados/tablas/residuos_estructuras.csv` y `residuos_resumen.csv`.

### Día de la semana: `R/07_dia_semana.R`

Por el pico de la ACF a los 7 días, al modelo de rezagos (K del 03, AR(1)) se le agrega el día de la semana de dos formas: una indicadora de día hábil (lunes a viernes) o un factor de 7 niveles con el domingo como referencia. Los feriados todavía cuentan como días hábiles.

| Serie | Hábil − fin de semana (°C, EE) | ΔAIC hábil | ΔBIC hábil | ΔAIC factor | ACF a 7 días: sin días → con hábil |
|---|---|---|---|---|---|
| MACN-Ar | **+1,08** (0,20) | −25,5 | −23,1 | −23,5 | 0,42 → 0,17 |
| MACN-Pv Ameghino | **+0,22** (0,04) | −33,7 | −30,1 | −35,2 | −0,03 → −0,17 |
| MACN-Ict Colección | **+0,13** (0,02) | −31,8 | −27,2 | −28,8 | 0,17 → 0,13 |
| MACN-In Arriba | +0,15 (0,06) | −5,0 | −3,0 | −1,7 | 0,06 → −0,02 |
| MACN-Pv (sin subsector) | +0,12 (0,05) | −4,1 | −0,3 | −0,9 | 0,05 → 0,03 |
| MACN-In Abajo | +0,13 (0,06) | −3,3 | −1,4 | −0,3 | 0,10 → 0,04 |
| BA-Pb | +0,05 (0,03) | −1,8 | +2,0 | −0,9 | 0,06 → 0,05 |
| MACN-He | +0,04 (0,04) | +0,8 | +4,1 | −0,2 | 0,00 → −0,01 |
| MACN-Pv Anexo | +0,08 (0,10) | +1,4 | +5,3 | +8,4 | −0,02 → −0,02 |
| MACN-Ict Laboratorio | +0,02 (0,03) | +1,6 | +5,6 | +0,9 | −0,01 → −0,01 |
| MACN-Ma | +0,01 (0,03) | +1,9 | +5,7 | +4,2 | 0,02 → 0,02 |
| BA-v | 0,00 (0,02) | +2,0 | +7,1 | +5,6 | 0,11 → 0,11 |

- **Hay un efecto claro de los días hábiles en tres salas**, con AIC y con BIC: MACN-Ar (+1,1 °C), Pv Ameghino (+0,2 °C) e Ict Colección (+0,1 °C). En MACN-In y Pv sin subsector es chico y dudoso. En el resto no aparece.
- **Alcanza con la indicadora de día hábil.** El factor de 7 días solo mejora algo el AIC en Ameghino, y el BIC siempre prefiere la indicadora. El factor muestra, eso sí, que en varias salas el calor se acumula a lo largo de la semana: en Ameghino pasa de +0,21 °C el lunes a +0,36 °C el viernes, y el sábado queda en un valor intermedio (+0,14 °C).
- **No cambia el efecto de la TMAX exterior:** el efecto acumulado se mueve menos de 0,05 en todas las series, porque el día de la semana no está correlacionado con el clima.
- **El efecto cambia a lo largo del año** (diferencia media de los residuos entre días hábiles y fines de semana, por mes):
  - MACN-Ar: 0,3 °C en abril, 1,2 °C en mayo y 2,3 °C en junio de 2016. Crece con el frío, lo que sugiere **calefacción en días hábiles**.
  - Pv Ameghino: positivo en todos los meses (0,06 a 0,38 °C), así que es robusto.
  - Ict Colección: casi nulo en 2022–2023 y claro en 2024 (0,14 a 0,50 °C). Los dos períodos vienen de archivos distintos, separados por 11 meses.
  - BA-v: el promedio es cero porque el signo cambia. Entre octubre de 2016 y febrero de 2017, los días hábiles fueron entre 0,3 y 0,5 °C **más frescos**, lo que podría ser aire acondicionado en días hábiles en verano. Es una hipótesis: habría que confirmarla con el museo.
- El pico a los 7 días baja mucho en MACN-Ar (de 0,42 a 0,17) y algo en Ict Colección. En BA-v no cambia, coherente con un efecto que cambia de signo según la estación.
- Figura: `resultados/figuras/fig_dia_semana.png`. Tablas: `resultados/tablas/dia_semana_resumen.csv`, `dia_semana_modelos.csv` y `dia_semana_coeficientes.csv`.

## 6. Siguientes pasos a probar

**Arreglos de datos**

1. Separar los dataloggers de MACN-Or y revisar las series que mezclan archivos (punto 2 de los problemas de la sección 2 y puntos 1 a 4 de la sección 4). Lo más robusto es que cada archivo sea su propia serie, salvo que se confirme que se trata del mismo datalogger. Hay que tener en cuenta el intervalo de muestreo de cada serie, o al menos indicarlo, cuando se comparan máximas diarias.
2. Descartar los días incompletos: quedarse solo con los días que tienen todos los registros esperados. *(En `modelo_rezagos.R` ya se descartan el primer y el último día de cada archivo.)*
3. Hacer que el script avise cuando una serie falla, en vez de descartarla. *(Hecho en `modelo_rezagos.R`.)*
4. Exportar los `.hobo` que faltan, sobre todo el de BA-v de 12/2019 a 03/2020, y corregir `general_trampas.csv`.

**Modelo**

5. **Medir el aporte real del exterior:** comparar cada modelo contra uno sin `TMAX_prom3d` (ΔAIC, o R² parcial). *(Hecho en `modelo_rezagos.R`, columna `dAIC_vs_sin_ext`.)*
6. **Rehacer la comparación de ventanas** con `method = "ML"` y ventanas de más de 5 días. *(Hecho: `R/03_modelo_rezagos.R` prueba hasta 14 días con ML. Queda como alternativa un modelo de inercia con la máxima interior de ayer como predictora, o una media móvil exponencial, cuyo parámetro de decaimiento mide directamente la inercia de cada sala.)*
7. **Modelo conjunto de todas las salas**, para comparar pendientes con su incertidumbre:
   ```r
   gls(temp_max ~ TMAX_prom3d * serie + sin_year + cos_year,
       correlation = corAR1(form = ~ t | serie),
       weights = varIdent(form = ~ 1 | serie))
   ```
   Otra opción es `lme` con pendiente aleatoria por serie.
8. **Comparar salas medidas en el mismo período**, para no confundir sala con año: 2019 (los 3 registros de MACN-Pv y BA-v), 2020 (MACN-Ma, MACN-Or y BA-v) y 2022–2023 (MACN-Ict colección y laboratorio, y BA-Pb).
9. **Otras respuestas:** temperatura media y mínima, amplitud diaria y, sobre todo, la **humedad relativa interior** contra la humedad exterior, que es la variable clave para la conservación.
10. **Residuos:** *(Hecho en `R/06_diagnostico_residuos.R`: ARMA(1,1) ajusta mejor en la mayoría, pero los resultados no cambian.)* *(Día de la semana: hecho en `R/07_dia_semana.R`.)* Faltan los feriados y una interacción entre día hábil y estación, para medir si la climatización calienta en invierno y enfría en verano. Lo ideal es preguntarle al museo cuándo se usa la climatización de cada sala.
11. **Robustez de la ventana con BIC.** *(Hecho en `R/03_modelo_rezagos.R`: la velocidad se mantiene, pero la cola larga de las aisladas depende del criterio. Ver la sección 3.)* El siguiente paso es un modelo con forma impuesta para los rezagos.

**Trampas**

12. `general_trampas.csv` ya tiene las capturas por trampa y fecha (BA, Ma, En y Or). El paso natural es relacionarlas con el microclima de la misma sala en los días previos a cada recolección, por ejemplo con un modelo de Poisson o binomial negativo para los conteos.
