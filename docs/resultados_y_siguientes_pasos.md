# Resultados y siguientes pasos

Estado al 11/11/2026. Las secciones 1 y 2 vienen de la revisión del commit `64969b7` ("modelo_temperatura tiene gls con tiempo defasado"). Las secciones 3 a 5 se agregaron después: el modelo con rezagos separados (`modelo_rezagos.R`), problemas nuevos en los datos y estructuras de correlación.

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

## 3. Rezagos separados: `modelo_rezagos.R`

En vez del promedio de 3 días, cada día anterior entra como predictora separada:

```r
gls(temp_max ~ TMAX_l0 + TMAX_l1 + ... + TMAX_lK + sin_year + cos_year,
    correlation = corAR1(form = ~ t), method = "ML")
```

`TMAX_lk` es la máxima exterior de hace k días. Para cada serie se prueba K = 0 a 7 y se elige el de menor AIC. Todos los modelos se ajustan por ML y sobre las mismas filas, así que los AIC son comparables. Además, el script:

- calcula los rezagos sobre la serie completa del SMN, así que no pierde los primeros días de cada serie;
- descarta el primer y el último día de cada archivo, que están incompletos;
- avisa cuando una serie falla, en vez de descartarla en silencio;
- usa rutas relativas y `read.csv2`.

**Efecto acumulado:** es la suma de los coeficientes, β₀ + … + β_K. Indica cuánto sube la máxima interior si la exterior sube 1 °C y se mantiene así K+1 días seguidos. Equivale también al efecto total, sumado sobre los días siguientes, de un solo día 1 °C más caluroso. Se puede comparar directamente con la pendiente del promedio de 3 días, que es el caso particular β₀ = β₁ = β₂. Los rezagos están muy correlacionados entre sí (r ≈ 0,85 entre días consecutivos), así que cada β por separado es inestable, pero la suma no: su error estándar es √(1ᵀV1).

| Serie | Rezagos elegidos | Efecto acumulado (EE) | Pendiente del promedio de 3 días | ΔAIC contra el promedio |
|---|---|---|---|---|
| BA-v | 0 a **7** | 0,63 (0,02) | 0,40 | −657 |
| MACN-Ar | 0 a 2 | 0,60 (0,10) | 0,61 | −5 |
| MACN-Pv Anexo | 0 a 2 | 0,55 (0,02) | 0,54 | −152 |
| MACN-Ma | 0 a **7** | 0,40 (0,03) | 0,20 | −92 |
| MACN-Ict Laboratorio | 0 a **7** | 0,40 (0,03) | 0,22 | −113 |
| MACN-Pv Ameghino | 0 a **7** | 0,32 (0,03) | 0,13 | −44 |
| MACN-Ict Colección | 0 a **7** | 0,31 (0,02) | 0,11 | −99 |
| MACN-In Arriba | 0 a 5 | 0,28 (0,05) | 0,14 | −7 |
| MACN-Pv (sin subsector) | 0 a 5 | 0,27 (0,03) | 0,17 | −23 |
| BA-Pb | 0 a **7** | 0,25 (0,03) | 0,09 | −74 |
| MACN-He | 0 a **7** | 0,23 (0,03) | 0,10 | −24 |
| MACN-In Abajo | 0 a 4 | 0,22 (0,04) | 0,13 | −8 |
| MACN-Or | se saltea por las fechas repetidas | — | — | — |

- **Los rezagos separados ajustan mejor que el promedio en todas las salas.** En casi todas, el test de razón de verosimilitud rechaza los pesos iguales (p < 0,05). Las excepciones son MACN-He, Pv Ameghino, Pv sin subsector e In Arriba, con p entre 0,07 y 0,10.
- **El promedio de 3 días subestima el efecto en las salas aisladas:** el efecto acumulado es el doble o el triple de la pendiente del promedio. La diferencia entre salas expuestas y aisladas es menor de lo que parecía: de 0,6 contra 0,1–0,2 pasa a 0,6 contra 0,2–0,4.
- **Las salas responden con perfiles distintos.** Las expuestas (MACN-Ar, Pv Anexo) reaccionan sobre todo el mismo día, y el efecto se termina a los 2 días. BA-v también tiene su mayor peso el mismo día, pero el efecto sigue presente hasta 7 días después. En las aisladas, el mayor peso está en el **día anterior** y baja de a poco a lo largo de una semana.
- **En 8 de 12 series ganó el máximo probado (7 días).** La inercia real puede ser mayor: hay que probar ventanas más largas o rezagos con forma impuesta (polinomial o exponencial, por ejemplo con `dlnm`).
- El efecto acumulado es condicional al ciclo anual: mide la respuesta a días más calurosos o más fríos de lo normal para la época, no todo el vínculo entre interior y exterior.
- Al correr, R muestra más de 50 warnings que todavía no se revisaron. Probablemente vienen de leer el xlsx del SMN.

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

- **Entre días consecutivos de una misma serie:** se modela con `corAR1(form = ~ t)`. Como `t` cuenta días, se respetan los huecos. Falta verificar que el AR(1) alcance (ver el punto 10 de los siguientes pasos).
- **Entre predictoras:** los rezagos de TMAX están muy correlacionados entre sí, y también con el seno y el coseno. Por eso se informa el efecto acumulado y no cada coeficiente por separado.
- **Entre salas medidas al mismo tiempo: no se tiene en cuenta.** Las salas medidas en el mismo período comparten el clima y factores que el modelo no ve, como el viento, el sol o la climatización del edificio. Pasa con Pv en 2019, Ict en 2022–2023 e In en 2025. Mientras cada sala se ajuste por separado no es un problema, pero para comparar salas hace falta un modelo conjunto con un efecto aleatorio de fecha o un modelo multivariado.
- **Trampas:** las capturas van a estar correlacionadas entre trampas de la misma sala, entre recolecciones sucesivas de la misma trampa y entre coleópteros y lepidópteros. Hace falta un modelo mixto para conteos, por ejemplo `glmmTMB(capturas ~ clima + (1 | sala/trampa), family = nbinom2)`.

## 6. Siguientes pasos a probar

**Arreglos de datos**

1. Separar los dataloggers de MACN-Or y revisar las series que mezclan archivos (punto 2 de los problemas de la sección 2 y puntos 1 a 4 de la sección 4). Lo más robusto es que cada archivo sea su propia serie, salvo que se confirme que se trata del mismo datalogger. Hay que tener en cuenta el intervalo de muestreo de cada serie, o al menos indicarlo, cuando se comparan máximas diarias.
2. Descartar los días incompletos: quedarse solo con los días que tienen todos los registros esperados. *(En `modelo_rezagos.R` ya se descartan el primer y el último día de cada archivo.)*
3. Hacer que el script avise cuando una serie falla, en vez de descartarla. *(Hecho en `modelo_rezagos.R`.)*
4. Exportar los `.hobo` que faltan, sobre todo el de BA-v de 12/2019 a 03/2020, y corregir `general_trampas.csv`.

**Modelo**

5. **Medir el aporte real del exterior:** comparar cada modelo contra uno sin `TMAX_prom3d` (ΔAIC, o R² parcial). *(Hecho en `modelo_rezagos.R`, columna `dAIC_vs_sin_ext`.)*
6. **Rehacer la comparación de ventanas** con `method = "ML"` y ventanas de más de 5 días. *(`modelo_rezagos.R` ya prueba hasta 7 días con ML. Falta probar ventanas de 10 y 14 días, o una media móvil exponencial, cuyo parámetro de decaimiento mide directamente la inercia de cada sala.)*
7. **Modelo conjunto de todas las salas**, para comparar pendientes con su incertidumbre:
   ```r
   gls(temp_max ~ TMAX_prom3d * serie + sin_year + cos_year,
       correlation = corAR1(form = ~ t | serie),
       weights = varIdent(form = ~ 1 | serie))
   ```
   Otra opción es `lme` con pendiente aleatoria por serie.
8. **Comparar salas medidas en el mismo período**, para no confundir sala con año: 2019 (los 3 registros de MACN-Pv y BA-v), 2020 (MACN-Ma, MACN-Or y BA-v) y 2022–2023 (MACN-Ict colección y laboratorio, y BA-Pb).
9. **Otras respuestas:** temperatura media y mínima, amplitud diaria y, sobre todo, la **humedad relativa interior** contra la humedad exterior, que es la variable clave para la conservación.
10. **Residuos:** con φ ≈ 0,9 conviene probar un ARMA(2,1) o similar, y revisar la ACF de los residuos normalizados.

**Trampas**

11. `general_trampas.csv` ya tiene las capturas por trampa y fecha (BA, Ma, En y Or). El paso natural es relacionarlas con el microclima de la misma sala en los días previos a cada recolección, por ejemplo con un modelo de Poisson o binomial negativo para los conteos.
