# Resultados y siguientes pasos

Estado al 03/11/2026, después de revisar el commit `64969b7` ("modelo_temperatura tiene gls con tiempo defasado").

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
- `Datos TP FINAL/plot_zoom_png?width=1470&height=923` tiene un `?` en el nombre. En Windows impide hacer `git pull`, así que hay que borrarlo o renombrarlo.
- El xlsx del SMN 1991–2020 quedó commiteado (72 MB). Funciona, pero cualquier cambio a ese archivo va a hacer crecer el repo.

## 3. Siguientes pasos a probar

**Arreglos de datos**

1. Separar los dataloggers de MACN-Or y revisar las series que mezclan archivos (punto 2 de los problemas). Lo más robusto es que cada archivo sea su propia serie, salvo que se confirme que se trata del mismo datalogger.
2. Descartar los días incompletos: quedarse solo con los días que tienen todos los registros esperados.
3. Hacer que el script avise cuando una serie falla, en vez de descartarla.

**Modelo**

4. **Medir el aporte real del exterior:** comparar cada modelo contra uno sin `TMAX_prom3d` (ΔAIC, o R² parcial).
5. **Rehacer la comparación de ventanas** con `method = "ML"` y ventanas de más de 5 días (por ejemplo, 7, 10 y 14). Otra opción es una media móvil exponencial, que no obliga a elegir un número de días y cuyo parámetro de decaimiento mide directamente la inercia de cada sala.
6. **Modelo conjunto de todas las salas**, para comparar pendientes con su incertidumbre:
   ```r
   gls(temp_max ~ TMAX_prom3d * serie + sin_year + cos_year,
       correlation = corAR1(form = ~ t | serie),
       weights = varIdent(form = ~ 1 | serie))
   ```
   Otra opción es `lme` con pendiente aleatoria por serie.
7. **Comparar salas medidas en el mismo período**, para no confundir sala con año: 2019 (los 3 registros de MACN-Pv y BA-v), 2020 (MACN-Ma, MACN-Or y BA-v) y 2022–2023 (MACN-Ict colección y laboratorio, y BA-Pb).
8. **Otras respuestas:** temperatura media y mínima, amplitud diaria y, sobre todo, la **humedad relativa interior** contra la humedad exterior, que es la variable clave para la conservación.
9. **Residuos:** con φ ≈ 0,9 conviene probar un ARMA(2,1) o similar, y revisar la ACF de los residuos normalizados.

**Trampas**

10. `general_trampas.csv` ya tiene las capturas por trampa y fecha (BA, Ma, En y Or). El paso natural es relacionarlas con el microclima de la misma sala en los días previos a cada recolección, por ejemplo con un modelo de Poisson o binomial negativo para los conteos.
