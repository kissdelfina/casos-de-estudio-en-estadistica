# Primer modelo: MACN-Ar (`modelo_lin_temp_MACN_Ar.R`)

## Enfoque

- **Respuesta:** temperatura máxima diaria interior de la sala MACN-Ar (HOBO, 05/04/2016 – 23/06/2016).
- **Predictoras:** temperatura máxima (`TMAX`) y humedad relativa (`HUM_REL`) exteriores del SMN, estación Buenos Aires Observatorio, y su interacción.
- **Modelo:** GLS con errores AR(1) sobre un índice entero de días, para que los días faltantes se respeten:

```r
gls(temp_max ~ TMAX * HUM_REL, correlation = corAR1(form = ~ t))
```

- **Alternativa anotada:** un GAM con spline cíclico sobre el mes, para la estacionalidad.

## Problemas

- **Datos exteriores fuera del repo:** los de 2016 vienen de `smn_datos_meteorologicos_1991_2020.xlsx` (72 MB, ignorado por git). Cada uno lo tiene que descargar del SMN y ponerlo en `Datos TP FINAL/`. Cubre los 80 días del período con `TMAX` y `HUM_REL` completos.
- **Rutas absolutas:** el script lee los datos desde `D:/...`. Hay que cambiarlas por rutas relativas al repo.
- **Un solo período:** son unos 80 días de otoño, así que no se puede estimar la estacionalidad.
- **`max()` sobre texto:** la temperatura se lee como texto por la coma decimal y `max()` la compara alfabéticamente. Hay que convertirla a número antes de resumir, o leer con `read.csv2`.
- **Datos faltantes:** si falta un día del SMN, `gls` falla (`na.fail`).
- **Sin rezagos:** el edificio puede tardar en responder al clima exterior. Conviene probar predictoras rezagadas.
- **Sin diagnóstico:** no se miran los residuos ni su ACF, ni se compara contra un modelo sin AR(1).
