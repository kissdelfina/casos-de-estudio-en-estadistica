# Ideas para la presentación oral

Notas de trabajo para armar la presentación. Las pautas (`docs/pautas_oral.pdf`) piden:

- 15 minutos de exposición más 5 de preguntas;
- entre 10 y 12 diapositivas, cada una con una sola idea;
- títulos de 32 pt o más y texto y etiquetas de 20 a 24 pt o más;
- gráficos con título informativo, ejes con unidades y color con un propósito;
- explicar los supuestos y la justificación de cada modelo, e interpretar cada número en contexto.

## Propuesta de diapositivas

| # | Sección | Mensaje clave (el título de la diapositiva) | Material | Estado |
|---|---|---|---|---|
| 1 | Carátula | Microclima de las salas del MACN y su relación con el clima exterior | — | Falta |
| 2 | Contexto | El clima exterior entra a las salas, pero ¿cuánto y con qué retraso? | Foto o esquema: SMN → edificio → sala → colección | Falta |
| 3 | Datos | 13 dataloggers, 2014–2025, contra la estación SMN Observatorio | Línea de tiempo de los períodos medidos por serie | Falta: se puede sacar de la tabla "Datos usados" del documento de rezagos |
| 4 | Limpieza | Cuatro decisiones de limpieza que cambian los resultados | Lista corta: `max()` sobre texto, días incompletos, dataloggers mezclados (Or), intervalos de muestreo distintos | Falta |
| 5 | Modelo | GLS con AR(1): los días consecutivos no son independientes | La ecuación con rezagos + una línea sobre φ | Falta |
| 6 | Por qué el AR(1) | Sin AR(1) todo parece significativo | AIC con y sin AR(1) en MACN-Ar (309 contra 270), o la ACF de los residuos | Falta la ACF |
| 7 | Resultado 1 | Hay salas expuestas y salas aisladas | Gráfico del efecto acumulado contra la pendiente del promedio (documento de rezagos) | Hecho, falta exportarlo con letra grande |
| 8 | Resultado 2 | Las expuestas reaccionan el mismo día; las aisladas, al día siguiente y de a poco | Perfil de coeficientes por rezago (documento de rezagos). Para la diapositiva conviene mostrar 3 o 4 salas, no las 12 | Hecho, falta una versión reducida |
| 9 | Resultado 3 | En la mayoría de las salas se puede pronosticar con información del día anterior | `resultados/figuras/fig_pronostico_rmse.png` (`R/05_modelo_pronostico.R`) | Hecho |
| 10 | Trampas | Más calor, más capturas: +17 % por °C | Gráfico de la tasa diaria contra la temperatura (`fig_trampas_temperatura.png`) | Existe; revisar las limitaciones (ver abajo) |
| 11 | Conclusiones | Qué salas priorizar y qué medir | 3 viñetas: salas expuestas a aislar, salas aisladas estables, recomendaciones de medición | Falta |
| 12 | Limitaciones | Lo que este análisis no puede decir | Series cortas, dataloggers mezclados, correlación entre salas sin modelar, solo temperatura | Falta |

## Números que conviene tener a mano

- Salas expuestas (BA-v, MACN-Ar, MACN-Pv Anexo): efecto acumulado de 0,55 a 0,63 °C por cada °C exterior.
- Salas aisladas: de 0,22 a 0,50, es decir, entre 2 y 4 veces la pendiente del promedio de 3 días. Varias (Pv Ameghino, Ict Laboratorio, Ma) reciben casi tanto como las expuestas, pero más lento.
- Ventanas elegidas (probando de 0 a 14 días): 2 días en Pv Anexo y MACN-Ar, 7 en BA-v, de 9 a 13 en las aisladas (Ameghino 13). Ninguna elige el máximo, así que 14 alcanza.
- Las expuestas reciben el 86–93 % del efecto en el mismo día y el siguiente; las aisladas, solo el 24–42 %.
- φ del AR(1): 0,4–0,5 en las salas expuestas y de 0,8 a 0,93 en las aisladas.
- Residuos: un ARMA(1,1) ajusta mejor que el AR(1) en 9 de 12 series, pero el efecto acumulado cambia como mucho 0,04 y el EE crece hasta 15 %.
- Pronóstico sin la TMAX de hoy (11 series con 60 días evaluados o más):
  - el error sube mucho solo en BA-v (+66 %) y Pv Anexo (+55 %), las salas rápidas;
  - le siguen Ict Laboratorio (+24 %), Pv Ameghino (+17 %) e Ict Colección (+17 %); en las otras 6 sube 12 % o menos, o incluso baja;
  - con máximo 7 días, Ameghino daba +49 %: era un artefacto de la ventana corta;
  - en todas, la TMAX de los días anteriores mejora entre 12 % y 40 % respecto de usar solo la época del año;
  - el error típico sin la TMAX de hoy va de 0,5 °C (Or Azara) a 1,6 °C (Pv Anexo).
- Trampas: +17 % de capturas por cada °C de temperatura media exterior (p = 0,017), con dispersión 21,7 y 46 observaciones.

## Preguntas probables en la ronda final

- **¿Por qué GLS y no una regresión lineal común?** Porque los residuos de días consecutivos están correlacionados (φ de hasta 0,93). Si se ignora, los errores estándar quedan demasiado chicos.
- **¿Por qué AIC y no BIC para elegir la ventana?** Con BIC las aisladas eligen 5–7 días en vez de 9–13, y su efecto acumulado baja (Ameghino de 0,50 a 0,29). La diferencia de velocidad entre salas se mantiene con los dos. Preferimos AIC porque buscamos predecir, y fuera de la muestra las ventanas largas pronosticaron mejor en 8 de 10 casos. Conviene presentar el efecto de las aisladas como un rango.
- **¿Revisaron que el AR(1) alcance?** Sí (`R/06_diagnostico_residuos.R`): en 9 de 12 series ajusta mejor un ARMA(1,1), pero el efecto acumulado y las ventanas largas de las salas aisladas no cambian. Queda un pico de autocorrelación a los 7 días en algunas salas, que sugiere un efecto del día de la semana.
- **¿Por qué no usaron el promedio de 3 días?** Fuerza a que los 3 días pesen igual y a que el efecto termine a los 3 días. El test de razón de verosimilitud rechaza los pesos iguales en 8 de 12 series, y en las salas aisladas el promedio subestima el efecto entre 2 y 4 veces.
- **¿Los coeficientes de cada rezago no son inestables por la colinealidad?** Sí (r ≈ 0,85 entre días consecutivos). Por eso informamos la suma, que tiene un error estándar chico porque las covarianzas negativas se compensan.
- **¿El R² de 0,95 significa que el exterior explica el 95 %?** No: lo domina el ciclo anual. Para medir el aporte del exterior comparamos contra un modelo sin datos del exterior (ΔAIC).
- **¿Se pueden comparar las salas entre sí?** Solo de forma descriptiva: se midieron en años distintos y se ajustaron por separado. Para compararlas formalmente haría falta un modelo conjunto.
- **¿Cómo evaluaron el pronóstico?** Fuera de la muestra, con origen móvil: se entrena con la primera mitad de la serie, se predicen los 30 días siguientes, se reentrena y se repite. La ventana de rezagos se elige solo con el primer período de entrenamiento, para no mirar el futuro.
- **¿Por qué no se incluyen MACN-Ar ni MACN-In en el pronóstico?** Tienen 78 y 52 días: el primer entrenamiento queda de unos 40 días, el seno y el coseno extrapolan mal, y el error no es comparable.
- **¿Por qué el modelo de trampas usa la temperatura exterior y no la de la sala?** Hay que tener la respuesta preparada, o cambiarlo antes de presentar.

## Cosas para resolver antes de presentar

- El informe dice "7 días" de inercia para algunas salas, pero la comparación de ventanas de `R/02_exploracion_datos.R` llega solo hasta 5.
- El informe usa la pendiente del promedio de 3 días, que subestima el efecto en las salas aisladas. Decidir si se pasa al efecto acumulado de los rezagos.
- El resumen del informe dice "14 sectores", pero la tabla muestra 9.
- Las figuras actuales usan `base_size` de 10 a 14: para las diapositivas hay que exportarlas con letra más grande (unos 20 pt).
