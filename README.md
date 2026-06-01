```markdown
# Modelado de Coste de Movimiento Terrestre y Marítimo - Can Ribes II

Este repositorio contiene el script principal en R desarrollado para el análisis de paleoconectividad, cálculo de superficies de coste acumulado y modelado bimodal (terrestre-marítimo) aplicado al yacimiento arqueológico de Can Ribes II. 

El análisis integra variables topográficas, distribución de recursos bióticos y geológicos, y penalizaciones físicas ambientales (viento y oleaje) para evaluar los patrones de movilidad y captación económica.

## Estructura del Análisis

El script (`movecost_canribes.R`) centraliza el flujo de trabajo en cinco bloques secuenciales:

1. **Preparación y Blindaje de Datos:** Carga del Modelo Digital del Terreno (MDT), definición del CRS del proyecto (ETRS89 / UTM 31N) y procesamiento de geometrías (conversión automática de multipuntos a puntos y control de valores NoData).
2. **Logística de Entrada (Áreas de Captación):** Generación de la superficie de coste acumulado en horas desde el yacimiento, cálculo de bandas isócronas personalizadas y extracción automatizada de estadísticas de recursos (hectáreas disponibles por recurso físico según intervalos de tiempo).
3. **Logística de Salida (Rutas de Menor Coste):** Trazado de los *Least Cost Paths* (LCP) terrestres desde el origen hasta los destinos especificados, exportando métricas de distancia (km) y tiempos estimados de trayecto.
4. **Análisis de Robustez:** Evaluación comparativa del comportamiento de las rutas utilizando seis funciones de coste algorítmicas diferentes (Tobler, Pandolf, Herzog, Uriarte, Llobera y Wheeler/wcs para transporte rodado) con el fin de contrastar los modelos analíticos.
5. **Modelo Marítimo Bimodal (Cabotaje):** Implementación de un modelo de fricción anisotrópica para la navegación costera basado en datos de reanálisis climático (viento y oleaje de archivos NetCDF). Incorpora límites operativos de la embarcación (umbral de oleaje de 1.5 m y vientos superiores a 15 nudos) y restricciones de visibilidad costera (fórmula de Andrew 2020).

## Requisitos del Sistema

Para ejecutar el script es necesario contar con un entorno de R con las siguientes librerías instaladas:

```R
install.packages(c("movecost", "raster", "sp", "sf", "gdistance", "ncdf4", "Matrix", "igraph
