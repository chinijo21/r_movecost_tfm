# =====================================================================
#  movecost para Can Ribes II  —  analisis completo (entrada + salida)
# 
#
#  Hace en una sola ejecucion:
#    - ENTRADA: captacion de recursos (isocronas + ha de recurso por banda
#               + tiempo minimo de acceso a cada area de recurso)
#    - SALIDA : rutas de menor coste (LCP) a los destinos, con km y tiempo
#    - ROBUSTEZ: repite con varias funciones de coste (Tobler, Pandolf,
#               Herzog, Uriarte, Llobera, vehiculo con ruedas)
#    - EXPORTA: rasters de coste, isocronas, LCP (shp) y tablas (csv)
#
#  Los PLACEHOLDERS estan marcados con  ###<-- EDITAR.
#  El script SALTA los recursos cuyo shapefile aun no exista, asi que
#  puedes ejecutarlo ya con lo que tengas e ir anadiendo poligonos.
# =====================================================================
# ---------------------------------------------------------------------
# install.packages(c("terra","raster","sp","sf","gdistance","igraph"))
# install.packages("movecost")
# print(packageVersion("movecost"))   # debe ser 2.x

# ---------------------------------------------------------------------
#  PASO 1 — LIBRERIAS, RUTAS Y CARGA
# ---------------------------------------------------------------------
library(movecost); library(raster); library(sp); library(sf)

setwd("RUTA/A/TU/CARPETA")            ###<-- EDITAR: carpeta con tus datos
crs_proyecto <- 25831                  # ETRS89 / UTM 31N (el de tu MDT)

ruta_mdt     <- "recortado.tif"        ###<-- EDITAR si cambia el nombre
ruta_origen  <- "can_ribes.shp"        ###<-- EDITAR: punto de Can Ribes II
ruta_destinos<- "destinos.shp"         ###<-- EDITAR: puntos Ebusus/EsCodolar/PlatjaBossa
                                       #     (debe tener un campo "name")
ruta_barrera <- NA                     ###<-- EDITAR (opcional): poligono salinas/costa
                                       #     deja NA si aun no lo tienes

# Poligonos de RECURSO (areas). Anade/comenta segun los vayas teniendo:
recursos <- c(
  pinar    = "pinar.shp",              ###<-- EDITAR (poligono)
  matorral = "matorral.shp",           ###<-- EDITAR (poligono)
  margas   = "margas_U7.shp"           ###<-- EDITAR (poligono)
)

# Parametros del analisis
funct_principal <- "t"                 # "t"=Tobler (peatonal). Ver codigos abajo.
paso_isocrona   <- 0.5                 # ancho de banda en horas (0.5 = 30 min)
dir.create("movecost_resultados", showWarnings = FALSE)

# --- Carga y blindaje (esto evita los cuelgues) ---
dtm <- raster(ruta_mdt)
cat("MDT:", ncell(dtm), "celdas | resolucion:", res(dtm)[1], "m\n")

cargar_pts <- function(ruta){
  p <- as(st_transform(st_read(ruta, quiet=TRUE), crs_proyecto), "Spatial")
  crs(p) <- crs(dtm); p
}
origen   <- cargar_pts(ruta_origen)
destinos <- cargar_pts(ruta_destinos)
if(is.null(destinos$name)) destinos$name <- paste0("destino_", seq_len(length(destinos)))

# Barrera (opcional): se pasa al argumento 'barrier' de movecost.
# IMPORTANTE: NO uses NA en el MDT para el mar; rompe la generacion de
# isocronas. Usa SIEMPRE un poligono via 'barrier'.
barrera <- NULL
if(!is.na(ruta_barrera)){
  barrera <- as(st_transform(st_read(ruta_barrera, quiet=TRUE), crs_proyecto), "Spatial")
  crs(barrera) <- crs(dtm)
}

# Comprobacion: que los puntos caigan sobre celdas validas del MDT
stopifnot(all(!is.na(extract(dtm, origen))),
          all(!is.na(extract(dtm, destinos))))
cat("Puntos OK dentro del MDT.\n")

# ---------------------------------------------------------------------
#  PASO 2 — LOGISTICA DE ENTRADA: captacion de recursos (tu §5)
#  Superficie de coste desde Can Ribes + isocronas + estadistica por area
# ---------------------------------------------------------------------
cat("\n>> ENTRADA: superficie de coste y captacion de recursos\n")
ent  <- movecost(dtm = dtm, origin = origen, funct = funct_principal,
                 time = "h", barrier = barrera, graph.out = FALSE)
cost <- ent$accumulated.cost.raster
writeRaster(cost, "movecost_resultados/coste_acumulado_h.tif", overwrite = TRUE)

# Isocronas propias por reclasificacion (robusto, control total)
topo   <- ceiling(cellStats(cost, max) / paso_isocrona) * paso_isocrona
cortes <- seq(0, topo, by = paso_isocrona)
rcl    <- cbind(cortes[-length(cortes)], cortes[-1], seq_len(length(cortes)-1))
bandas <- reclassify(cost, rcl)
writeRaster(bandas, "movecost_resultados/isocronas_banda.tif", overwrite = TRUE)
cell_ha <- prod(res(cost)) / 10000     # 25 m -> 0.0625 ha/celda

# Para cada recurso (area): tiempo minimo de acceso + ha por banda
tabla_recursos <- data.frame()
for(nm in names(recursos)){
  f <- recursos[[nm]]
  if(!file.exists(f)){ cat("  [salto]", nm, "- falta", f, "\n"); next }
  poly <- as(st_transform(st_read(f, quiet=TRUE), crs_proyecto), "Spatial"); crs(poly) <- crs(dtm)
  v_cost  <- unlist(raster::extract(cost,   poly))
  v_banda <- unlist(raster::extract(bandas, poly))
  tmin    <- round(min(v_cost, na.rm = TRUE), 2)
  por_banda <- round(table(v_banda) * cell_ha, 1)
  cat(sprintf("  %s: acceso minimo = %.2f h | ha por banda: %s\n",
              nm, tmin, paste(names(por_banda), por_banda, sep="=", collapse=", ")))
  tabla_recursos <- rbind(tabla_recursos,
    data.frame(recurso = nm, acceso_min_h = tmin,
               banda = names(por_banda), ha = as.numeric(por_banda)))
}
if(nrow(tabla_recursos)) write.csv(tabla_recursos,
   "movecost_resultados/captacion_recursos.csv", row.names = FALSE)

# ---------------------------------------------------------------------
#  PASO 3 — LOGISTICA DE SALIDA: rutas de exportacion (tu §6)
#  LCP de Can Ribes a cada destino TERRESTRE (distancia, coste, tiempo)
#  Nota: el tramo MARITIMO (cabotaje) NO lo hace movecost; se suma aparte
#        como distancia / velocidad de navegacion (ver tu §4.5).
# ---------------------------------------------------------------------
cat("\n>> SALIDA: rutas de menor coste (LCP)\n")
sal <- movecost(dtm = dtm, origin = origen, destin = destinos,
                funct = funct_principal, time = "h",
                barrier = barrera, graph.out = TRUE, export = TRUE)
salida <- sal$dest.loc.w.cost@data
salida$dist_km <- round(sal$LCPs@data$length / 1000, 2)
salida$tiempo_h <- round(salida$cost, 2)
print(salida[, c("name","dist_km","tiempo_h","cost_hms")])
write.csv(salida[, c("name","dist_km","tiempo_h","cost_hms")],
          "movecost_resultados/rutas_salida.csv", row.names = FALSE)
shapefile(sal$LCPs, "movecost_resultados/LCP_salida.shp", overwrite = TRUE)

# ---------------------------------------------------------------------
#  PASO 4 — CONTRASTE DE ROBUSTEZ (responde a Agustin p.35 y p.40)
#  Repite los LCP con varias funciones de coste y compara.
#  Codigos verificados (movecost 2.2):
#    "t"  Tobler         "p"   Pandolf 1977 (p.40)   "hrz" Herzog
#    "ug" Uriarte (p.35) "ls"  Llobera (p.35)        "wcs" vehiculo ruedas
#         (wcs conecta con el umbral 8-10% de Raepsaet, util para el carro)
# ---------------------------------------------------------------------
cat("\n>> ROBUSTEZ: varias funciones de coste\n")
funciones <- c("t","p","hrz","ug","ls","wcs")
robustez <- data.frame()
for(f in funciones){
  r <- tryCatch(movecost(dtm=dtm, origin=origen, destin=destinos,
                         funct=f, time="h", barrier=barrera, graph.out=FALSE),
                error = function(x){ cat("  [error]", f, "\n"); NULL })
  if(is.null(r)) next
  robustez <- rbind(robustez, data.frame(
    funcion = f, destino = r$dest.loc.w.cost@data$name,
    tiempo_h = round(r$dest.loc.w.cost@data$cost, 2),
    dist_km  = round(r$LCPs@data$length/1000, 2)))
}
print(robustez)
write.csv(robustez, "movecost_resultados/robustez_funciones.csv", row.names = FALSE)

cat("\n=== LISTO ===\n",
    "Resultados en la carpeta 'movecost_resultados':\n",
    " - coste_acumulado_h.tif / isocronas_banda.tif  (rasters para QGIS)\n",
    " - captacion_recursos.csv     (tu §5: ha de recurso por isocrona)\n",
    " - rutas_salida.csv + LCP_salida.shp  (tu §6: Tablas 3-4 y Mapa 6)\n",
    " - robustez_funciones.csv     (contraste de algoritmos, Agustin p.35/p.40)\n", sep="")

# ---------------------------------------------------------------------
#  NOTAS
#  1. El cabotaje (tramo marino) queda fuera de movecost: suma el tiempo
#     marino como distancia/velocidad de navegacion y compon el bimodal.
#  2. El escenario "asno" no tiene funcion calibrada: usa "wcs" o ajusta
#     parametros y declaralo como modelo de CONTRASTE, no empirico.
#  3. Para el mar usa barrera (poligono), nunca NA en el MDT.
# ---------------------------------------------------------------------
