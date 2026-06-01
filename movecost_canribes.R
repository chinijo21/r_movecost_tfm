Aquí tienes el código con los comentarios reescritos en primera persona, manteniendo intacta la lógica, estructura y formato del código original:

```R
# =====================================================================
#   - ENTRADA: capto mis recursos (isocronas, mis hectáreas por banda
#              y el tiempo mínimo de acceso a mis áreas de recurso)
#   - SALIDA : calculo mis rutas de menor coste (LCP) a mis destinos, con km y tiempo
#   - ROBUSTEZ: repito mi proceso con varias de mis funciones de coste (Tobler, 
#               Pandolf, Herzog, Uriarte, Llobera, vehículo con ruedas)
#   - EXPORTA: guardo mis rasters de coste, isocronas, LCP (shp) y tablas (csv)
#
# =====================================================================
# ---------------------------------------------------------------------
#install.packages(c("terra","raster","sp","sf","gdistance","igraph"))
#install.packages("movecost")
#print(packageVersion("movecost"))   # compruebo que esté usando mi versión 2.x

# ---------------------------------------------------------------------
#  PASO 1 — AQUÍ CARGO MIS LIBRERÍAS, RUTAS Y DATOS
# ---------------------------------------------------------------------
library(movecost); library(raster); library(sp); library(sf)

setwd("/home/chinijo/Documentos/Claude/Projects/TFM/capas_movecost")
crs_proyecto <- 25831                  # este es el sistema de coordenadas de mi proyecto

ruta_mdt     <- "mdt_unido.tif"        # mi MDT recortado y unido a 25 m
ruta_origen  <- "can_ribes.shp"        # mi punto de origen en el yacimiento
ruta_destinos<- "destinos.shp"         # mis puntos de destino
ruta_barrera <- "barrera_mar.shp"      # mi polígono para establecer la barrera del mar

# Mis polígonos de RECURSO. Los iré añadiendo o comentando según los vaya teniendo:
recursos <- c(
  pinar    = "pinar.shp",              # mi capa de pinar reproyectada a UTM 31N
  matorral = "matorral.shp",           # mi capa de matorral reproyectada a UTM 31N
  margas   = "margas_U7.shp"           # mi capa de margas con calizas margosas
)

# Defino los parámetros de mi análisis
funct_principal <- "t"                 # utilizo la función de Tobler peatonal
paso_isocrona   <- 0.5                 # mi ancho de banda en horas
dir.create("movecost_resultados", showWarnings = FALSE)

# --- Mi proceso de carga y blindaje (para evitar que se me cuelgue) ---
dtm <- raster(ruta_mdt)
cat("MDT:", ncell(dtm), "celdas | resolucion:", res(dtm)[1], "m\n")

cargar_pts <- function(ruta){
  x <- st_read(ruta, quiet = TRUE)
  # Convierto mis puntos MULTIPOINT a POINT por si los guardé así en QGIS
  if (any(grepl("MULTI", st_geometry_type(x)))) {
    x <- suppressWarnings(st_cast(x, "POINT"))
  }
  x <- st_transform(x, crs_proyecto)
  p <- as(x, "Spatial")
  crs(p) <- crs(dtm); p
}
origen   <- cargar_pts(ruta_origen)
destinos <- cargar_pts(ruta_destinos)
if(is.null(destinos$name)) destinos$name <- paste0("destino_", seq_len(length(destinos)))

# Configuro mi barrera opcional para el mar.
# IMPORTANTE: Me aseguro de no usar NA en mi MDT para el mar, ya que
# me rompe las isocronas. Siempre uso mi polígono aquí.
barrera <- NULL
if(!is.na(ruta_barrera)){
  barrera <- as(st_transform(st_read(ruta_barrera, quiet=TRUE), crs_proyecto), "Spatial")
  crs(barrera) <- crs(dtm)
}

# Compruebo que mis puntos caigan sobre celdas válidas de mi MDT
stopifnot(all(!is.na(extract(dtm, origen))),
          all(!is.na(extract(dtm, destinos))))
cat("Puntos OK dentro del MDT.\n")

# ---------------------------------------------------------------------
#  PASO 2 — MI LOGÍSTICA DE ENTRADA: así calculo mi captación de recursos
#  Genero mi superficie de coste desde mi origen, mis isocronas y mis estadísticas
# ---------------------------------------------------------------------
cat("\n>> ENTRADA: superficie de coste y captacion de recursos\n")
ent  <- movecost(dtm = dtm, origin = origen, funct = funct_principal,
                 time = "h", barrier = barrera, graph.out = FALSE)
cost <- ent$accumulated.cost.raster
writeRaster(cost, "movecost_resultados/coste_acumulado_h.tif", overwrite = TRUE)

# Creo mis propias isocronas reclasificando el raster para tener un control total.
# Como movecost me devuelve Infinito en zonas inalcanzables, me quedo solo con mis valores finitos.
vals_cost <- values(cost)
max_finite <- max(vals_cost[is.finite(vals_cost)], na.rm = TRUE)
topo   <- ceiling(max_finite / paso_isocrona) * paso_isocrona
cortes <- seq(0, topo, by = paso_isocrona)
rcl    <- cbind(cortes[-length(cortes)], cortes[-1], seq_len(length(cortes)-1))
bandas <- reclassify(cost, rcl)
cat(sprintf("Coste max alcanzable: %.2f h (rampa hasta %.1f h, %d bandas)\n",
            max_finite, topo, length(cortes)-1))
writeRaster(bandas, "movecost_resultados/isocronas_banda.tif", overwrite = TRUE)
cell_ha <- prod(res(cost)) / 10000     # calculo mis hectáreas por celda

# Para cada uno de mis recursos, extraigo el tiempo mínimo de acceso y mis hectáreas
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
#  PASO 3 — MI LOGÍSTICA DE SALIDA: aquí trazo mis rutas de exportación
#  Calculo las rutas (LCP) desde mi origen a cada destino TERRESTRE.
#  El tramo marítimo no lo calculo con movecost, así que lo sumaré
#  más adelante como distancia / velocidad de navegación.
# ---------------------------------------------------------------------
cat("\n>> SALIDA: rutas de menor coste (LCP)\n")
sal <- movecost(dtm = dtm, origin = origen, destin = destinos,
                funct = funct_principal, time = "h",
                barrier = barrera, graph.out = FALSE, export = TRUE)
salida <- sal$dest.loc.w.cost@data
salida$dist_km <- round(sal$LCPs@data$length / 1000, 2)
salida$tiempo_h <- round(salida$cost, 2)
print(salida[, c("name","dist_km","tiempo_h","cost_hms")])
write.csv(salida[, c("name","dist_km","tiempo_h","cost_hms")],
          "movecost_resultados/rutas_salida.csv", row.names = FALSE)
shapefile(sal$LCPs, "movecost_resultados/LCP_salida.shp", overwrite = TRUE)

# ---------------------------------------------------------------------
#  PASO 4 — MI CONTRASTE DE ROBUSTEZ
#  Voy a repetir el cálculo de mis rutas con distintas funciones de coste
#  para comparar y asegurarme de que mis resultados son sólidos.
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
    " - coste_acumulado_h.tif / isocronas_banda.tif  (mis rasters para QGIS)\n",
    " - captacion_recursos.csv     (mis hectáreas de recurso por isocrona)\n",
    " - rutas_salida.csv + LCP_salida.shp  (mis tablas y mapas de salida)\n",
    " - robustez_funciones.csv     (mi contraste de algoritmos)\n", sep="")

# ---------------------------------------------------------------------
#  MIS NOTAS
#  1. El cabotaje lo calculo en mi PASO 5 usando gdistance.
#  2. No tengo una función calibrada para el asno, así que usaré "wcs"
#     como mi modelo de contraste.
#  3. Recuerdo usar siempre un polígono de barrera para el mar, nunca NAs.
# ---------------------------------------------------------------------


# =====================================================================
#  PASO 5 — MI LOGÍSTICA MARÍTIMA Y MODELO BIMODAL
#  Aquí calculo el cabotaje integrando la fricción de mi viento y oleaje,
#  y luego lo sumo a mi tramo terrestre. Utilizo la metodología calibrada
#  para mi embarcación Gyptis.
# =====================================================================
library(gdistance); library(ncdf4); library(Matrix)

# ---- MIS INPUTS ----
ruta_nc          <- "datos_marinos.nc"     # mi NetCDF con los datos marinos
ruta_embarque    <- "pto_embarque.shp"     # mi punto de origen en la costa
ruta_desembarque <- "pto_desembarque.shp"  # mi punto de destino en la costa
nombre_pto_terr  <- "EsCodolar"            # mi destino terrestre cuyo tiempo sumaré al marítimo
var_u  <- "u"                              # los nombres de mis variables
var_v  <- "v"
var_hs <- "hs"

# ---- MIS PARÁMETROS FÍSICOS (basados en la Gyptis) ----
vel_nudos      <- 3.5                            # mi velocidad media
vel_ms         <- vel_nudos * 0.5144             # lo paso a m/s
coste_base_h_m <- 1 / vel_ms / 3600              # mi tiempo base por metro
umbral_ola_m   <- 1.5                            # mi límite operativo de oleaje
mult_ola       <- 5                              # mi multiplicador si supero el umbral

# Mi límite operativo de viento: mi embarcación no soporta más de 15 nudos
vmax_viento_kn <- 15
vmax_viento_ms <- vmax_viento_kn * 0.5144        # = 7.72 m/s

# Calculo mi visibilidad costera. El alcance útil desde mi embarcación
# es la suma de mi visibilidad y la visibilidad de la costa.
h_observador_m <- 1                              # la altura de mis ojos sobre la cubierta
h_costa_m      <- 50                             ###<-- EDITAR: la altura media de mi costa visible
vis_observador_km <- sqrt(12.756 * h_observador_m)
vis_costa_km      <- sqrt(12.756 * h_costa_m)
dist_max_m        <- (vis_observador_km + vis_costa_km) * 1000

# Defino mi rampa de penalización por la dirección del viento
ramp_angulos <- c(0, 45, 90, 135, 180)
ramp_mult    <- c(1, 1.5, 2.5, 3.5, 4)

# ---- 1) NETCDF: cargo y reproyecto mis datos a la rejilla de mi MDT ----
u_r  <- raster(ruta_nc, varname = var_u)
v_r  <- raster(ruta_nc, varname = var_v)
hs_r <- raster(ruta_nc, varname = var_hs)
if (is.na(crs(u_r))) {                            # me aseguro de asignar WGS84 si vienen sin proyección
  crs(u_r) <- crs(v_r) <- crs(hs_r) <- CRS("+init=epsg:4326")
}
target <- raster(dtm)                              # utilizo la misma rejilla que en mi MDT terrestre
u_r  <- projectRaster(u_r,  target, method = "bilinear")
v_r  <- projectRaster(v_r,  target, method = "bilinear")
hs_r <- projectRaster(hs_r, target, method = "bilinear")

# ---- 2) FRICCIÓN BASE + penalización por mi oleaje ----
fric_base <- target;  values(fric_base) <- coste_base_h_m
penal_ola <- calc(hs_r, fun = function(x)
              ifelse(!is.na(x) & x > umbral_ola_m, mult_ola, 1))
fric_olas <- fric_base * penal_ola                 # mi penalización de horas por metro en cada celda

# ---- 3) MI MÁSCARA DE CABOTAJE: evalúo la visibilidad y el viento --
# Solo navego donde es mar, la costa es visible y el viento me lo permite.
tierra <- dtm
values(tierra) <- ifelse(values(dtm) > 0, 1, NA)
dist_costa <- raster::distance(tierra)             # calculo mi distancia hacia la costa
mar_ok <- is.na(values(tierra)) & (values(dist_costa) <= dist_max_m)
# Aplico mi tope operativo de viento
wind_speed_ms <- overlay(u_r, v_r, fun = function(u, v) sqrt(u^2 + v^2))
mar_ok <- mar_ok & (!is.na(values(wind_speed_ms))) & (values(wind_speed_ms) <= vmax_viento_ms)
mask_mar <- target;  values(mask_mar) <- ifelse(mar_ok, 1, NA)
fric_olas <- fric_olas * mask_mar                  # aplico NA si salgo de mi corredor navegable
cat(sprintf("Corredor cabotaje: visibilidad max %.1f km (obs=%.0f m + costa=%.0f m); viento max %.1f kn\n",
            dist_max_m/1000, h_observador_m, h_costa_m, vmax_viento_kn))

# ---- 4) DIRECCIÓN DEL VIENTO (hacia dónde sopla en mi modelo) ----
# Calculo mi penalización respecto a esta dirección.
wind_to <- overlay(u_r, v_r,
                   fun = function(u, v) (atan2(u, v) * 180 / pi) %% 360)

# ---- 5) MI MATRIZ DE TRANSICIÓN ANISOTRÓPICA -----------------------------
# 5.1 Configuro mi base isotrópica
conductancia <- 1 / fric_olas
tr_base <- transition(conductancia, transitionFunction = mean,
                      directions = 16, symm = FALSE)
tr_base <- geoCorrection(tr_base, type = "c")      # corrijo mi distancia entre celdas

# 5.2 Aplico mi ajuste vectorizado por el viento
m   <- transitionMatrix(tr_base)
tri <- summary(m)                                  # extraigo mis valores no nulos
xy_from <- xyFromCell(fric_olas, tri$i)
xy_to   <- xyFromCell(fric_olas, tri$j)
# Calculo el rumbo de mi barco
bearing <- (atan2(xy_to[,1] - xy_from[,1],
                  xy_to[,2] - xy_from[,2]) * 180 / pi) %% 360
# Extraigo el viento en el punto medio
wd <- raster::extract(wind_to, (xy_from + xy_to) / 2)
# Reduzco mi diferencia angular
diff_ang <- abs(((bearing - wd + 180) %% 360) - 180)
# Aplico mi multiplicador de coste por viento
bins_diff <- c(-Inf, 22.5, 67.5, 112.5, 157.5, Inf)
idx_bin   <- findInterval(diff_ang, bins_diff)
mult      <- ramp_mult[idx_bin]
mult[is.na(mult)] <- max(ramp_mult)                 # si no tengo datos, soy conservador y penalizo al máximo
# Ajusto mi conductancia
m_aniso <- sparseMatrix(i = tri$i, j = tri$j,
                        x = tri$x / mult, dims = dim(m))
tr_aniso <- tr_base
transitionMatrix(tr_aniso) <- m_aniso

# ---- 6) MI LCP MARÍTIMO (calculo el tiempo y la geometría) ---------------------
emb <- as(st_transform(st_read(ruta_embarque,    quiet = TRUE), crs_proyecto), "Spatial")
des <- as(st_transform(st_read(ruta_desembarque, quiet = TRUE), crs_proyecto), "Spatial")

# Si mi punto cae fuera del agua, uso esta función para engancharlo
# a mi celda de mar válida más cercana.
snap_to_sea <- function(pts, fric_r) {
  ok_idx <- which(!is.na(values(fric_r)))
  if (!length(ok_idx)) stop("No hay celdas de mar validas en fric_olas.")
  ok_xy  <- xyFromCell(fric_r, ok_idx)
  pts_xy <- coordinates(pts)
  for (i in seq_len(nrow(pts_xy))) {
    d <- sqrt((ok_xy[,1] - pts_xy[i,1])^2 + (ok_xy[,2] - pts_xy[i,2])^2)
    j <- which.min(d)
    if (d[j] > 0) {
      cat(sprintf("  Snap %s -> celda mar mas cercana (desplazamiento %.0f m)\n",
                  deparse(substitute(pts)), d[j]))
    }
    pts@coords[i, ] <- ok_xy[j, ]
  }
  pts
}
emb <- snap_to_sea(emb, fric_olas)
des <- snap_to_sea(des, fric_olas)

tiempo_mar_h <- as.numeric(costDistance(tr_aniso, emb, des))
if (!is.finite(tiempo_mar_h))
  stop("costDistance devuelve Inf: mi rejilla marina no conecta mis puntos.\n",
       "  Posibles causas: mi barrera ahoga las celdas, o mi tope de viento aísla la ruta.")
lcp_mar      <- shortestPath(tr_aniso, emb, des, output = "SpatialLines")
crs(lcp_mar) <- crs(dtm)
cat(sprintf(">> Tiempo maritimo (cabotaje): %.2f h\n", tiempo_mar_h))

# ---- 7) MI INTEGRACIÓN BIMODAL (sumo lo terrestre y lo marítimo) --------------------
tiempo_tier_h <- salida$tiempo_h[salida$name == nombre_pto_terr]
if (length(tiempo_tier_h) == 0)
  stop("No encuentro el destino terrestre '", nombre_pto_terr,
       "' en mi salida (Paso 3). Reviso mi shape de destinos.")

coste_bimodal <- data.frame(
  embarque          = nombre_pto_terr,
  tramo_terrestre_h = round(tiempo_tier_h, 2),
  tramo_maritimo_h  = round(tiempo_mar_h,  2),
  total_bimodal_h   = round(tiempo_tier_h + tiempo_mar_h, 2)
)
print(coste_bimodal)

# ---- 8) EXPORTACIÓN DE MIS RESULTADOS ---------------------------------------------------
write.csv(coste_bimodal, "movecost_resultados/coste_bimodal.csv", row.names = FALSE)
shapefile(lcp_mar,       "movecost_resultados/LCP_maritimo.shp", overwrite = TRUE)
cat("Listo: he guardado mis resultados finales en mi carpeta movecost_resultados/\n")

```
