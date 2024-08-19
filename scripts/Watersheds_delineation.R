#install.packages("whitebox")
#install.packages("terra")
#install.packages("stars")
library(whitebox)
#whitebox::install_whitebox()
wbt_init()
library(terra)
library(stars)

#archivo raster  del DEM (30 mts. resolucion)
dem_tif <- "E:/DATA_MAIN/laDelicia/cuencas/dem/1_clip.tif"
#dem_tif <- "E:/DATA_MAIN/laDelicia/dem/5mts_v2/mosaic_filt_5.tif"
#cortar?
cut_dem <- F
#shapefile del area de estudio paa cortar
areaEstudio_shp <- "E:/DATA_MAIN/laDelicia/areaEstudio/areaEstudio.shp"
#preparar DEM
prepare_DEM <- F
#computar cuencas?
compute_watersheds_streams <- T
#area en pixeles para cerrar cuencas
watershed_pix <- 10000
stream_pix <- 10000
#directorio de salida para el procesamiento del DEM
output_DEM <- "E:/DATA_MAIN/laDelicia/cuencas/dem"
#directorio de salida para los parametros de las cuencas
output_parameters <- "E:/DATA_MAIN/laDelicia/cuencas/w10k_s10k"

#### PREPARAR AREA ESTUDIO ####

#preparar directorios
dir.create(output_DEM,showWarnings=F)
dir.create(output_parameters,showWarnings=F)
#cortar
if(cut_dem){
  #leer area de estudio y calcular area
  areaEst.map <- st_read(areaEstudio_shp,quiet=T)
  #leer raster
  dem.ras <- rast(dem_tif)
  if(crs(dem.ras,proj=T)!=st_crs(areaEst.map)$proj4string){
    areaEst.map <- st_transform(areaEst.map,crs(dem.ras,proj=T))
  }
  #extraer
  areaEst.map <- vect(areaEst.map)
  dem.ras <- crop(dem.ras, ext(areaEst.map))
  dem.ras <- mask(dem.ras, areaEst.map)
  #proyectar raster
  dem.ras <- project(x=dem.ras,y=crs(vect(areaEstudio_shp),proj=T),method="bilinear")
  #grabar
  out.file <- paste0(output_DEM,"/1_clip.tif")
  writeRaster(dem.ras,out.file,NAflag=0,datatype = "INT2U",overwrite=T)
  in.file <- out.file
}else{
  in.file <- dem_tif
}

#### PREPARAR DEM ####

if(prepare_DEM){
  #rellenar 1 (breach)
  out.file <- paste0(output_DEM,"/2_clip_breach.tif")
  wbt_breach_depressions_least_cost(
    dem = in.file,
    output = out.file,
    dist = 5,
    fill = TRUE)
  #rellenar 2 (fill)
  in.file <- out.file
  out.file <- paste0(output_DEM,"/3_clip_breach_fill.tif")
  wbt_fill_depressions_wang_and_liu(
    dem = in.file,
    output = out.file
  )
  filled.dem <- out.file
  #derivar mapa de sombras
  out.file <- paste0(output_DEM,"/4_hillshade.tif")
  wbt_hillshade(dem = filled.dem,
                output = out.file,
                azimuth = 115)
  #derivar acumulacion de flujo
  fac.dem <- paste0(output_DEM,"/5_fac.tif")
  wbt_d8_flow_accumulation(input = filled.dem,
                           output = fac.dem)
  #derivar punteros de flujo
  pd8.dem <- paste0(output_DEM,"/6_pd8.tif")
  wbt_d8_pointer(dem = filled.dem,
                 output = pd8.dem)
  #derivar cuencas + vectorizar
  basins.dem <- paste0(output_DEM,"/7_basins.tif")
  wbt_basins(d8_pntr = pd8.dem,
             output = basins.dem)
  basins.vec <- rast(basins.dem)
  basins.vec <- st_as_stars(basins.vec) %>% st_as_sf(merge = T)
  names(basins.vec)[1] <- "basins"
  out.file <- paste0(output_DEM,"/7_basins.shp")
  if(file.exists(out.file)){
    st_write(basins.vec,out.file,quiet=T,delete_dsn=T)
  }else{
    st_write(basins.vec,out.file,quiet=T,delete_dsn=F)
  }
}

#### DERIVAR CUENCAS ####

if(compute_watersheds_streams){
  #leer inputs
  filled.dem <- paste0(output_DEM,"/3_clip_breach_fill.tif")
  fac.dem <- paste0(output_DEM,"/5_fac.tif")
  pd8.dem <- paste0(output_DEM,"/6_pd8.tif")
  #derivar cuencas homogeneas + vectorizar
  isoBasins.dem <- paste0(output_parameters,"/8_isoBasins.tif")
  wbt_isobasins(dem = filled.dem,
                size=watershed_pix, #area en pixeles
                output = isoBasins.dem)
  basins.vec <- rast(isoBasins.dem)
  basins.vec <- st_as_stars(basins.vec) %>% st_as_sf(merge = T)
  names(basins.vec)[1] <- "isoBasins"
  out.file <- paste0(output_parameters,"/8_isoBasins.shp")
  if(file.exists(out.file)){
    st_write(basins.vec,out.file,quiet=T,delete_dsn=T)
  }else{
    st_write(basins.vec,out.file,quiet=T,delete_dsn=F)
  }
  #derivar red de escurrimiento 
  streams.dem <- paste0(output_parameters,"/9_streams.tif")
  wbt_extract_streams(flow_accum = fac.dem,
                      output = streams.dem,
                      threshold = stream_pix)
  #encontrar canal principal + vectorizar
  mainStream.dem <- paste0(output_parameters,"/10_mainStream.tif")
  wbt_find_main_stem(d8_pntr= pd8.dem,
                     streams = streams.dem,
                     output = mainStream.dem)
  out.file <- paste0(output_parameters,"/10_mainStream.shp")
  wbt_raster_streams_to_vector(streams = mainStream.dem,
                               d8_pntr = pd8.dem,
                               output = out.file)
  #derivar order por Strahler + vectorizar
  strahler.dem <- paste0(output_parameters,"/11_strahler.tif")
  wbt_strahler_stream_order(d8_pntr= pd8.dem,
                            streams = streams.dem,
                            output = strahler.dem)
  out.file <- paste0(output_parameters,"/11_strahler.shp")
  wbt_raster_streams_to_vector(streams = strahler.dem,
                               d8_pntr = pd8.dem,
                               output = out.file)
}
