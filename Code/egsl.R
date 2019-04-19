# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(magrittr)
library(sf)
library(sp)
library(raster)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   DOWNLOAD DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For now the data is local, will need to be uploaded to Zenodo
# Output location for downloaded data
output <- './Data/RawData'


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    IMPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File name
fileName <- dir(output, pattern = '.zip')

# Unzip kmz file
unzip(zipfile = paste0(output, '/', fileName),
      exdir = output)

# Identify newly extracted files
fileName <- dir(output, pattern = 'shp$', full.names = T)

# Import shapefile file
egsl <- rgdal::readOGR(dsn = output, layer = "egsl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               SIMPLIFIED CONTOUR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Selected plygons in egsl (remove small islands that add many vertex to the dataset)
j <- c(1, 9127, 8999, 9038, 9005, 8998, 131, 18, 17)
p <- list()
for(i in 1:length(j)) {
  p[[i]] <- Polygon(egsl@polygons[[1]]@Polygons[[j[i]]]@coords)
}
ps <- Polygons(p, 1)
sps <- SpatialPolygons(list(ps))
proj4string(sps) <- proj4string(egsl)

# Transform to sf object and simplify
egslSimple <- sps %>%
              SpatialPolygonsDataFrame(. , data.frame(ID = 1)) %>%
              st_as_sf() %>%
              st_transform(crs = 32198) %>%
              st_simplify(preserveTopology = T, dTolerance = 1000)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               SIMPLIFIED CONTOUR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
egslCoast <- st_boundary(egslSimple)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           1KM^2 HEXAGONAL VECTOR GRID
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cell size
areaKM2 = 1 # area of cells in grid
A <- areaKM2 * 1000000 # area in m^2
cellsize <- sqrt(2*A)/3^(1/4) # Short diagonal used to calculate cell size

# Create grid
egslGrid <- egsl %>%
            # Buffer around whole area to make sure that the complete area is selected
            rgeos::gBuffer(width = cellsize) %>%

            # Sample points for hexagonal grid
            sp::spsample(type = "hexagonal",
                         cellsize = cellsize,
                         offset = c(0,0)) %>%

            # Create polygons
            sp::HexPoints2SpatialPolygons() %>%

            # Transform as sf object
            sf::st_as_sf()

# Add ID
egslGrid$ID <- paste0('ID',1:nrow(egslGrid))

# Transform egsl as sf object from here
egsl <- sf::st_as_sf(egsl)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                1KM RASTER GRID
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters
extGrid <- extent(egsl)
crsGrid <- projection(egsl)
resGrid <- 1000

# Grid
rasterGrid <- raster(ext = extGrid,
                     crs = crsGrid,
                     res = c(resGrid, resGrid))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         1.5KM RASTER GRID (shiny app)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters
extGrid <- extent(egsl)
crsGrid <- projection(egsl)
resGrid <- 1500

# Grid
appGrid <- raster(ext = extGrid,
                  crs = crsGrid,
                  res = c(resGrid, resGrid))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                     EXPORT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export
save(egslSimple, file = './data/egslSimple.RData')
save(egslGrid, file = './data/egslGrid.RData')
save(egslCoast, file = './data/egslCoast.RData')
save(rasterGrid, file = './data/rasterGrid.RData')
save(appGrid, file = './data/appGrid.RData')
