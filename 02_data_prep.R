

# package and environment setup
library(tidyverse)
library(raster)
library(geosphere)
library(caret)
library(ggmap)

select <- dplyr::select

setwd("e:/seeds")


###### restoration parameters #####

# a hypothetical restoration site in the Presidio
site <- data.frame(species=NA, lat=37.801064, lon=-122.478557)
coordinates(site) <- c("lon", "lat")

# list of focal species
spp <- read.csv("focal_species.csv", stringsAsFactors = F)
spp <- spp$Species


###### background data #####

# species occurrence data
#bien <- read.csv("e:/BIEN/bien_with_chelsa.csv", stringsAsFactors = F) %>%
#      filter(scrubbed_species_binomial %in% spp)
jeps <- read.csv("E:/phycon/data/occurrences/California_Species_clean_All_epsg_3310.csv", stringsAsFactors=F) %>%
      filter(current_name_binomial %in% spp) %>%
      select(current_name_binomial, longitude, latitude) %>%
      rename(species=current_name_binomial)
coordinates(jeps) <- c("longitude", "latitude")

# projections
wgs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
crs(site) <- crs(jeps) <- wgs 

# soils data
soil <- list.files("F:/SoilGrids", full.names=T, pattern=".tif") %>%
      stack() %>%
      crop(extent(jeps))
sv <- values(soil)
na <- apply(sv, 1, function(x) any(is.na(x)))
sva <- sv[!na,]
pca <- prcomp(sva, center=T, scale.=T)
npcs <- 5
svpc <- sv[,1:npcs]
svpc[!na,] <- pca$x[,1:npcs]
soil_pc <- soil[[1:npcs]]
soil_pc[] <- svpc

jeps_soil <- raster::extract(soil_pc, jeps)



# climate data

clim <- data.frame(path=list.files("F:/seeds/chelsa_derived", recursive=T, full.names=T),
                   stringsAsFactors=F) 

ext <- extent(spTransform(jeps, crs(raster(clim$path[1]))))

jeps_clim <- raster::extract(stack(clim$path[1]), jeps)
colnames(jeps_clim) <- c("PPT", "PET", "AET", "CWD", "RAR", "DJF", "JJA")
jeps_clim[,"PPT"] <- log10(jeps_clim[,"PPT"])
jeps_clim <- jeps_clim[,c("PPT", "AET", "CWD", "DJF", "JJA")]


lapply(clim$path, function(x){
      message(x)
      s <- stack(x) %>% crop(ext)
      names(s) <- c("PPT", "PET", "AET", "CWD", "RAR", "DJF", "JJA")
      s[["PPT"]] <- log10(s[["PPT"]])
      s <- s[[c("PPT", "AET", "CWD", "DJF", "JJA")]]
      writeRaster(s, paste0("e:/seeds/cloud/assets/seedsource_data/climate ", basename(x)),
                  overwrite=T)
})



md <- map_data("state", "California")

# save prepped data
d <- list(spp=spp,
          jeps=jeps,
          jeps_clim=jeps_clim,
          jeps_soil=jeps_soil,
          md=md)
saveRDS(d, "e:/seeds/cloud/assets/seedsource_data/species_data.rds")
writeRaster(soil_pc, "e:/seeds/cloud/assets/seedsource_data/soil.tif", overwrite=T)



stop("donezo")








########## deprecated below #########



####
clim <- data.frame(path=list.files("F:/BCM/Normals_30years", recursive=T, full.names=T)) %>%
      separate(path, c("junk", "md"), sep="BCM2014_", remove=F) %>%
      select(-junk) %>%
      separate(md, c("var", "md"), sep=3) %>%
      separate(md, c("year", "model"), sep="_wy_ave_") %>%
      mutate(model=gsub("\\.Rdata|_2", "", model)) %>%
      separate(model, c("model", "scenario"), sep="_") %>%
      mutate(scenario = ifelse(is.na(scenario), "HST", scenario)) %>%
      filter(grepl("HST|rcp", scenario))



clim <- list.files("E:/phycon/data/bcm_normals/Normals_30years", full.names=T, pattern="BCM2014_") %>%
      lapply(readRDS) %>%
      stack()
names(clim) <- c("AET", "CWD", "DJF", "JJA", "PPT", "TMN", "TMX")
clim1 <- clim[[c("CWD", "DJF", "JJA", "PPT")]]
clim1$PPT <- log10(clim1$PPT)

clim1 <- crop(clim1, extent(spTransform(jeps, crs(clim1))))

clim2 <- clim1 + c(200, 2, 1.5, 0) # fake climate change data
clim3 <- clim2 + c(200, 3, 2.5, 0)
clim <- list(clim1, clim2, clim3)

jeps_clim <- extract(clim1, jeps)





md <- map_data("state", "California")





# save prepped data
d <- list(spp=spp,
          jeps=jeps,
          jeps_clim=jeps_clim,
          jeps_soil=jeps_soil,
          #soil=soil_pc,
          #clim=clim,
          md=md)
saveRDS(d, "e:/seeds/app/species_data.rds")
writeRaster(soil_pc, "e:/seeds/app/soil.tif", overwrite=T)
writeRaster(clim1, "e:/seeds/app/clim1.tif")
writeRaster(clim2, "e:/seeds/app/clim2.tif")
writeRaster(clim3, "e:/seeds/app/clim3.tif")


stop("bingo")


function(sp){
      
      # geographic distances from restoration site
      j <- jeps[jeps$species==sp,]
      gdist <- as.vector(distm(site, j))
      
      # climate and soils data for restoration site
      sc <- lapply(clim, extract, y=site)
      ss <- extract(soil, site)
      
      # climate and soil data for species occurrences
      jc <- jeps_clim[jeps$species==sp,]
      js <- jeps_soil[jeps$species==sp,]
      
      # 3 climate PCs
      ctrans <- preProcess(jc, c("center", "scale", "pca"), pcaComp=3)
      tag <- function(d, tag){colnames(d) <- paste0(tag, colnames(d)); return(d)}
      sc <- lapply(sc, function(x) predict(ctrans, x) %>% tag("clim"))
      jc <- predict(ctrans, jc) %>% tag("clim")
      
      # 3 soil PCs
      strans <- preProcess(js, c("center", "scale", "pca"), pcaComp=3)
      ss <- predict(strans, ss) %>% tag("soil")
      js <- predict(strans, js) %>% tag("soil")
      
      # combine soil and climate data
      je <- cbind(jc, js)
      se <- lapply(sc, function(x) cbind(x, ss))
      
      # placeholder for spatial smoothing to represent gene flow?
      
      # environmental distances
      edist <- lapply(se, function(x){
            rbind(x, je) %>% dist() %>% as.matrix() %>% "["((2:nrow(.)), 1) }) %>%
            do.call("cbind", .)
      colnames(edist) <- paste0("env_dist", 1:ncol(edist))
      
      
      
      
      
      # convert to data frames for plotting
      sd <- as.data.frame(site)
      jd <- as.data.frame(j) %>%
            mutate(geo_dist = gdist) %>%
            cbind(edist)
      
      wm <- apply(jd, 1, function(x) as.integer(which.min(tail(x, 3))))
      wm[sapply(wm, length)==0] <- NA
      wm <- unlist(wm)
      jd$wm <- wm
      
      jdl <- jd %>%
            gather(time, env_dist, env_dist1:env_dist3)
      
      
      
      ggplot(jd, aes(longitude, latitude, color=env_dist1)) + 
            geom_point() +
            geom_point(data=sd, aes(lon, lat), color="red") +
            scale_color_viridis_c() +
            theme_void() +
            coord_fixed()
      
      ggplot(jd, aes(geo_dist, env_dist1)) + 
            geom_point() +
            annotate(geom="point", x=0, y=0, color="red") +
            theme_minimal()
      
      ggplot(jd, aes(longitude, latitude, color=factor(wm))) +
            geom_point() +
            theme_void() +
            coord_fixed()
      
      
      ggplot(jdl, aes(longitude, latitude, color=env_dist)) + 
            geom_point() +
            geom_point(data=sd, aes(lon, lat), color="red") +
            scale_color_viridis_c() +
            theme_void() +
            facet_wrap(~time, nrow=1) +
            coord_fixed()
      
      ggplot(jdl, aes(geo_dist, env_dist, color=wm)) + 
            geom_path(aes(geo_dist, env_dist, group=geo_dist)) +
            geom_point() +
            annotate(geom="point", x=0, y=0, color="black") +
            theme_minimal()
      
      
      ggplot(jd, aes(env_dist1, env_dist3, color=geo_dist)) +
            geom_abline(slope=1, intercept=0) +
            geom_point() +
            theme_minimal() +
            coord_fixed() +
            xlim(0, NA) + ylim(0, NA)
      
      
}
