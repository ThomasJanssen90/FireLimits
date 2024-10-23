#### This R script contains the analysis loop to attribute fire stops to fire perimeter points in North East Siberia
#### Created by Thomas Janssen, contact: t.a.j.janssen@vu.nl

### load required libraries
library(reshape2)
library(lubridate)
library(strucchange)
library(data.table)
library(terra)
library(spatstat)
library(parallel)
library(data.table)
library(dplyr)
library(sf)

rm(list = ls()) ### empty everything

sourcedir = "/gpfs/home6/janssent/Rscripts/Firebarriers/Analysis/SourceFunctions" ## directory with source functions
datadir = "/gpfs/work2/0/einf3869/Firebarriers/" ## directory containing datasets

source(file.path(sourcedir, "Bearing_and_Directions.R"))## functions calculating bearing and directions
source(file.path(sourcedir, "Wind_Functions.R"))## functions calculating wind fire spread index
source(file.path(sourcedir, "VPD_Functions.R"))## functions calculating vapor pressure deficit change
source(file.path(sourcedir, "SoilMoisture_Functions.R"))## functions calculating soil moisture change
source(file.path(sourcedir, "Landscape_Functions.R"))## functions calculating landscape drivers
    
print("finished reading in functions, continue with files")

files = list.files(file.path(datadir,"/BurnedArea/Derived/PerimeterPoints/CSV/"), full.names = T, pattern = ".csv$") ### Read file names with all perimeter points

filesvpd = data.frame(files = list.files(file.path(datadir,"/ECMWF/VPD/"), full.names = T, pattern = ".tif$"))
filesvpd$year = sapply(strsplit(basename(filesvpd$files), "_"), `[`, 3)
filesvpd$month = sapply(strsplit(basename(filesvpd$files), "_"), `[`, 4)

filesws = data.frame(files =list.files(file.path(datadir, "/ECMWF/Wind/Windspeed/"), full.names = T, pattern = ".tif$"))
filesws$year = sapply(strsplit(basename(filesws$files), "_"), `[`, 3)
filesws$month = sapply(strsplit(sapply(strsplit(basename(filesws$files), "_"), `[`, 4), ".tif"), `[`, 1)

fileswd = data.frame(files =list.files(file.path(datadir,"/ECMWF/Wind/Winddirection"), full.names = T, pattern = ".tif$"))
fileswd$year = sapply(strsplit(basename(fileswd$files), "_"), `[`, 3)
fileswd$month = sapply(strsplit(sapply(strsplit(basename(fileswd$files), "_"), `[`, 4), ".tif"), `[`, 1)

filessm = data.frame(files =list.files(file.path(datadir,"/ECMWF/SoilMoisture"), full.names = T, pattern = ".tif$"))
filessm$year = sapply(strsplit(basename(filessm$files), "_"), `[`, 3)
filessm$month = sapply(strsplit(sapply(strsplit(basename(filessm$files), "_"), `[`, 4), ".tif"), `[`, 1)

print("start loop")

for(j in 1:length(files)){
  files = list.files(file.path(datadir, "/BurnedArea/Derived/PerimeterPoints/CSV/"), full.names = T, pattern = ".csv$")## Read file names with all perimeter points
  
  ### Read in the perimeter point locations and add the values to the points
  data = fread(files[j])
  data = data[order(as.Date(data$burndate)),] ## order by burn date
  data$year = format(as.Date(data$burndate), "%Y") ## define year
  year = unique(data$year) 
  print(paste("start year ", year, sep = ""))
  
  data$month = format(as.Date(data$burndate), "%m") ## define month
  perimeters = vect(data, geom = "geometry", crs = "EPSG:3576")## convert data table to spatial vector
  perimeters = project(perimeters, "EPSG:4326") ## project points to EPSG:4326 WGS84 - World Geodetic System 1984
    
  perimeters$lat = crds(perimeters)[,2] ## get perimeter latitude
  perimeters$lon = crds(perimeters)[,1] ## get perimeter longitude
  
  xs = seq(80,175,5) ## define longitude range in 5 degree tiles
  ys = seq(40,70,5) ## define latitude range in 5 degree tiles
  tiles = expand.grid(xs, ys) ## expand grid to include all combinations
  colnames(tiles) = c("xmin", "ymin")
  tiles$xmax = tiles$xmin+5
  tiles$ymax = tiles$ymin+5 
  tiles$xminb = tiles$xmin-10 ## define buffer around tile of 10 degrees
  tiles$yminb = tiles$ymin-10
  tiles$xmaxb = tiles$xmax+10
  tiles$ymaxb = tiles$ymax+10
  tiles$tile = paste(tiles$xmin,tiles$xmax,tiles$ymin,tiles$ymax, sep = "_") ## define tile name
  
  ## start loop for every tile  
  for(i in 1:nrow(tiles)){
    #i = 43 #use example tile 43
    print(paste("start tile",i,"out of",nrow(tiles),sep =" "))
    
    #filesexist = list.files(file.path(datadir, "/Perimeter/Attributed/"))
    #if(paste("Landscape_attributed_year_",year, "_tile_",tiles$tile[i],".csv",sep = "") %in% filesexist) {next}
    
    #subset perimeters inside tile i
    s = subset(perimeters, perimeters$lat >= tiles$ymin[i] & perimeters$lat <= tiles$ymax[i] & perimeters$lon >= tiles$xmin[i] & perimeters$lon <= tiles$xmax[i])
    if(nrow(s) < 1) {next} #if there is no data, skip to next
    s = data.table(as.data.frame(s)) ## create data table from perimeters subset 
    s$burndate = as.Date(s$burndate)  ## burn date in Date format
        
    n.cores = 128 ## define number of cores for parallel processing
    
    print(paste("start doing wfsi analysis",year, tiles$tile[i],sep = " " ))
    data = mclapply(1:nrow(s), FUN = funws, mc.cores = n.cores) ## wind fire spread index analysis in parallel
    print(paste("finished wfsi analysis",year, tiles$tile[i],sep = " " ))
    
    ### convert returned values to data table
    wfsi_list = vector("list", nrow(s))
    metdata_list = vector("list", nrow(s))
    pval_list = vector("list", nrow(s))
    for(l in 1:nrow(s)){
      wfsi_list[[l]] = data[[l]][[2]]
      metdata_list[[l]] = data[[l]][[1]]
      pval_list[[l]] = data[[l]][[3]]}
      
    # Update 's$wfsi' with the list of wfsi values
    s$wfsi = unlist(wfsi_list)
    s$wfsip = unlist(pval_list) 
    # Combine the list of metdata elements into a single data.table
    metdata = rbindlist(metdata_list)
    
    fwrite(metdata, paste(datadir, "/Perimeter/Extracted/WFSI/Extracted_WFSI_", year, "_tile_",tiles$tile[i], ".csv", sep = ""))
    rm(data,metdata, wfsi_list,metdata_list)
    print(paste("finished wfsi extracting",year, tiles$tile[i],sep = " " ))
             
    print(paste("start doing vpd analysis",year, tiles$tile[i],sep = " " ))
    data = mclapply(1:nrow(s), FUN = funvpd, mc.cores = n.cores) ## vapor pressure deficit analysis in parallel
    print(paste("finished vpd analysis",year, tiles$tile[i],sep = " " ))
    
    vpd_list = vector("list", nrow(s))
    metdata_list = vector("list", nrow(s))
    pval_list = vector("list", nrow(s))
    for(l in 1:nrow(s)){
      vpd_list[[l]] = data[[l]][[2]]
      metdata_list[[l]] = data[[l]][[1]]
      pval_list[[l]] = data[[l]][[3]]}
      
    # Update 's$wfsi' with the list of wfsi values
    s$vpd = unlist(vpd_list)
    s$vpdp = unlist(pval_list) 
    # Combine the list of metdata elements into a single data.table
    metdata = rbindlist(metdata_list)
    
    fwrite(metdata, paste(datadir, "/Perimeter/Extracted/VPD/Extracted_VPD_", year, "_tile_",tiles$tile[i], ".csv", sep = ""))
    rm(data,metdata,vpd_list,metdata_list)
    print(paste("finished VPD extracting",year, tiles$tile[i],sep = " " ))
           
    print(paste("start doing sm analysis",year, tiles$tile[i],sep = " " ))
    data = mclapply(1:nrow(s), FUN = funsm, mc.cores = n.cores) # soil moisture analysis in parallel
    print(paste("finished sm analysis",year, tiles$tile[i],sep = " " ))
    
    sm_list = vector("list", nrow(s))
    metdata_list = vector("list", nrow(s))
    pval_list = vector("list", nrow(s))
    for(l in 1:nrow(s)){
      sm_list[[l]] = data[[l]][[2]]
      metdata_list[[l]] = data[[l]][[1]]
      pval_list[[l]] = data[[l]][[3]]}
      
    # Update 's$wfsi' with the list of wfsi values
    s$sm = unlist(sm_list)
    s$smp = unlist(pval_list) 
    # Combine the list of metdata elements into a single data.table
    metdata = rbindlist(metdata_list)
    
    fwrite(metdata, paste(datadir, "/Perimeter/Extracted/SM/Extracted_SM_", year, "_tile_",tiles$tile[i], ".csv", sep = ""))
    rm(data,metdata, sm_list, metdata_list)
    print(paste("finished SM extracting",year, tiles$tile[i],sep = " " ))
        
    ### This completes the fire weather (temporal) analysis, from here we start the landscape drivers (spatial) analysis    
      
    ### Load burned area data
    files = list.files(file.path(datadir, "/BurnedArea/Harmonized/"), full.names = T)
    db = vrt(files[grep(year, files)])
             
    year = as.numeric(year)
    
    ### Load ESA CCI above-ground biomass data
    if(year>2021) {filelist = data.frame(files = list.files(paste(datadir, "/Biomass/ESACCI-BIOMASS/Raw_years/",year-2, sep = ""), full.names = T))
    filelist$ymax = as.numeric(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1),start = 2, stop =3))
    filelist$ymax = ifelse(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1),start = 1, stop =1) == "S", filelist$ymax*-1, filelist$ymax)
    filelist$ymin = filelist$ymax-10
    filelist$xmin = as.numeric(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1), start = 5, stop = 7))
    filelist$xmin = ifelse(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1),start = 4, stop =4) == "W", filelist$xmin*-1, filelist$xmin)
    filelist$xmax = filelist$xmin+10}
    
    ### Load ESA CCI above-ground biomass data for years in between data
    if(year<2018) {filelist = data.frame(files = list.files(paste(datadir, "/Biomass/ESACCI-BIOMASS/Calculated_years/",year-1, sep = ""), full.names = T))
    filelist$ymax = as.numeric(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 3),start = 2, stop =3))
    filelist$ymin = filelist$ymax-10
    filelist$xmin = as.numeric(sapply(strsplit(sapply(strsplit(basename(filelist$files), "_"),`[`, 3), "E"),`[`, 2))
    filelist$xmax = filelist$xmin+10}
    
    if(year>=2018 & year <= 2021){
    filelist = data.frame(files = list.files(paste(datadir, "/Biomass/ESACCI-BIOMASS/Raw_years/",year-1, sep = ""), full.names = T))
    filelist$ymax = as.numeric(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1),start = 2, stop =3))
    filelist$ymax = ifelse(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1),start = 1, stop =1) == "S", filelist$ymax*-1, filelist$ymax)
    filelist$ymin = filelist$ymax-10
    filelist$xmin = as.numeric(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1), start = 5, stop = 7))
    filelist$xmin = ifelse(substr(sapply(strsplit(basename(filelist$files), "_"),`[`, 1),start = 4, stop =4) == "W", filelist$xmin*-1, filelist$xmin)
    filelist$xmax = filelist$xmin+10}
    
    ### Pick spatial subset of biomass files
    pickfiles = subset(filelist, filelist$xmin >= tiles$xminb[i] & filelist$xmax <= tiles$xmaxb[i] & filelist$ymin >= tiles$yminb[i] & filelist$ymax <= tiles$ymaxb[i]) 
    agb = vrt(pickfiles$files, overwrite = T, filename="agb") ### Read biomass data
    
    ### Read Tree cover files        
    files = list.files(file.path(datadir, "/Treecover/Hansen/TreeCoverCalculated"), full.names = T)
    filelist = data.table(files = files[grep(as.character(as.numeric(year)-1), files)]) ### Pick only files in the previous year
    filelist$xmin = as.numeric(sapply(strsplit(sapply(strsplit(basename(filelist$files), "_"),`[`, 6), "E"),`[`, 1))
    filelist$xmax = filelist$xmin+10
    filelist$ymax = as.numeric(sapply(strsplit(sapply(strsplit(basename(filelist$files), "_"),`[`, 5), "N"),`[`, 1))
    filelist$ymin = filelist$ymax-10
    pickfiles = subset(filelist, filelist$xmin >= tiles$xminb[i] & filelist$xmax <= tiles$xmaxb[i] & filelist$ymin >= tiles$yminb[i] & filelist$ymax <= tiles$ymaxb[i]) 
    ### Pick only files within the tile + buffer zone
    treecoverr = vrt(pickfiles$files, overwrite = T, filename="treecoverr")
    
    ### List and pick the road data files
    files = list.files(file.path(datadir,"/Roads/Merged/Raster/"), full.names = T)
    latlon = sapply(strsplit(sapply(strsplit(basename(pickfiles$files), "year_"),`[`, 2), "tree"),`[`, 1)
    files = files[grep(paste(latlon, collapse = "|"),files)]
    
    roads = vrt(files, overwrite = T, filename="roads") ### Read road data
    
    ### List and pick  the global surface water files
    if(year<2022) {filelist = fread(file.path(datadir, "/GSW/YearlyClassification/GSW_tile_list.csv"))
    filelist = filelist[grep(year,filelist$name),]
    pickfiles = subset(filelist, filelist$xmin >= tiles$xminb[i] & filelist$xmax <= tiles$xmaxb[i] & filelist$ymin >= tiles$yminb[i] & filelist$ymax <= tiles$ymaxb[i])}
    
    if(year == 2022) {filelist = fread(file.path(datadir, "/GSW/YearlyClassification/GSW_tile_list.csv"))
    filelist = filelist[grep(year-1,filelist$name),]
    pickfiles = subset(filelist, filelist$xmin >= tiles$xminb[i] & filelist$xmax <= tiles$xmaxb[i] & filelist$ymin >= tiles$yminb[i] & filelist$ymax <= tiles$ymaxb[i])}
    
    gsw = vrt(pickfiles$name, overwrite = T, filename="gsw") ### Read GSW data
    
    ### List and pick files DEM
    filelist = fread(file.path(datadir, "/JAXA/DEM/ALOS_DEM_tile_list.csv"))
    pickfiles = subset(filelist, filelist$xmin >= tiles$xminb[i] & filelist$xmax <= tiles$xmaxb[i] & filelist$ymin >= tiles$yminb[i] & filelist$ymax <= tiles$ymaxb[i])     
    demr = vrt(pickfiles$name, overwrite = T, filename="demr") ### Read DEM
    
    #### List and pick files Population
    filelist = fread(file.path(datadir, "/Population/GHS_POP_E2015_GLOBE/Pop_tile_list.csv"))
    pickfiles = subset(filelist, filelist$xmin >= tiles$xminb[i] & filelist$xmax <= tiles$xmaxb[i] & filelist$ymin >= tiles$yminb[i] & filelist$ymax <= tiles$ymaxb[i])
    population = vrt(pickfiles$name, overwrite = T, filename="pop")
    
    ### Read in burn history and land cover data
    filesbh = list.files(file.path(datadir, "/BurnedArea/Derived/BurnHistory/"), full.names = T)
    bh = vrt(filesbh[grep(as.character(as.numeric(year)-1), filesbh)])    
    lc = vrt(file.path(datadir, "/LandCover/Bartalev/Projected/Projected_Bartalev_Merged_Landcover.tif"))
  
    print("start doing the analysis")
    n.cores = 128
    data = mclapply(1:nrow(s), FUN = extrdata, mc.cores = n.cores) ### Do the landscape drivers analysis
    print("finished the analysis")
    data = do.call(rbind, data)
                   
    s = cbind(s, data[,2:ncol(data)])
    
    ### Write the attributed perimeter points to a csv file
    fwrite(s, paste(datadir,"/Perimeter/Attributed/Landscape_attributed_year_", year, "_tile_",tiles$tile[i], ".csv", sep = ""))
          
    print(paste(i, " out of ", nrow(tiles),sep = ""))
    gc()
    tmpFiles(remove = T) 
    }}
  
  
  
  
  

