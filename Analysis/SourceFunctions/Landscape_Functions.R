#### This R script contains the code to attribute landscape drivers as fire stops to fire perimeter points in North East Siberia
#### Created by Thomas Janssen, contact: t.a.j.janssen@vu.nl
    
### Function to extract data
extrdata = function(k){
  ### Select first spatial point, create extent around it and crop the burned area data to the extent
  sub = s[k,]
  geom = sub$geometry ### Extract the geometry from the selected spatial data
  
  ### Transform the spatial data to a specific projection ("EPSG:4326" to "EPSG:3576")
  sub = vect(sub,geom=c("lon","lat"),"EPSG:4326")
  sub = project(sub, "EPSG:3576")
  
  ### Create an extent of 2 km around the transformed spatial data   
  extent = ext(ext(sub)[1]-2000,ext(sub)[2]+2000,ext(sub)[3]-2000,ext(sub)[4]+2000)
  
  ### Crop the burn area raster data to the extent
  cdb = crop(db, extent)
  cdb = as.points(cdb) ### Convert the cropped data to spatial points
  
  ### Find the nearest burn pixels with a 2 pixel search window
  near = nearby(sub, cdb, distance = sqrt((res(db)[1]*2.0001)**2+(res(db)[1]*2.0001)**2))
  ### Extract the coordinates of these nearest points and store them in a data frame
  nearpoints = cdb[c(near[,2])]
  ### Extract the coordinates of these nearest points and store them in a data frame
  near = data.table(crds(nearpoints))
  ### Remove temporary objects to free up memory
  rm(cdb,extent)
  
  ### Calculate the angle of the spatial point to the nearest burned points
  near$angle = mapply(bearing, a1 = crds(sub)[1], a2 = crds(sub)[2], b1 = near$x, b2 = near$y)
  sub$angle = circ.mean(near$angle) ## Calculate mean of angles
  
  ### Calculate the center of the boundary grid cell
  # Calculate the new center coordinates based on the angle and resolution
  centre = as.numeric(crds(sub))
  centre = centre + (res(db)[1]/2) * c(sin(sub$angle*pi/180), cos(sub$angle*pi/180))
  
  # Create an ellipse with specified parameters, minor axis = 1.5 times the spatial resolution, major axis = 3 times the spatial resolution  
  e = ellipse(a = res(db)[1]*1.5, b = res(db)[1]*3, centre = centre, phi = (360-sub$angle)*pi/180)
  # Make a spatial polygon of the ellipse
  p = vect(cbind(e$bdry[[1]][[1]], e$bdry[[1]][[2]]), "polygons", crs = "EPSG:3576")
  # Remove temporary objects to free up memory
  rm(centre,e,sub)
    
  ### Create a grid from the rasterized polygon with a 100 meter resolution
  grid = rast(ext(p), resolution = 100, crs = "EPSG:3576")
  grid = centroids(as.polygons(grid)) 
  grid = crop(grid, p)
  rm(p)
  
  ### Extract burned data to grid, indicating whether the grid point burned or not
  grid$burned = extract(db, grid)[,2]
  grid$burned = ifelse(is.na(grid$burned),0,1)
  
  # Additional checks and calculations
  # If the number of burned and unburned grid points is less than 6 (~10%), return no data (-9999) for all variables
  if(nrow(subset(grid, grid$burned == 1)) < 6) {return(c(geom,1, rep(-9999, 27)))}
  
  if(nrow(subset(grid, grid$burned == 0)) < 6) {return(c(geom,1, rep(-9999, 27)))}
  
  ### Extract Bartalev land cover to grid
  grid$lc = factor(extract(lc, grid)[,2])
  ### Extract burn history to grid, whether the grid point burned in the past 
  grid$bhistory = extract(bh, grid)[,2]
  
  #optional: plot grid and write grid to a shapefile to inspect
  #plot(grid, "burned")
  #plot(grid, "lc")

  #writeVector(sub, "~sub.shp", overwrite = T)
  #writeVector(grid, "~grid.shp", overwrite = T)
  #writeVector(cdb, "~cdb.shp", overwrite = T)
      
  ### Reproject grid to Lat/Lon to extract spatial data to grid points
  grid = project(grid, "EPSG:4326")
  grid$agb = extract(agb, grid)[,2] ## extract biomass
  grid$treecover = extract(treecoverr, grid)[,2]## extract tree cover
  grid$water = extract(gsw, grid)[,2]## extract surface water
  grid$road = extract(roads, grid)[,2]## extract data on roads
  grid$dem = extract(demr, grid)[,2]## extract elevation
  grid$pop = extract(population, grid)[,2]## extract population density
  
  grid = as.data.frame(grid) ## convert grid to data frame
  
  #### subset grid points with land cover class 26 which is unclassified
  grids = grid[grid$lc != 26, ]
  grids$lc = factor(grids$lc, levels = unique(grids$lc))
  # Calculate statistics for land cover in burned and unburned grid point
  a = table(grids$lc[grids$burned == 1]) #burned sites
  b = table(grids$lc[grids$burned == 0]) #unburned sites
  #perform chi square test on land cover
  maxlength = max(length(a),length(b))
  length(a) = maxlength
  length(b) = maxlength
  tlc = data.table(cbind(a,b))
  rownames(tlc) = names(a)
  tlc[is.na(tlc)] = 0
  tlc[,2] = round((tlc[,2]/sum(tlc[,2])*100),0)
  tlc[,1] = round((tlc[,1]/sum(tlc[,1])*100),0)
  lctest = try(chisq.test(tlc), silent = T)
  tlc$c = (tlc[,2]-tlc[,1])
  
  ## if the test returns an error for whatever reason, return no data for this test (-9999)
  if(class(lctest) == "try-error") {
  lac = -9999
  lacdif = -9999
  lacp = -9999} else{
  lac = ifelse(lctest$p.value < 0.05, 1,0)
  lacdif = as.character(paste(rownames(tlc)[which.min(tlc$c)], rownames(tlc)[which.max(tlc$c)],sep = "-"))
  lacp = lctest$p.value}  
  
  ### perform chi square test on surface water
  a= table(subset(grid, grid$burned == 1)$water)
  b= table(subset(grid, grid$burned == 0)$water)
  maxlength = max(length(a),length(b))
  length(a) = maxlength
  length(b) = maxlength
  tlc = data.table(cbind(a,b))
  tlc[is.na(tlc)] = 0
  tlc[,2] = round((tlc[,2]/sum(tlc[,2])*100),0)
  tlc[,1] = round((tlc[,1]/sum(tlc[,1])*100),0)
  lctest = try(chisq.test(tlc), silent = T)
  tlc$c = (tlc[,2]-tlc[,1])
  
  ## if the test returns an error for whatever reason, return no data for this test (-9999)
  if(class(lctest) == "try-error") {
  water = -9999
  waterb = -9999
  waterub = -9999
  waterdif = -9999
  waterp = -9999} else if(mean(subset(grid, grid$burned == 1)$water) > 0.25) {
  water = 1
  waterb = mean(subset(grid, grid$burned == 1)$water,na.rm = T)
  waterub = mean(subset(grid, grid$burned == 0)$water,na.rm = T)
  waterdif = waterub - waterb 
  waterp = 0.00001} else  {
  water = ifelse(lctest$p.value < 0.05 & mean(subset(grid, grid$burned == 1)$water,na.rm = T) < mean(subset(grid, grid$burned == 0)$water,na.rm = T), 1,0)
  waterb = mean(subset(grid, grid$burned == 1)$water,na.rm = T)
  waterub = mean(subset(grid, grid$burned == 0)$water,na.rm = T)
  waterdif = waterub - waterb 
  waterp = lctest$p.value}  
  
  ### perform chi square test on roads
  a= table(subset(grid, grid$burned == 1)$road)
  b= table(subset(grid, grid$burned == 0)$road)
  maxlength = max(length(a),length(b))
  length(a) = maxlength
  length(b) = maxlength
  tlc = data.table(cbind(a,b))
  tlc[is.na(tlc)] = 0
  tlc[,2] = round((tlc[,2]/sum(tlc[,2])*100),0)
  tlc[,1] = round((tlc[,1]/sum(tlc[,1])*100),0)
  lctest = try(chisq.test(tlc), silent = T)
  tlc = data.table(tlc)
  tlc$c = (tlc[,2]-tlc[,1])
  
  ## if the test returns an error for whatever reason, return no data for this test (-9999)
  if(class(lctest) == "try-error") {
  road = -9999
  roadb = -9999
  roadub = -9999
  roaddif = -9999
  roadp = -9999} else if (mean(subset(grid, grid$burned == 1)$road) > 0.25 ){
  road = 1
  roadb = mean(subset(grid, grid$burned == 1)$road,na.rm = T)
  roadub = mean(subset(grid, grid$burned == 0)$road,na.rm = T)
  roaddif = roadub - roadb
  roadp = 0.0001} else {
  road = ifelse(lctest$p.value < 0.05 & mean(subset(grid, grid$burned == 1)$road,na.rm = T) < mean(subset(grid, grid$burned == 0)$road,na.rm = T), 1,0)
  roadb = mean(subset(grid, grid$burned == 1)$road,na.rm = T)
  roadub = mean(subset(grid, grid$burned == 0)$road,na.rm = T)
  roaddif = roadub - roadb 
  roadp = lctest$p.value}  
  
  ### perform t-test on burn history values
  ttest = try(t.test(subset(grid, grid$burned == 1)$bhistory,subset(grid, grid$burned == 0)$bhistory), silent = T)
  if(class(ttest) == "try-error") {
  bhistory = -9999
  bhistorydif = -9999
  bhistoryp = -9999
  bhistoryb = -9999
  bhistoryub = -9999} else{
  bhistory = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) < as.numeric(ttest$estimate[2]), 1,0)
  bhistoryp = ttest$p.value
  bhistoryb = mean(subset(grid, grid$burned == 1)$bhistory,na.rm = T)
  bhistoryub = mean(subset(grid, grid$burned == 0)$bhistory,na.rm = T)
  bhistorydif = bhistoryub-bhistoryb}  
  
  ### perform t-test on elevation and terrain slope
  ttest = try(t.test(subset(grid, grid$burned == 1)$dem,subset(grid, grid$burned == 0)$dem), silent = T)
  if(class(ttest) == "try-error") {
  dem = -9999
  demb = -9999
  demub = -9999
  demdif = -9999
  demp = -9999
  slope = -9999} else{
  dem = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) > as.numeric(ttest$estimate[2]), 1,0)
  demb = mean(subset(grid, grid$burned == 1)$dem,na.rm = T)
  demub = mean(subset(grid, grid$burned == 0)$dem,na.rm = T)
  demdif = demub-demb
  demp = ttest$p.value
  slope = round(atan(demdif/300)*(180/pi),2)}    
  
  ### perform t-test on tree cover
  ttest = try(t.test(subset(grid, grid$burned == 1)$treecover,subset(grid, grid$burned == 0)$treecover), silent = T)
  if(class(ttest) == "try-error") {
  treecover = -9999
  treecoverb = -9999
  treecoverub = -9999
  treecoverdif = -9999
  treecoverp = -9999} else{
  treecover = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) > as.numeric(ttest$estimate[2]), 1,0)
  treecoverb = mean(subset(grid, grid$burned == 1)$treecover,na.rm = T)
  treecoverub = mean(subset(grid, grid$burned == 0)$treecover,na.rm = T)
  treecoverdif = treecoverub-treecoverb
  treecoverp = ttest$p.value}  
  
  ### perform t-test on above-ground biomass
  ttest = try(t.test(subset(grid, grid$burned == 1)$agb,subset(grid, grid$burned == 0)$agb), silent = T)
  if(class(ttest) == "try-error") {
  agbt = -9999
  agbb = -9999
  agbub = -9999
  agbdif = -9999
  agbp = -9999} else{
  agbt = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) > as.numeric(ttest$estimate[2]), 1,0)
  agbb = mean(subset(grid, grid$burned == 1)$agb,na.rm = T)
  agbub = mean(subset(grid, grid$burned == 0)$agb,na.rm = T)
  agbdif = agbub-agbb
  agbp = ttest$p.value}  
  
  ### perform t-test on population density
  ttest = try(t.test(subset(grid, grid$burned == 1)$pop,subset(grid, grid$burned == 0)$pop), silent = T)
  if(class(ttest) == "try-error") {
  pop = -9999
  popdif = -9999
  popp = -9999} else{
  pop = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) < as.numeric(ttest$estimate[2]), 1,0)
  popdif = round(as.numeric(ttest$estimate[2])-as.numeric(ttest$estimate[1]),2)
  popp = ttest$p.value}
  
  dat = data.frame(geom = geom, landscapefail = 0,lac = lac,lacdif = lacdif,lacp = round(lacp,5),treecover= treecover,treecoverb =round(treecoverb,3),treecoverub = round(treecoverub,3),treecoverdif = round(treecoverdif,3),treecoverp = round(treecoverp,5),agb = agbt,agbb = round(agbb,3),agbub = round(agbub,3), agbdif = round(agbdif,3),agbp = round(agbp,5),bhistory =  bhistory,bhistoryb = round(bhistoryb,4),bhistoryub= round(bhistoryub,4),bhistorydif= round(bhistorydif,3),bhistoryp =round(bhistoryp,5),water = water,waterb = round(waterb,3),waterub = round(waterub,3),waterdif= round(waterdif,3),waterp = round(waterp,5),road = road,roadb = round(roadb,3),roadub = round(roadub,3),roaddif= round(roaddif,3),roadp = round(roadp,5),dem = dem,demb = round(demb,3),demub= round(demub,3),demdif= round(demdif,3),demp=round(demp,5),slope= slope,pop= pop, popdif = round(popdif,3), popp = round(popp,5))
  
  ### return all values, with rounded decimals
  return(dat)
  gc()}