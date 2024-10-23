wfsi = function(winddir, firedir, windsp) {ifelse(sqrt((winddir-firedir)**2) > pi, windsp*(0.5-sin(((2*pi)-sqrt((winddir-firedir)**2))/2)), windsp*(0.5-sin((sqrt((winddir-firedir)**2))/2)))}


funws = function(k){
    sub = s[k,]
    
    # Determine date range
    days = 5 #5 days of burning and 5 days without burning
    dates = seq.Date(from = as.Date(sub$burndate) - (days-1), to = as.Date(sub$burndate) + days, by = "day") 
    ### Select n number of days before and after last burn date
    months = unique(format(dates,"%m")) ### Select number of unique months in this range
        
    # Load wind speed data    
    sfilesws = filesws[filesws$year == year & filesws$month %in% months,]
    ws = rast(sfilesws$files)
    ews = data.table(ws= as.numeric(extract(ws,vect(sub,geom=c("lon","lat"),"EPSG:4326")))[2:(nlyr(ws)+1)])
    ews$datetime = time(ws)
    rm(ws)
    
    # Load wind direction data
    sfileswd = fileswd[fileswd$year == year & fileswd$month %in% months,]
    wd = rast(sfileswd$files)
    ewd = data.table(wd= as.numeric(extract(wd,vect(sub,geom=c("lon","lat"),"EPSG:4326")))[2:(nlyr(wd)+1)])
    ewd$datetime = time(wd)
    rm(wd)    
      
    ews = merge(ews,ewd, by = "datetime", all = T)
    ews = ews[order(ews$datetime),]
    rownames(ews) = seq(1,nrow(ews),1)
    ews$date = as.Date(ews$datetime)
    subs = ews[ews$date %in% dates,]
    rm(ews,ewd,sfileswd,sfilesws)
    
    
    if(is.na(mean(subs$wd))) {return(list(data.table("ID" = sub$ID, "datetime" = subs$datetime, "burndate" = sub$burndate, "wfsi" = rep(NA)),"Missing Wind data",9999))}
    if(nrow(subs) < 20) {return(list(data.table("ID" = sub$ID, "datetime" = subs$datetime, "burndate" = sub$burndate, "wfsi" = rep(NA)),"Missing Wind data - end of December", 9999))}
    
    #subs$wfsi = wfsi(winddir = subs$wd*pi/180, firedir = sub$firedirection*pi/180, windsp = subs$ws)
    subs$wfsi = ifelse(sqrt((subs$wd * pi / 180 - sub$firedirection * pi / 180)^2) > pi,
                      subs$ws * (0.5 - sin((2 * pi - sqrt((subs$wd * pi / 180 - sub$firedirection * pi / 180)^2)) / 2)),
                      subs$ws * (0.5 - sin(sqrt((subs$wd * pi / 180 - sub$firedirection * pi / 180)^2) / 2)))
                      
    so = data.table("ID" = sub$ID, "datetime" = subs$datetime, "burndate" = sub$burndate, "wfsi" = round(subs$wfsi,3))
    
            
    ttest <- try(t.test(subs$wfsi[73:120], subs$wfsi[121:168]), silent = TRUE)
      if (inherits(ttest, "try-error")) {
      return(list(so, "Error in t-test: Data contains missing or invalid values", 9999))}
    windstop = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) > as.numeric(ttest$estimate[2]), 1,0)
    
    if(windstop == 0) {return(list(so,"No change in wfsi", ttest$p.value))}
    
    ts = ts(subs$wfsi, start = -5, frequency = 24)
    is = ts.intersect(y = ts, y1 = stats::lag(ts, k = -24))
    sc = efp(y ~ y1, data = is, type = "Score-CUSUM")
    
    bp = try(breakpoints(y~y1, data = is), silent = T)
    
    #plot(ts, col = "lightgray", lwd = 2)
    #lines(fitted(bp))
    #lines(confint(bp))
    rm(subs,ts,is)
        
    if(is.character(bp) & windstop == 1) {return(list(so,"Lower wfsi but no breakpoint", ttest$p.value))}
    if(is.character(bp) & windstop == 0) {return(list(so,"No change in wfsi", ttest$p.value))}
    if(is.na(bp$breakpoints)[1] & windstop == 1) {return(list(so,"Lower wfsi but no breakpoint", ttest$p.value))}
    if(is.na(bp$breakpoints)[1] & windstop == 0) {return(list(so,"No change in wfsi", ttest$p.value))}
    
    ### Select the breakpoint closest to day 0
    breaks = tryCatch({
     suppressWarnings(confint(bp)[[1]][which.min(abs(confint(bp)[[1]][, 2] - 96)),])}, error = function(e) {})
        
    rm(bp)
    
    if(is.na(breaks[1])& windstop == 1 ) {return(list(so,"Lower wfsi but no breakpoint", ttest$p.value))}
    if(is.na(breaks[1])& windstop == 0 ) {return(list(so,"No change in wfsi", ttest$p.value))}
    
    breaks = c(breaks[1]:breaks[length(breaks)])-96 ### Calculate confidence range of the break point
    brkwfsi = ifelse(min(breaks/24) >= -1 | max(breaks/24) <= 2, 1, 0) ## If the mean uncertainty of the breaks falls on last burn day up to two days later, count breakpoint
    
    if(brkwfsi == 0 & windstop == 1 ) {return(list(so,"Lower wfsi but no breakpoint", ttest$p.value))}
    if(brkwfsi == 0 & windstop == 0 ) {return(list(so,"No change in wfsi", ttest$p.value))}
    if(brkwfsi == 1 & windstop == 1 ) {return(list(so,"Lower wfsi and breakpoint", ttest$p.value))}
    if(brkwfsi == 1 & windstop == 0 ) {return(list(so,"No change in wfsi", ttest$p.value))}}
