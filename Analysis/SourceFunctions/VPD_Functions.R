funvpd = function(k){
    sub = s[k,]
    
    # Determine date range
    days = 5
    dates = seq.Date(from = sub$burndate - (days-1), to = sub$burndate + days, by = "day") ### Select n number of days before and after last burn date
    months = unique(format(dates,"%m")) ### Select number of unique months in this range
    
    # Load VPD data
    sfilesvpd = filesvpd[filesvpd$year == year & filesvpd$month %in% months,]
    vpd = rast(sfilesvpd$files)
    evpd = data.table(vpd = as.numeric(extract(vpd,vect(sub,geom=c("lon","lat"),"EPSG:4326")))[2:(nlyr(vpd)+1)])
    evpd$datetime = time(vpd)
    rm(vpd)
                
    evpd$date = as.Date(evpd$datetime)
    subs = evpd[evpd$date %in% dates,]
    rm(evpd,sfilesvpd)
        
    so = data.table("ID" = sub$ID, "datetime" = subs$datetime, "burndate" = sub$burndate, "vpd" = round(subs$vpd,3))
            
    subs = subs[order(subs$datetime),]
    rownames(subs) = seq(1,nrow(subs),1)
    
    if(is.na(mean(subs$vpd))) {return(list(so,"Missing VPD data", 9999))}
    if(nrow(subs) < 20) {return(list(so,"Missing VPD data - end of December",9999))}
    
    ttest <- try(t.test(subs$vpd[73:120], subs$vpd[121:168]), silent = TRUE)
      if (inherits(ttest, "try-error")) {
      return(list(so, "Error in t-test: Data contains missing or invalid values", 9999))}
    wetter = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) > as.numeric(ttest$estimate[2]), 1,0)
    if(wetter == 0) {return(list(so,"No change in VPD", ttest$p.value))}
    
    ts = ts(subs$vpd, start = -5, frequency = 24)
    is = ts.intersect(y = ts, y1 = stats::lag(ts, k = -24))
    sc = efp(y ~ y1, data = is, type = "Score-CUSUM")
    
    bp = try(breakpoints(y~y1, data = is), silent = T)
    
    #plot(ts, col = "lightgray", lwd = 2)
    #lines(fitted(bp))
    #lines(confint(bp))
    rm(subs,ts,is)
    
    if(is.character(bp) & wetter == 1) {return(list(so,"Lower VPD but no breakpoint", ttest$p.value))}
    if(is.character(bp) & wetter == 0) {return(list(so,"No change in VPD", ttest$p.value))}
    if(is.na(bp$breakpoints)[1] & wetter == 1) {return(list(so,"Lower VPD but no breakpoint", ttest$p.value))}
    if(is.na(bp$breakpoints)[1] & wetter == 0) {return(list(so,"No change in VPD", ttest$p.value))}
        
    ### Select the breakpoint closest to day 0
    breaks = tryCatch({
     suppressWarnings(confint(bp)[[1]][which.min(abs(confint(bp)[[1]][, 2] - 96)),])}, error = function(e) {})
     
    rm(bp)
    
    if(is.na(breaks[1])& wetter == 1 ) {return(list(so,"Lower VPD but no breakpoint", ttest$p.value))}
    if(is.na(breaks[1])& wetter == 0 ) {return(list(so,"No change in VPD", ttest$p.value))}
    
    breaks = c(breaks[1]:breaks[length(breaks)])-96 ### Calculate confidence range of the break point
    brkvpd = ifelse(min(breaks)/24 >= -1 | max(breaks) <= 2, 1, 0) ## If the mean uncertainty of the breaks falls on last burn day up to two days later, count breakpoint
    
    if(brkvpd == 0 & wetter == 1 ) {return(list(so,"Lower VPD but no breakpoint", ttest$p.value))}
    if(brkvpd == 0 & wetter == 0 ) {return(list(so,"No change in VPD", ttest$p.value))}
    if(brkvpd == 1 & wetter == 1 ) {return(list(so,"Lower VPD and breakpoint", ttest$p.value))}
    if(brkvpd == 1 & wetter == 0 ) {return(list(so,"No change in VPD", ttest$p.value))}}