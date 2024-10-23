
funsm = function(k){
    sub = s[k,]
    
    # Determine date range
    days = 5
    dates = seq.Date(from = sub$burndate - (days-1), to = sub$burndate + days, by = "day") ### Select n number of days before and after last burn date
    months = unique(format(dates,"%m")) ### Select number of unique months in this range
    
    # Load sm data
    sfilessm = filessm[filessm$year == year & filessm$month %in% months,]
    sm = rast(sfilessm$files)
    esm = data.table(sm = as.numeric(extract(sm,vect(sub,geom=c("lon","lat"),"EPSG:4326")))[2:(nlyr(sm)+1)])
    esm$datetime = time(sm)
    rm(sm)
            
    esm$date = as.Date(esm$datetime)
    subs = esm[esm$date %in% dates,]
    rm(esm,sfilessm)  
        
    so = data.table("ID" = sub$ID, "datetime" = subs$datetime, "burndate" = sub$burndate, "sm" = round(subs$sm,3))   
    
    subs = subs[order(subs$datetime),]
    rownames(subs) = seq(1,nrow(subs),1)
    
    if(is.na(mean(subs$sm))) {return(list(so,"Missing sm data", 9999))}
    if(nrow(subs) < 20) {return(list(so,"Missing sm data - end of December", 9999))}
    if(sqrt((mean(subs$sm[73:120])-mean(subs$sm[121:168]))**2) < 0.001) {return(list(so,"No change in sm", 9999))}
    
    ttest <- try(t.test(subs$sm[73:120], subs$sm[121:168]), silent = TRUE)
      if (inherits(ttest, "try-error")) {
      return(list(so, "Error in t-test: Data contains missing or invalid values", 9999))}
     
    wetter = ifelse(ttest$p.value < 0.05 & as.numeric(ttest$estimate[1]) < as.numeric(ttest$estimate[2]), 1,0)
    if(wetter == 0) {return(list(so,"No change in sm", ttest$p.value))}
    
    ts = ts(subs$sm, start = -5, frequency = 24)
    is = ts.intersect(y = ts, y1 = stats::lag(ts, k = -24))
    sc = efp(y ~ y1, data = is, type = "Score-CUSUM")
   
    bp = try(breakpoints(y~y1, data = is), silent = T)
    
    #plot(ts, col = "lightgray", lwd = 2)
    #lines(fitted(bp))
    #lines(confint(bp))
    rm(subs,ts,is)
    
    if(is.character(bp) & wetter == 1) {return(list(so,"Higher sm but no breakpoint", ttest$p.value))}
    if(is.character(bp) & wetter == 0) {return(list(so,"No change in sm", ttest$p.value))}
    if(is.na(bp$breakpoints)[1] & wetter == 1) {return(list(so,"Higher sm but no breakpoint", ttest$p.value))}
    if(is.na(bp$breakpoints)[1] & wetter == 0) {return(list(so,"No change in sm", ttest$p.value))}
        
    breaks = tryCatch({
     suppressWarnings(confint(bp)[[1]][which.min(abs(confint(bp)[[1]][, 2] - 96)),])}, error = function(e) {})
  
    rm(bp)
    
    if(length(breaks) == 0 & wetter == 1 ) {return(list(so,"Higher sm but no breakpoint", ttest$p.value))}
    if(length(breaks) == 0 & wetter == 0 ) {return(list(so,"No change in sm", ttest$p.value))}
    if(is.na(breaks[1])& wetter == 1 ) {return(list(so,"Higher sm but no breakpoint", ttest$p.value))}
    if(is.na(breaks[1])& wetter == 0 ) {return(list(so,"No change in sm", ttest$p.value))}
    
    breaks = c(breaks[1]:breaks[length(breaks)])-96 ### Calculate confidence range of the break point
    brksm = ifelse(min(breaks)/24 >= -1 | max(breaks)/24 <= 2, 1, 0) ## If the mean uncertainty of the breaks falls on last burn day up to two days later, count breakpoint
    
    if(brksm == 0 & wetter == 1 ) {return(list(so,"Higher sm but no breakpoint", ttest$p.value))}
    if(brksm == 0 & wetter == 0 ) {return(list(so,"No change in sm", ttest$p.value))}
    if(brksm == 1 & wetter == 1 ) {return(list(so,"Higher sm and breakpoint", ttest$p.value))}
    if(brksm == 1 & wetter == 0 ) {return(list(so,"No change in sm", ttest$p.value))}}