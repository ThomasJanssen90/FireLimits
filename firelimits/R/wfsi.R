wfsi = function(winddir, firedir, windsp) {
  ifelse(sqrt((winddir-firedir)**2) > pi, windsp*(0.5-sin(((2*pi)-sqrt((winddir-firedir)**2))/2)),
         windsp*(0.5-sin((sqrt((winddir-firedir)**2))/2)))}