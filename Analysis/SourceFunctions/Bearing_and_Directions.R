### Function to calculate bearing in degrees (0 and 360 is North) between two points
bearing = function(a1, a2, b1,b2) {
  twopi = pi*2
  RAD2DEG = 57.2957795130823209
  theta = atan2(b1 - a1, b2 - a2)
  if(theta < 0.0) {theta = theta+twopi}
  return(RAD2DEG * theta)}

### Function to average directions in degrees
circ.mean = function(angle){
  rad.m    = angle*(pi/180)
  mean.cos = mean(cos(rad.m))
  mean.sin = mean(sin(rad.m))
  x.deg    = atan2(mean.sin,mean.cos)*(180/pi)
  if(x.deg < 0) {x.deg=x.deg+360}
  return(x.deg)}