insideHull <- function(p, HULL){
  # Checks whether a point is inside a convex hull HULL.
  # Really bad implementation.
  # step 1: polar coordinates:
  HULL <- unique(HULL) # in case last coordinate is repeated
  HULL <- HULL - colMeans(HULL)
  p <- p - colMeans(HULL)
  # https://en.wikipedia.org/wiki/Polar_coordinate_system#Converting_between_polar_and_Cartesian_coordinates
  tHULL <- cbind(sqrt(HULL[,1]^2 + HULL[,2]^2),
                 atan2(HULL[,2], HULL[,1]))
  tp <- cbind(sqrt(p[,1]^2 + p[,2]^2),
              atan2(p[,2], p[,1]))
  # Sort by angle
  # Find distance between point and edge
}
