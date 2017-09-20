injectSuff <- function(theta1, theta2, z1, z2){
  
  PEN <- sum(rhoPen(apply(theta1, 1, diff), z1, z2)) + # Differences by row
    sum(rhoPen(apply(theta1, 2, diff), z1, z2)) + # Differences by column
    sum(rhoPen(apply(theta2, 1, diff), z1, z2)) + # Differences by row
    sum(rhoPen(apply(theta2, 2, diff), z1, z2)) # Differences by column

  return(PEN)

}