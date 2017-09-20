rhoPen <- function(x, zeta1, zeta2) {
  ifelse(x <= -zeta1, 0.5*(x+zeta1)^2, ifelse(x >= zeta2, 0.5*(x-zeta2)^2, 0))
}
# Look: knot spacing m_q