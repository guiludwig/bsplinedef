# bsplinedef
## Spatial deformation using b-splines
Computes coordinate transformations of the form 
    $(y1, y2) = (f1(x1, x2), f2(x1, x2))$ for spatial regression, 
    where a spatial process $Y$ on $(y1, y2)$ has known stationary 
    covariance function. The functions $f1$ and $f2$ are obtained
    via the tensor product of B-splines, with a regularization
    penalty to ensure $f1$, $f2$ are injective functions. The case 
    for $Y$ Gaussian with general covariance function is implemented, 
    as well as documentation for extensions to different spatial 
    covariance functions.
