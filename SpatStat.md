# SpatStat

This demo will explore point process features in the `spatstat` package.
Note the spatstat “package” actually contains a set of R packages.

For additional `spatstat` details run the following code
`vignette('getstart')` in your R console.

#### 1. Generate Point Process Data

Note this will require both a $n$ and $\boldsymbol{s}$

``` r
lambda <- 100
n <- rpois(1, lambda)

x <- rbeta(n, 1, 1)
y <- rbeta(n, 1, 1)

sim_pp <- ppp(x, y, window = owin(xrange = c(0,1), yrange = c(0,1)))
plot(sim_pp, main = 'simulated huckleberry map')
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-1-1.png)

There are also a collection of methods in spatstat to simulation random
point process, see `rpoispp` or many others. These can more closely link
to our assumed generative models when considering, for example trend
surfaces or the inclusion of raster data.

``` r
CSR <- rpoispp(100)
plot(CSR)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-2-1.png)

**Activity 1:** Generate a point process, where the intensity is a
non-uniform function of x and y.

``` r
pp <- rpoispp(lambda = function(x,y) { exp(4 * x + 3 *y)})
plot(pp)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-3-1.png)

#### 2. Summarize Point Process

``` r
summary(sim_pp)
```

    Planar point pattern:  93 points
    Average intensity 93 points per square unit

    Coordinates are given to 16 decimal places

    Window: rectangle = [0, 1] x [0, 1] units
    Window area = 1 square unit

``` r
summary(pp)
```

    Planar point pattern:  84 points
    Average intensity 84 points per square unit

    Coordinates are given to 16 decimal places

    Window: rectangle = [0, 1] x [0, 1] units
    Window area = 1 square unit

#### 3. Assess Spatial Structure - Does it exhibit CSP?

Remember, the $F$ function corresponds to empty space and $G$
corresponds to nearest neighbor space. The $K$ function explores the
number of points in a specified radius.

``` r
plot(envelope(sim_pp,Fest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-6-1.png)

``` r
plot(envelope(sim_pp,Gest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-6-2.png)

``` r
plot(envelope(sim_pp,Kest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-6-3.png)

``` r
plot(envelope(CSR,Fest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-6-4.png)

``` r
plot(envelope(CSR,Gest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-6-5.png)

``` r
plot(envelope(CSR,Kest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-6-6.png)

``` r
plot(envelope(pp,Fest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-7-1.png)

``` r
plot(envelope(pp,Gest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-7-2.png)

``` r
plot(envelope(pp,Kest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-7-3.png)

*Q:* Think about the interpretation when observed lines are above and/or
below the gray curves.

#### 4. How does surface differ from CSR? Smoothed density surfaces

``` r
plot(density(sim_pp))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-8-1.png)

``` r
plot(density(CSR))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-8-2.png)

``` r
plot(density(pp))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-9-1.png)

### Model Fitting

The `ppm` function can be used for model fitting with a point process.
Consider an example that uses the location of *Beilschmiedia pendula*
rainforest trees.

``` r
plot(bei)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-10-1.png)

The `bei` dataset contains locations of trees in a tropical rain forest.
The point pattern is clearly non-homogenous

``` r
plot(envelope(bei, Kest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
plot(envelope(bei, Fest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-11-2.png)

``` r
plot(envelope(bei, Gest, verbose = F))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-11-3.png)

``` r
plot(density.ppp(bei))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-11-4.png)

The pattern in the intensity of the trees may be related to elevation
and the elevation gradient.

``` r
elev <- bei.extra$elev
grad <- bei.extra$grad
plot(elev)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
plot(grad)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-12-2.png)



The `ppm` function allows model fitting

``` r
tree.model <- ppm(bei ~ elev + grad);
tree.model
```

    Nonstationary Poisson process
    Fitted to point pattern dataset 'bei'

    Log intensity:  ~elev + grad

    Fitted trend coefficients:
    (Intercept)        elev        grad 
    -8.56355220  0.02143995  5.84646680 

                   Estimate        S.E.     CI95.lo     CI95.hi Ztest       Zval
    (Intercept) -8.56355220 0.341113849 -9.23212306 -7.89498134   *** -25.104675
    elev         0.02143995 0.002287866  0.01695581  0.02592408   ***   9.371155
    grad         5.84646680 0.255781018  5.34514522  6.34778838   ***  22.857313

``` r
plot(tree.model)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-13-1.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-13-2.png)

For more complicated models, `kppm` can be used for clustering behavior.

## Now lets fit our simulated data

``` r
sim_model <- ppm(pp ~ x + y);
sim_model
```

    Nonstationary Poisson process
    Fitted to point pattern dataset 'pp'

    Log intensity:  ~x + y

    Fitted trend coefficients:
    (Intercept)           x           y 
      0.5032604   3.5962761   2.7030856 

                 Estimate      S.E.    CI95.lo  CI95.hi Ztest     Zval
    (Intercept) 0.5032604 0.4983305 -0.4734493 1.479970       1.009893
    x           3.5962761 0.4966681  2.6228245 4.569728   *** 7.240803
    y           2.7030856 0.4458140  1.8293063 3.576865   *** 6.063259

``` r
plot(sim_model)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-14-1.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-14-2.png)

### Marked Point Patterns

Marked point process data contains meta data for each point. Rather than
just $\boldsymbol{s}$, we have $(\boldsymbol{s},m)$.

The marked information can either be categorical (multi-type) or
continuous.

The `lansing` data set contains locations of six types of trees.

``` r
plot(lansing, cols = 1:6)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-15-1.png)

``` r
plot(split(lansing))
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-16-1.png)



To analyze this data, consider the following model.

``` r
lansing.model <- ppm(lansing ~  marks - 1)
lansing.model
```

    Stationary multitype Poisson process
    Fitted to point pattern dataset 'lansing'

    Possible marks: 'blackoak', 'hickory', 'maple', 'misc', 'redoak' and 'whiteoak'

    Log intensity:  ~marks - 1

    Intensities:
    beta_blackoak  beta_hickory    beta_maple     beta_misc   beta_redoak 
              135           703           514           105           346 
    beta_whiteoak 
              448 

                  Estimate       S.E.  CI95.lo  CI95.hi Ztest      Zval
    marksblackoak 4.905275 0.08606630 4.736588 5.073962   ***  56.99414
    markshickory  6.555357 0.03771571 6.481435 6.629278   *** 173.80970
    marksmaple    6.242223 0.04410811 6.155773 6.328674   *** 141.52099
    marksmisc     4.653960 0.09759001 4.462687 4.845233   ***  47.68890
    marksredoak   5.846439 0.05376033 5.741070 5.951807   *** 108.75005
    markswhiteoak 6.104793 0.04724556 6.012194 6.197393   *** 129.21412

In contrast with this model, we can also include

``` r
lansing.model2 <- ppm(lansing ~  marks * polynom(x,y,3))
#lansing.model2
plot(lansing.model2)
```

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-1.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-2.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-3.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-4.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-5.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-6.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-7.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-8.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-9.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-10.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-11.png)

![](SpatStat_files/figure-commonmark/unnamed-chunk-18-12.png)

Similarly continuous marked data can be included as a predictor i the
`ppm` framework, potentially with interactions with spatially referenced
data.

Marked point process data can also be used for spatial-temporal point
patterns, where the year corresponds to the mark.



### More advanced point pattern models

##### Cluster processes

Clustering is not well defined. In general the idea is that the point
distances are shorter than expected. However, there “is a fundamental
ambiguity between heterogeneity and clustering” (Diggle 2007).

**Neyman-Scott Process**: This is a two stage process.

. Generate parents . For each parent, generate a set of offspring

**The shot noise processes** are variations on the Neyman-Scott process,
also with a two stage process.

**Strauss Process**: contains a term that allows repulsion by adjusting
the intensity in a vicinity of an existing point. The “hardcore” process
will make the intensity 0 for any pair of points less than a specified
distance $d_0$.
